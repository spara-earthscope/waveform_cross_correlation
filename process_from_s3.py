#!/usr/bin/env python3
"""
Parallel preprocessing of MiniSEED files in S3 using ONLY memory (no temp files):
- list objects in s3://ssa2026/waveforms/
- for each *.mseed:
    * download into memory (bytes)
    * ObsPy read from BytesIO
    * detrend("demean"), detrend("linear"), taper(...)
    * write processed MiniSEED to a local output directory

Usage:
  python preprocess_s3_mseed_parallel_mem.py OUTPUT_DIR

Optional env vars:
  AWS_REGION      e.g. us-west-2 (optional)
  N_WORKERS       default: cpu_count-1
  TAPER_PCT       default: 0.05
  TAPER_TYPE      default: cosine
  OVERWRITE       "1" to overwrite existing outputs (default: 0)
  MAX_MB          skip files larger than this many MB (default: 200)

Requirements:
  pip install obspy boto3
"""

from __future__ import annotations

import io
import os
import sys
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import boto3
from botocore.config import Config
from botocore.exceptions import ClientError

from obspy import read


# --------------------------
# S3 settings
# --------------------------
S3_BUCKET = "ssa2026"
S3_PREFIX = "waveforms"  # normalized to "waveforms/" internally


# --------------------------
# Processing defaults
# --------------------------
DEFAULT_TAPER_PCT = float(os.getenv("TAPER_PCT", "0.05"))
DEFAULT_TAPER_TYPE = os.getenv("TAPER_TYPE", "cosine")
DEFAULT_OVERWRITE = os.getenv("OVERWRITE", "0") == "1"
MAX_MB = float(os.getenv("MAX_MB", "200"))  # safety limit for all-in-memory downloads


def default_workers() -> int:
    n = os.cpu_count() or 2
    return max(1, n - 1)


N_WORKERS = int(os.getenv("N_WORKERS", str(default_workers())))


def make_s3_client():
    region = os.getenv("AWS_REGION")
    cfg = Config(region_name=region) if region else None
    return boto3.client("s3", config=cfg) if cfg else boto3.client("s3")


def list_mseed_objects(bucket: str, prefix: str) -> list[tuple[str, int]]:
    """
    List (.mseed) keys under prefix and include object sizes.
    Returns list of (key, size_bytes).
    """
    s3 = make_s3_client()
    paginator = s3.get_paginator("list_objects_v2")

    out: list[tuple[str, int]] = []
    for page in paginator.paginate(Bucket=bucket, Prefix=prefix):
        for obj in page.get("Contents", []):
            key = obj.get("Key")
            size = int(obj.get("Size", 0) or 0)
            if key and key.lower().endswith(".mseed"):
                out.append((key, size))
    return out


def key_to_output_path(key: str, prefix: str, output_dir: Path) -> Path:
    rel = key[len(prefix):] if key.startswith(prefix) else key
    return output_dir / rel


def process_one_s3_object_in_memory(
    bucket: str,
    key: str,
    size_bytes: int,
    prefix: str,
    output_dir_str: str,
    taper_pct: float,
    taper_type: str,
    overwrite: bool,
    max_mb: float,
) -> dict:
    """
    Worker function:
      - fetch object bytes into memory (get_object)
      - read MiniSEED from BytesIO
      - detrend + taper
      - write processed MiniSEED to local output path
    """
    output_dir = Path(output_dir_str)
    out_path = key_to_output_path(key, prefix, output_dir)

    if out_path.exists() and not overwrite:
        return {"key": key, "status": "skipped", "reason": "exists", "out": str(out_path)}

    if size_bytes > int(max_mb * 1024 * 1024):
        return {
            "key": key,
            "status": "skipped",
            "reason": f"size {size_bytes/1024/1024:.1f} MB exceeds MAX_MB={max_mb}",
        }

    s3 = make_s3_client()

    try:
        # Download ALL bytes into memory
        resp = s3.get_object(Bucket=bucket, Key=key)
        data = resp["Body"].read()  # bytes

        # Read with ObsPy from memory
        bio = io.BytesIO(data)
        st = read(bio)

        # Process
        st.detrend("demean")
        st.detrend("linear")
        st.taper(max_percentage=taper_pct, type=taper_type)

        out_path.parent.mkdir(parents=True, exist_ok=True)
        st.write(str(out_path), format="MSEED")

        return {"key": key, "status": "ok", "traces": len(st), "out": str(out_path)}

    except ClientError as e:
        return {"key": key, "status": "error", "error": f"S3 error: {e}", "traceback": traceback.format_exc(limit=5)}
    except Exception as e:
        return {"key": key, "status": "error", "error": str(e), "traceback": traceback.format_exc(limit=5)}


def main():
    if len(sys.argv) != 2:
        print("Usage: python preprocess_s3_mseed_parallel_mem.py OUTPUT_DIR", file=sys.stderr)
        sys.exit(2)

    output_dir = Path(sys.argv[1]).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    prefix = S3_PREFIX
    if prefix and not prefix.endswith("/"):
        prefix += "/"

    print(f"S3 bucket : {S3_BUCKET}")
    print(f"S3 prefix : {prefix}")
    print(f"Output dir: {output_dir}")
    print(f"Workers   : {N_WORKERS}")
    print(f"Detrend   : demean + linear")
    print(f"Taper     : {DEFAULT_TAPER_PCT*100:.1f}% {DEFAULT_TAPER_TYPE}")
    print(f"Overwrite : {DEFAULT_OVERWRITE}")
    print(f"MAX_MB    : {MAX_MB} (skip larger objects)")

    objs = list_mseed_objects(S3_BUCKET, prefix)
    if not objs:
        print(f"No .mseed objects found under s3://{S3_BUCKET}/{prefix}")
        sys.exit(0)

    print(f"Found {len(objs)} .mseed object(s). Starting parallel processing...")

    ok = skipped = err = 0

    with ProcessPoolExecutor(max_workers=N_WORKERS) as ex:
        futures = [
            ex.submit(
                process_one_s3_object_in_memory,
                S3_BUCKET,
                key,
                size_bytes,
                prefix,
                str(output_dir),
                DEFAULT_TAPER_PCT,
                DEFAULT_TAPER_TYPE,
                DEFAULT_OVERWRITE,
                MAX_MB,
            )
            for key, size_bytes in objs
        ]

        for fut in as_completed(futures):
            res = fut.result()
            status = res.get("status")
            key = res.get("key")

            if status == "ok":
                ok += 1
                print(f"OK    {key} -> {res.get('out')} (traces={res.get('traces')})")
            elif status == "skipped":
                skipped += 1
                print(f"SKIP  {key} ({res.get('reason')})")
            else:
                err += 1
                print(f"ERR   {key}: {res.get('error')}")
                tb = res.get("traceback")
                if tb:
                    print(tb)

    print(f"\nDone. ok={ok} skipped={skipped} error={err}")


if __name__ == "__main__":
    main()


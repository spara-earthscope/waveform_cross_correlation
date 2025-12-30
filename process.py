#!/usr/bin/env python3
"""
Parallel preprocessing of MiniSEED files in a directory using ObsPy:
- detrend (demean)
- detrend (linear)
- taper

Writes processed files to an output directory, preserving filenames.

Usage:
  python preprocess_mseed_parallel.py input_dir output_dir

Optional env vars (or edit defaults below):
  N_WORKERS        number of parallel workers (default: cpu_count-1)
  TAPER_PCT        taper max_percentage (default: 0.05)
  TAPER_TYPE       taper type (default: "cosine")
  OVERWRITE        "1" to overwrite outputs (default: 0)

Requires:
  pip install obspy
"""

from __future__ import annotations

import os
import sys
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

from obspy import read


# --------------------------
# Defaults (can override via env vars)
# --------------------------
DEFAULT_TAPER_PCT = float(os.getenv("TAPER_PCT", "0.05"))
DEFAULT_TAPER_TYPE = os.getenv("TAPER_TYPE", "cosine")
DEFAULT_OVERWRITE = os.getenv("OVERWRITE", "0") == "1"

def default_workers() -> int:
    n = os.cpu_count() or 2
    return max(1, n - 1)

N_WORKERS = int(os.getenv("N_WORKERS", str(default_workers())))


def process_one(in_path: str, out_path: str, taper_pct: float, taper_type: str, overwrite: bool) -> dict:
    """
    Read one MiniSEED file, detrend (demean + linear), taper, write output.
    Returns a small status dict for logging.
    """
    inp = Path(in_path)
    outp = Path(out_path)

    if outp.exists() and not overwrite:
        return {"file": inp.name, "status": "skipped", "reason": "exists"}

    try:
        st = read(str(inp))

        # Detrend: remove mean then linear trend
        st.detrend("demean")
        st.detrend("linear")

        # Taper edges (per trace)
        st.taper(max_percentage=taper_pct, type=taper_type)

        outp.parent.mkdir(parents=True, exist_ok=True)
        st.write(str(outp), format="MSEED")

        return {"file": inp.name, "status": "ok", "traces": len(st)}
    except Exception as e:
        return {
            "file": inp.name,
            "status": "error",
            "error": str(e),
            "traceback": traceback.format_exc(limit=5),
        }


def main():
    if len(sys.argv) != 3:
        print("Usage: python preprocess_mseed_parallel.py input_dir output_dir", file=sys.stderr)
        sys.exit(2)

    input_dir = Path(sys.argv[1]).expanduser().resolve()
    output_dir = Path(sys.argv[2]).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    if not input_dir.is_dir():
        raise ValueError(f"Input is not a directory: {input_dir}")

    files = sorted(input_dir.rglob("*.mseed"))
    if not files:
        print(f"No .mseed files found under {input_dir}")
        sys.exit(0)

    print(f"Found {len(files)} file(s). Using {N_WORKERS} worker(s).")
    print(f"Detrend: demean + linear | Taper: {DEFAULT_TAPER_PCT*100:.1f}% {DEFAULT_TAPER_TYPE} | Overwrite: {DEFAULT_OVERWRITE}")

    futures = []
    ok = skipped = err = 0

    with ProcessPoolExecutor(max_workers=N_WORKERS) as ex:
        for in_path in files:
            # Preserve directory structure relative to input_dir
            rel = in_path.relative_to(input_dir)
            out_path = output_dir / rel
            futures.append(
                ex.submit(
                    process_one,
                    str(in_path),
                    str(out_path),
                    DEFAULT_TAPER_PCT,
                    DEFAULT_TAPER_TYPE,
                    DEFAULT_OVERWRITE,
                )
            )

        for fut in as_completed(futures):
            res = fut.result()
            if res["status"] == "ok":
                ok += 1
                print(f"OK      {res['file']}  (traces={res.get('traces')})")
            elif res["status"] == "skipped":
                skipped += 1
                print(f"SKIP    {res['file']}  ({res.get('reason')})")
            else:
                err += 1
                print(f"ERROR   {res['file']}: {res.get('error')}")
                tb = res.get("traceback")
                if tb:
                    print(tb)

    print(f"\nDone. ok={ok} skipped={skipped} error={err}")


if __name__ == "__main__":
    main()

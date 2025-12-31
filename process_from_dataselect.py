#!/usr/bin/env python3
"""
EarthScope dataselect (queryauth) -> async bulk downloads -> waveform processing
(demean + linear detrend + taper) -> write MiniSEED, in parallel across files.

This combines:
- Read metadata from local MiniSEED (headonly)
- Use metadata to build FDSN bulk dataselect requests
- Download from EarthScope via *async I/O* (aiohttp) to /queryauth
- Cache OAuth token *across worker processes* (Manager dict + lock)
- Process waveforms and write output MiniSEED per input file

USAGE
  python earthscope_async_parallel_process.py INPUT_DIR OUTPUT_DIR

REQUIRES
  pip install obspy earthscope-sdk aiohttp

OPTIONAL ENV VARS
  N_WORKERS        processes (default cpu_count-1)
  ASYNC_CONCURRENCY  max concurrent HTTP requests per worker (default 6)
  CHUNK_TRACES     max bulk lines per request (default 2000)
  PAD_START_S      pad start time (seconds) (default 0)
  PAD_END_S        pad end time (seconds) (default 0)
  TAPER_PCT        taper max_percentage (default 0.05)
  TAPER_TYPE       taper type (default cosine)
  OVERWRITE        "1" to overwrite outputs (default 0)
  BASE_URL         default https://service.earthscope.org/fdsnws/dataselect/1
"""

from __future__ import annotations

import asyncio
import io
import os
import sys
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Iterable, List, Tuple

import aiohttp
from obspy import Stream, UTCDateTime, read
from obspy.clients.fdsn.header import FDSNNoDataException  # just for semantics
from earthscope_sdk import EarthScopeClient

import multiprocessing as mp


# --------------------------
# Config
# --------------------------
BASE_URL = os.getenv("BASE_URL", "https://service.earthscope.org/fdsnws/dataselect/1").rstrip("/")
QUERYAUTH_URL = f"{BASE_URL}/queryauth"

TAPER_PCT = float(os.getenv("TAPER_PCT", "0.05"))
TAPER_TYPE = os.getenv("TAPER_TYPE", "cosine")
OVERWRITE = os.getenv("OVERWRITE", "0") == "1"

CHUNK_TRACES = int(os.getenv("CHUNK_TRACES", "2000"))
PAD_START_S = float(os.getenv("PAD_START_S", "0"))
PAD_END_S = float(os.getenv("PAD_END_S", "0"))

ASYNC_CONCURRENCY = int(os.getenv("ASYNC_CONCURRENCY", "6"))

def default_workers() -> int:
    n = os.cpu_count() or 2
    return max(1, n - 1)

N_WORKERS = int(os.getenv("N_WORKERS", str(default_workers())))


# --------------------------
# Token cache (shared across processes)
# --------------------------
def refresh_token() -> str:
    """
    Refresh (if necessary) and return EarthScope OAuth access token using earthscope_sdk.
    """
    client = EarthScopeClient()
    client.ctx.auth_flow.refresh_if_necessary()
    return client.ctx.auth_flow.access_token


def ensure_token_cached(token_cache, token_lock) -> str:
    """
    Return a cached token, refreshing under lock if missing.
    """
    tok = token_cache.get("access_token")
    if tok:
        return tok
    with token_lock:
        tok = token_cache.get("access_token")
        if tok:
            return tok
        tok = refresh_token()
        token_cache["access_token"] = tok
        token_cache["refresh_count"] = int(token_cache.get("refresh_count", 0)) + 1
        return tok


def refresh_token_cached(token_cache, token_lock) -> str:
    """
    Force refresh under lock and update cache (used after 401).
    """
    with token_lock:
        tok = refresh_token()
        token_cache["access_token"] = tok
        token_cache["refresh_count"] = int(token_cache.get("refresh_count", 0)) + 1
        return tok


# --------------------------
# Metadata -> bulk request lines
# --------------------------
BulkTuple = Tuple[str, str, str, str, UTCDateTime, UTCDateTime]

def iter_mseed_files(input_dir: Path) -> List[Path]:
    return sorted(p for p in input_dir.rglob("*.mseed") if p.is_file())


def trace_to_bulk_tuple(tr) -> BulkTuple:
    net = (tr.stats.network or "").strip()
    sta = (tr.stats.station or "").strip()
    loc = (tr.stats.location or "").strip()  # '' is valid for blank loc in FDSN
    cha = (tr.stats.channel or "").strip()
    t1 = tr.stats.starttime - PAD_START_S
    t2 = tr.stats.endtime + PAD_END_S
    return (net, sta, loc, cha, t1, t2)


def bulk_tuple_to_line(bt: BulkTuple) -> str:
    net, sta, loc, cha, t1, t2 = bt
    # Many bulk examples use "--" for blank location; keep explicit for readability.
    loc_out = loc if loc else "--"
    return f"{net} {sta} {loc_out} {cha} {t1.isoformat()} {t2.isoformat()}"


def chunked(items: List[BulkTuple], n: int) -> Iterable[List[BulkTuple]]:
    for i in range(0, len(items), n):
        yield items[i:i + n]


def output_path_for_input(input_file: Path, input_root: Path, output_root: Path) -> Path:
    rel = input_file.relative_to(input_root)
    return output_root / rel.with_name(f"{rel.stem}_earthscope_processed.mseed")


# --------------------------
# Async dataselect client
# --------------------------
class Unauthorized(Exception):
    pass


async def fetch_mseed_bytes(
    session: aiohttp.ClientSession,
    url: str,
    token: str,
    bulk_lines: str,
) -> bytes:
    """
    POST bulk request to EarthScope queryauth endpoint and return MiniSEED bytes.
    """
    headers = {
        "Authorization": f"Bearer {token}",
        "Content-Type": "text/plain",
        "Accept": "application/vnd.fdsn.mseed",
    }
    async with session.post(url, data=bulk_lines.encode("utf-8"), headers=headers) as resp:
        if resp.status == 401:
            raise Unauthorized("401 Unauthorized")
        # 204 No Content is common when no data matches
        if resp.status in (204, 404):
            return b""
        if resp.status != 200:
            txt = await resp.text()
            raise RuntimeError(f"HTTP {resp.status}: {txt[:300]}")
        return await resp.read()


async def download_bulk_chunks_async(
    url: str,
    token_cache,
    token_lock,
    chunks: List[str],
    concurrency: int,
) -> List[bytes]:
    """
    Download all bulk chunks concurrently (bounded by concurrency).
    If any request gets 401, refresh token once and retry that chunk.
    """
    sem = asyncio.Semaphore(concurrency)

    async def one_chunk(chunk_text: str) -> bytes:
        async with sem:
            token = ensure_token_cached(token_cache, token_lock)
            async with aiohttp.ClientSession() as session:
                try:
                    return await fetch_mseed_bytes(session, url, token, chunk_text)
                except Unauthorized:
                    # refresh token and retry once
                    token = refresh_token_cached(token_cache, token_lock)
                    return await fetch_mseed_bytes(session, url, token, chunk_text)

    tasks = [asyncio.create_task(one_chunk(ct)) for ct in chunks]
    return await asyncio.gather(*tasks)


def bytes_to_stream(mseed_bytes: bytes) -> Stream:
    if not mseed_bytes:
        return Stream()
    bio = io.BytesIO(mseed_bytes)
    return read(bio)


# --------------------------
# Worker (process) function
# --------------------------
def process_one_file_earthscope_async(
    input_file_str: str,
    input_root_str: str,
    output_root_str: str,
    token_cache,
    token_lock,
) -> dict:
    """
    Per input MiniSEED file:
      - read headonly metadata
      - build bulk lines, chunked
      - async-download each chunk (concurrently)
      - parse bytes into ObsPy Stream, process, write output
    """
    input_file = Path(input_file_str)
    input_root = Path(input_root_str)
    output_root = Path(output_root_str)

    out_path = output_path_for_input(input_file, input_root, output_root)
    if out_path.exists() and not OVERWRITE:
        return {"file": input_file.name, "status": "skipped", "reason": "exists", "out": str(out_path)}

    try:
        st_meta = read(str(input_file), headonly=True)
        bulk_tuples = [trace_to_bulk_tuple(tr) for tr in st_meta]
        if not bulk_tuples:
            return {"file": input_file.name, "status": "skipped", "reason": "no_traces"}

        # Build chunk texts (bulk lines joined with newlines)
        chunk_texts: List[str] = []
        for sub in chunked(bulk_tuples, CHUNK_TRACES):
            lines = [bulk_tuple_to_line(bt) for bt in sub]
            chunk_texts.append("\n".join(lines) + "\n")

        # Async download within this worker
        mseed_bytes_list = asyncio.run(
            download_bulk_chunks_async(
                QUERYAUTH_URL,
                token_cache,
                token_lock,
                chunk_texts,
                ASYNC_CONCURRENCY,
            )
        )

        st_out = Stream()
        empty_chunks = 0
        for b in mseed_bytes_list:
            if not b:
                empty_chunks += 1
                continue
            st_out += bytes_to_stream(b)

        if len(st_out) == 0:
            # mimic the earlier behavior
            raise FDSNNoDataException(f"No data returned (empty_chunks={empty_chunks})")

        # Process waveforms
        st_out.detrend("demean")
        st_out.detrend("linear")
        st_out.taper(max_percentage=TAPER_PCT, type=TAPER_TYPE)

        out_path.parent.mkdir(parents=True, exist_ok=True)
        st_out.write(str(out_path), format="MSEED")

        return {
            "file": input_file.name,
            "status": "ok",
            "out": str(out_path),
            "traces": len(st_out),
            "chunks": len(chunk_texts),
            "empty_chunks": empty_chunks,
        }

    except FDSNNoDataException as e:
        return {"file": input_file.name, "status": "skipped", "reason": str(e)}
    except Exception as e:
        return {
            "file": input_file.name,
            "status": "error",
            "error": str(e),
            "traceback": traceback.format_exc(limit=8),
        }


# --------------------------
# Main
# --------------------------
def main():
    if len(sys.argv) != 3:
        print("Usage: python earthscope_async_parallel_process.py INPUT_DIR OUTPUT_DIR", file=sys.stderr)
        sys.exit(2)

    input_dir = Path(sys.argv[1]).expanduser().resolve()
    output_dir = Path(sys.argv[2]).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    if not input_dir.is_dir():
        raise ValueError(f"INPUT_DIR is not a directory: {input_dir}")

    files = iter_mseed_files(input_dir)
    if not files:
        print(f"No .mseed files found under {input_dir}")
        sys.exit(0)

    # Shared token cache across workers
    ctx = mp.get_context("spawn")
    manager = ctx.Manager()
    token_cache = manager.dict()
    token_lock = manager.RLock()

    # Prime token once in parent (optional but reduces initial stampede)
    ensure_token_cached(token_cache, token_lock)

    print(f"Input dir : {input_dir}")
    print(f"Output dir: {output_dir}")
    print(f"EarthScope: {QUERYAUTH_URL}")
    print(f"Workers   : {N_WORKERS} (processes)")
    print(f"Async I/O : per-worker concurrency={ASYNC_CONCURRENCY}")
    print(f"Chunking  : CHUNK_TRACES={CHUNK_TRACES}  PAD_START_S={PAD_START_S}  PAD_END_S={PAD_END_S}")
    print(f"Process   : detrend(demean) + detrend(linear) + taper({TAPER_PCT*100:.1f}% {TAPER_TYPE})")
    print(f"Overwrite : {OVERWRITE}")
    print(f"Files     : {len(files)}")

    ok = skipped = err = 0

    with ProcessPoolExecutor(max_workers=N_WORKERS, mp_context=ctx) as ex:
        futures = [
            ex.submit(
                process_one_file_earthscope_async,
                str(p),
                str(input_dir),
                str(output_dir),
                token_cache,
                token_lock,
            )
            for p in files
        ]

        for fut in as_completed(futures):
            res = fut.result()
            status = res.get("status")
            if status == "ok":
                ok += 1
                print(f"OK    {res['file']} -> {res['out']} (traces={res.get('traces')}, chunks={res.get('chunks')}, empty={res.get('empty_chunks')})")
            elif status == "skipped":
                skipped += 1
                print(f"SKIP  {res['file']} ({res.get('reason')})")
            else:
                err += 1
                print(f"ERR   {res['file']}: {res.get('error')}")
                tb = res.get("traceback")
                if tb:
                    print(tb)

    refresh_count = int(token_cache.get("refresh_count", 0))
    print(f"\nDone. ok={ok} skipped={skipped} error={err} | token_refreshes={refresh_count}")


if __name__ == "__main__":
    main()

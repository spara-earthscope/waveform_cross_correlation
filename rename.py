#!/usr/bin/env python3
"""
Inspect directories for 'waveforms.mseed', rename it using the directory name,
and write it to a specified output directory.

Usage:
  python collect_mseed.py input_parent_dir output_dir

Example:
  input_parent_dir/
    event_40831111/waveforms.mseed
    event_40831112/waveforms.mseed

  output_dir/
    event_40831111_waveforms.mseed
    event_40831112_waveforms.mseed
"""

import sys
import shutil
from pathlib import Path


TARGET_FILENAME = "waveforms.mseed"


def main():
    if len(sys.argv) != 3:
        print("Usage: python collect_mseed.py input_parent_dir output_dir")
        sys.exit(1)

    input_parent = Path(sys.argv[1]).expanduser().resolve()
    output_dir = Path(sys.argv[2]).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    if not input_parent.is_dir():
        raise ValueError(f"Input path is not a directory: {input_parent}")

    count = 0

    for subdir in sorted(p for p in input_parent.iterdir() if p.is_dir()):
        mseed_path = subdir / TARGET_FILENAME

        if not mseed_path.exists():
            continue

        out_name = f"{subdir.name}_{TARGET_FILENAME}"
        out_path = output_dir / out_name

        shutil.copy2(mseed_path, out_path)
        print(f"Copied: {mseed_path} → {out_path}")
        count += 1

    print(f"\nDone. Copied {count} file(s).")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Convert fixed-width station metadata file to CSV.

Usage:
    python stations_to_csv.py input.txt output.csv
"""

import csv
import sys
from pathlib import Path

# Fixed-width column specification (1-based start, length)
SPECS = [
    ("net",       1,  2),
    ("sta",       4,  5),
    ("cha",      10,  3),
    ("loc",      14,  2),
    ("staname",  17, 30),
    ("lat",      48,  9),
    ("lon",      58, 10),
    ("elev_m",   69,  6),
    ("ondate",   76, 10),
    ("offdate",  87, 10),
    ("edepth",   98,  7),
    ("realtime",106,  8),
]

FIELDNAMES = [name for name, _, _ in SPECS]


def slice_field(line, start, length):
    """Extract fixed-width substring (1-based indexing)."""
    return line[start - 1:start - 1 + length].strip()


def parse_line(line: str) -> dict | None:
    if not line.strip():
        return None

    net = line[0:2].strip()
    sta = line[3:8].strip()
    cha = line[9:12].strip()
    loc = line[13:15].strip()

    rest = line[15:].rstrip()

    # Try with realtime (8 pieces total after rsplit)
    parts = rest.rsplit(None, 7)
    if len(parts) == 8:
        staname, lat, lon, elev, ondate, offdate, edepth, realtime = parts
    else:
        # Try without realtime (7 pieces total)
        parts = rest.rsplit(None, 6)
        if len(parts) != 7:
            raise ValueError(f"Unexpected format: {line}")
        staname, lat, lon, elev, ondate, offdate, edepth = parts
        realtime = ""

    return {
        "net": net,
        "sta": sta,
        "cha": cha,
        "loc": loc,
        "staname": staname.strip(),
        "lat": float(lat),
        "lon": float(lon),
        "elev_m": int(elev),
        "ondate": ondate,
        "offdate": offdate,
        "edepth": int(edepth),
        "realtime": realtime,
    }

def main():
    if len(sys.argv) != 3:
        print("Usage: python stations_to_csv.py input.txt output.csv")
        sys.exit(1)

    infile = Path(sys.argv[1])
    outfile = Path(sys.argv[2])

    with infile.open("r", encoding="utf-8", errors="replace") as fin, \
         outfile.open("w", newline="", encoding="utf-8") as fout:

        writer = csv.DictWriter(fout, fieldnames=FIELDNAMES)
        writer.writeheader()

        for line in fin:
            # Skip header line
            if line.startswith("NET STA"):
                continue

            row = parse_line(line)
            if row:
                writer.writerow(row)

    print(f"Wrote CSV file: {outfile}")


if __name__ == "__main__":
    main()

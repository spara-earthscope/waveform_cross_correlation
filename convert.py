#!/usr/bin/env python3
"""
Convert a fixed-width earthquake catalog (per provided COL/LEN spec)
to a comma-delimited (CSV) file.

Usage:
  python fixedwidth_to_csv.py input.txt output.csv
"""

import csv
import sys
from pathlib import Path

# Fixed-width specification: (field_name, start_col_1based, length, type)
# Note: start columns are 1-based as in your table.
SPECS = [
    ("year",    1,  4, "str"),
    ("mon",     6,  2, "str"),
    ("day",     9,  2, "str"),
    ("hour",   13,  2, "str"),
    ("min",    16,  2, "str"),
    ("sec",    19,  5, "float"),
    ("lat_deg",25,  3, "int"),
    ("lat_min",29,  5, "float"),
    ("lon_deg",34,  4, "int"),
    ("lon_min",39,  5, "float"),
    ("quality",45,  1, "str"),
    ("magnitude",47,3, "float"),
    ("depth_km",54, 6, "float"),
    ("nph",    60,  3, "int"),
    ("rms",    67,  5, "float"),
    ("eventid",73,  8, "str"),
]

FIELDNAMES = [name for name, *_ in SPECS] + ["lat_dd", "lon_dd", "iso_utc"]


def _slice(line: str, start_1based: int, length: int) -> str:
    """Return substring using 1-based start and fixed length."""
    i0 = start_1based - 1
    return line[i0:i0 + length]


def _to_number(s: str, kind: str):
    s = s.strip()
    if s == "":
        return None
    try:
        if kind == "int":
            return int(s)
        if kind == "float":
            return float(s)
        return s
    except ValueError:
        return None


def parse_line(line: str) -> dict | None:
    # Skip blank/comment lines
    if not line.strip():
        return None

    row = {}
    for name, start, length, kind in SPECS:
        raw = _slice(line, start, length)
        val = _to_number(raw, kind)
        if kind == "str" and val is not None:
            val = str(val).strip()
        row[name] = val

    # Compute decimal degrees (minutes -> degrees)
    # Assumes degrees may be signed; minutes are positive.
    lat_deg = row["lat_deg"]
    lat_min = row["lat_min"]
    lon_deg = row["lon_deg"]
    lon_min = row["lon_min"]

    def dd(deg, minutes):
        if deg is None or minutes is None:
            return None
        sign = -1 if deg < 0 else 1
        return sign * (abs(deg) + (minutes / 60.0))

    row["lat_dd"] = dd(lat_deg, lat_min)
    row["lon_dd"] = dd(lon_deg, lon_min)

    # Build ISO UTC timestamp (no timezone offset applied; values are UTC per comment)
    y, mo, d = row["year"], row["mon"], row["day"]
    hh, mm, ss = row["hour"], row["min"], row["sec"]
    if all(v is not None for v in (y, mo, d, hh, mm, ss)):
        # sec can be like 12.34; format to 2 decimals if present
        row["iso_utc"] = f"{int(y):04d}-{int(mo):02d}-{int(d):02d}T{int(hh):02d}:{int(mm):02d}:{float(ss):05.2f}Z"
    else:
        row["iso_utc"] = None

    return row


def main():
    if len(sys.argv) != 3:
        print("Usage: python fixedwidth_to_csv.py input.txt output.csv", file=sys.stderr)
        sys.exit(2)

    in_path = Path(sys.argv[1])
    out_path = Path(sys.argv[2])

    with in_path.open("r", encoding="utf-8", errors="replace") as fin, \
         out_path.open("w", newline="", encoding="utf-8") as fout:

        writer = csv.DictWriter(fout, fieldnames=FIELDNAMES)
        writer.writeheader()

        for line_no, line in enumerate(fin, start=1):
            row = parse_line(line.rstrip("\n"))
            if row is None:
                continue
            writer.writerow(row)

    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()

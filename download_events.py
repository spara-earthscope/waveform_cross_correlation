#!/usr/bin/env python3
"""
Match events against the SCEDC catalog (FDSN event service) using (lat, lon, iso_utc)
from a CSV, then download waveforms and write MiniSEED files.

Usage:
  python csv_scedc_match_to_mseed.py events.csv output_dir

Input CSV (minimum columns):
  lat_dd, lon_dd, iso_utc
Optional columns:
  eventid, magnitude

Outputs per row:
  output_dir/event_<eventid_or_row>/matched_event.xml   (QuakeML)
  output_dir/event_<eventid_or_row>/stations.xml        (StationXML)
  output_dir/event_<eventid_or_row>/waveforms.mseed     (MiniSEED)
"""

from __future__ import annotations

import csv
import math
import sys
from pathlib import Path

from obspy import Stream, UTCDateTime
from obspy.clients.fdsn import Client
from obspy.clients.fdsn.header import FDSNNoDataException
from obspy.core.event import Catalog


# --------------------------
# Providers (SCEDC)
# --------------------------
EVENT_PROVIDER = "SCEDC"     # match the SCEDC event catalog :contentReference[oaicite:1]{index=1}
WAVEFORM_PROVIDER = "SCEDC"  # download waveforms from SCEDC (can change to "IRIS"/"EARTHSCOPE")

# --------------------------
# Event matching (time+space)
# --------------------------
TIME_WINDOW_S = 30.0      # search +/- seconds around iso_utc
RADIUS_KM = 30.0          # search radius around CSV lat/lon
MAX_EVENTS_PER_ROW = 100  # safety cap

# --------------------------
# Waveform download
# --------------------------
T_BEFORE_S = 60.0
T_AFTER_S = 300.0
CHANNEL = "HH?,BH?"       # adjust as needed (e.g., "EH?,HH?,BH?")
LOCATION = "*"            # "*" any; "--" for blank only
NETWORK = "*"             # e.g., "CI,CE"
STATION = "*"
STATION_RADIUS_KM = 200.0
MAX_STATIONS = 50


def km_to_degrees(km: float) -> float:
    return km / 111.19


def parse_iso_utc(s: str) -> UTCDateTime:
    # Handles "2025-01-01T01:54:51.66Z"
    return UTCDateTime(s.strip())


def haversine_km(lat1, lon1, lat2, lon2) -> float:
    r = 6371.0
    p1, p2 = math.radians(lat1), math.radians(lat2)
    dphi = math.radians(lat2 - lat1)
    dl = math.radians(lon2 - lon1)
    a = math.sin(dphi / 2) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dl / 2) ** 2
    return 2 * r * math.asin(math.sqrt(a))


def read_events_csv(csv_path: Path) -> list[dict]:
    rows = []
    with csv_path.open("r", encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f)
        required = {"lat_dd", "lon_dd", "iso_utc"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"CSV missing required columns: {sorted(missing)}")

        for i, row in enumerate(reader, start=1):
            if not any((v or "").strip() for v in row.values()):
                continue
            evid = str(row.get("eventid", f"row{i}")).strip() or f"row{i}"
            rows.append(
                {
                    "row_index": i,
                    "eventid": evid,
                    "lat": float(row["lat_dd"]),
                    "lon": float(row["lon_dd"]),
                    "time": parse_iso_utc(row["iso_utc"]),
                    "magnitude": (row.get("magnitude") or "").strip(),
                }
            )
    return rows


def find_best_event_match(event_client: Client, lat: float, lon: float, t0: UTCDateTime, mag_str: str):
    """
    Query SCEDC events near (lat, lon, t0). Return (catalog, best_event, best_score)
    where best_score = (dt_seconds, dist_km).
    """
    t1 = t0 - TIME_WINDOW_S
    t2 = t0 + TIME_WINDOW_S
    maxradius_deg = km_to_degrees(RADIUS_KM)

    kwargs = dict(
        starttime=t1,
        endtime=t2,
        latitude=lat,
        longitude=lon,
        maxradius=maxradius_deg,
        limit=MAX_EVENTS_PER_ROW,
        includearrivals=False,
    )

    # If CSV magnitude exists, use as a loose hint to reduce unrelated matches.
    try:
        if mag_str and mag_str.lower() not in ("none", "nan"):
            mag = float(mag_str)
            kwargs["minmagnitude"] = max(mag - 0.2, -10.0)
    except ValueError:
        pass

    cat = event_client.get_events(**kwargs)

    best = None
    best_score = None  # (dt, dist_km)

    for ev in cat:
        origin = ev.preferred_origin() or (ev.origins[0] if ev.origins else None)
        if origin is None or origin.time is None or origin.latitude is None or origin.longitude is None:
            continue
        dt = abs(origin.time - t0)
        dist = haversine_km(lat, lon, origin.latitude, origin.longitude)
        score = (dt, dist)
        if best_score is None or score < best_score:
            best_score = score
            best = ev

    return cat, best, best_score


def main():
    if len(sys.argv) != 3:
        print("Usage: python csv_scedc_match_to_mseed.py events.csv output_dir", file=sys.stderr)
        sys.exit(2)

    csv_path = Path(sys.argv[1]).expanduser().resolve()
    out_dir = Path(sys.argv[2]).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    event_client = Client(EVENT_PROVIDER)
    wf_client = Client(WAVEFORM_PROVIDER)

    rows = read_events_csv(csv_path)
    if not rows:
        print("No events found in CSV.", file=sys.stderr)
        sys.exit(1)

    station_radius_deg = km_to_degrees(STATION_RADIUS_KM)

    for r in rows:
        evid = r["eventid"]
        lat, lon, t0 = r["lat"], r["lon"], r["time"]

        print(f"\nRow {r['row_index']}: eventid={evid}")
        print(f"  CSV: {t0.isoformat()}Z  lat={lat:.5f} lon={lon:.5f}")

        # 1) Match SCEDC catalog event near CSV row
        try:
            _cat, best, best_score = find_best_event_match(event_client, lat, lon, t0, r["magnitude"])
        except FDSNNoDataException:
            print("  SCEDC: no events returned.")
            continue

        if best is None:
            print("  SCEDC: returned events, but none had a usable origin. Skipping.")
            continue

        origin = best.preferred_origin() or best.origins[0]
        ot = origin.time
        dt, dist = best_score if best_score is not None else (None, None)
        print(f"  Matched SCEDC origin: {ot.isoformat()}Z  dt={dt:.2f}s  dist={dist:.2f} km")

        # 2) Per-event output folder
        ev_dir = out_dir / f"event_{evid}"
        ev_dir.mkdir(parents=True, exist_ok=True)

        # Save matched catalog object
        Catalog(events=[best]).write(str(ev_dir / "matched_event.xml"), format="QUAKEML")

        # 3) Find stations near matched event and download waveforms
        t1 = ot - T_BEFORE_S
        t2 = ot + T_AFTER_S

        try:
            inv = wf_client.get_stations(
                network=NETWORK,
                station=STATION,
                location=LOCATION,
                channel=CHANNEL,
                latitude=float(origin.latitude),
                longitude=float(origin.longitude),
                maxradius=station_radius_deg,
                starttime=t1,
                endtime=t2,
                level="channel",
            )
        except FDSNNoDataException:
            print("  No stations found for waveform request.")
            continue

        inv.write(str(ev_dir / "stations.xml"), format="STATIONXML")

        stations = sorted({(net.code, sta.code) for net in inv for sta in net})
        if len(stations) > MAX_STATIONS:
            print(f"  Found {len(stations)} stations; capping to {MAX_STATIONS}.")
            stations = stations[:MAX_STATIONS]
        else:
            print(f"  Found {len(stations)} stations.")

        st = Stream()
        for net_code, sta_code in stations:
            try:
                st += wf_client.get_waveforms(
                    network=net_code,
                    station=sta_code,
                    location=LOCATION,
                    channel=CHANNEL,
                    starttime=t1,
                    endtime=t2,
                    attach_response=False,
                )
            except FDSNNoDataException:
                continue
            except Exception as e:
                print(f"  Warning: waveform fetch failed for {net_code}.{sta_code}: {e}")

        if len(st) == 0:
            print("  No waveform traces returned.")
            continue

        # Optional: merge overlaps/gaps
        try:
            st.merge(method=1, fill_value="interpolate")
        except Exception:
            pass

        mseed_path = ev_dir / "waveforms.mseed"
        st.write(str(mseed_path), format="MSEED")
        print(f"  Wrote MiniSEED: {mseed_path}  (traces={len(st)})")

    print("\nDone.")


if __name__ == "__main__":
    main()

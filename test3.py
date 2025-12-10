# ============================================================
# 0) HV matched filtering parameters (edit these)
# ============================================================
from obspy import UTCDateTime

# Region & scan window (2018 Kīlauea eruption period is a good example)
CENTER_LAT, CENTER_LON = 19.3, -155.5
STARTTIME = "2018-05-01"   # scan start (inclusive)
ENDTIME   = "2018-05-03"   # scan end (exclusive) - try a couple of days first

# FDSN search for stations (vertical components around the volcano)
MAXRADIUS_DEG = 2.5        # tighten if you want only summit/East Rift Zone
N_STATIONS     = 6         # how many stations to include (top N that work)

# LP band typical for volcanoes (tune as needed)
FREQMIN = 0.5   # Hz
FREQMAX = 2.0   # Hz

# Template definition (YOU MUST SET THIS to a known LP event time)
# Example placeholder — set to the origin/peak time of a representative LP.
TEMPLATE_START   = "2018-05-02T00:15:10"  # <-- CHANGE to a known LP time
TEMPLATE_LENGTHS = 20.0                   # seconds of waveform to include per station

# Processing & detection knobs
TARGET_FS         = 20.0     # resample rate for templates & data (Hz)
TIME_NORMALIZATION = "onebit"  # 'onebit', 'ram', or None
RAM_WINDOW_SEC     = 5
APPLY_WHITENING    = True

# Detection thresholds
NET_SIGMA = 8.0    # network score threshold in units of MAD above median
STA_SIGMA = 5.0    # per-station NCC threshold in units of MAD
COINC_MIN = 3      # require at least this many stations above their own threshold
PEAK_SEP  = 0.5    # minimum separation between detections, fraction of template length

# Response removal (to velocity); keep broad pre-filter for deconvolution stability
PRE_FILT    = (0.005, 0.006, 45.0, 50.0)
WATER_LEVEL = 60

# ============================================================
# 1) Imports & small utilities
# ============================================================
import numpy as np
from obspy.clients.fdsn import Client
from obspy.geodetics import gps2dist_azimuth
from obspy.signal.util import next_pow_2
from scipy.signal import fftconvolve, find_peaks
import matplotlib.pyplot as plt
from tqdm import tqdm

client = Client("IRIS")  # EarthScope Data Services (FDSN)
EPS = 1e-12

# ============================================================
# 2) Station discovery for HV (vertical channels)
# ============================================================
def find_hv_stations(center_lat, center_lon, start, end, maxradius_deg=2.5,
                     n_stations=6, preferred=("HHZ", "BHZ")):
    """
    Return a list of up to n_stations entries: (net, sta, loc, cha).
    """
    t0, t1 = UTCDateTime(start), UTCDateTime(end)
    inv = client.get_stations(network="HV", starttime=t0, endtime=t1,
                              latitude=center_lat, longitude=center_lon,
                              maxradius=maxradius_deg,
                              channel="HHZ,BHZ", level="channel")
    cand = []
    for net in inv:
        for sta in net:
            # pick preferred channel (HHZ then BHZ)
            chans = sorted([ch for ch in sta.channels if ch.code in preferred],
                           key=lambda ch: preferred.index(ch.code))
            if not chans:
                continue
            ch = chans[0]
            loc = ch.location_code or ""
            cand.append((net.code, sta.code, loc, ch.code, sta.latitude, sta.longitude))
    if not cand:
        raise RuntimeError("No HV stations found—expand radius or adjust dates.")
    # Sort by distance to center and keep top n_stations
    cand.sort(key=lambda r: gps2dist_azimuth(center_lat, center_lon, r[4], r[5])[0])
    picked = [(c[0], c[1], c[2], c[3]) for c in cand[:n_stations]]
    return picked

# ============================================================
# 3) Fetch & preprocessing
# ============================================================
def fetch_stream(net, sta, loc, cha, t0, t1, attach_resp=True):
    st = client.get_waveforms(net, sta, loc, cha, t0, t1, attach_response=attach_resp)
    st.merge(fill_value="interpolate")
    return st

def remove_instrument(st, pre_filt, water_level, output="VEL"):
    st = st.copy()
    for tr in st:
        tr.remove_response(pre_filt=pre_filt, output=output, water_level=water_level)
    return st

def basic_preprocess(st, freqmin, freqmax, target_fs):
    st = st.copy()
    st.detrend("linear")
    st.detrend("demean")
    st.taper(max(10.0 / st[0].stats.sampling_rate, 0.05))
    st.filter("bandpass", freqmin=freqmin, freqmax=freqmax, corners=4, zerophase=True)
    if abs(st[0].stats.sampling_rate - target_fs) > 1e-6:
        st.resample(target_fs)
    st[0].data = st[0].data.astype(np.float64)
    return st

# ============================================================
# 4) Amplitude stabilization and spectral whitening
# ============================================================
def time_normalize(x, fs, mode="onebit", ram_win_sec=5):
    y = x.astype(np.float64).copy()
    if mode is None:
        return y
    if mode.lower() == "onebit":
        return np.sign(y)
    if mode.lower() == "ram":
        n = max(1, int(ram_win_sec * fs))
        kern = np.ones(n, dtype=np.float64) / n
        denom = np.convolve(np.abs(y), kern, mode="same")
        denom[denom == 0] = 1.0
        return y / denom
    raise ValueError("Unknown time normalization mode")

def spectral_whiten(x, fs, fmin, fmax, smooth_bins=5):
    n = len(x)
    nfft = next_pow_2(2 * n)
    X = np.fft.rfft(x, n=nfft)
    freqs = np.fft.rfftfreq(nfft, d=1.0/fs)
    amp = np.abs(X)
    if smooth_bins > 1:
        w = np.ones(smooth_bins) / smooth_bins
        amp = np.convolve(amp, w, mode="same")
    band = (freqs >= fmin) & (freqs <= fmax)
    Xw = np.zeros_like(X)
    Xw[band] = X[band] / (amp[band] + EPS)
    yw = np.fft.irfft(Xw, n=nfft)[:n]
    return yw

# ============================================================
# 5) Build per-station templates
# ============================================================
def build_templates(stations, tpl_start, tpl_len_s, fmin, fmax, target_fs,
                    time_norm_mode, ram_win_sec, apply_whiten=True):
    """
    Returns dict station_key -> template vector (numpy array).
    station_key is (net, sta, loc, cha).
    """
    t0 = UTCDateTime(tpl_start)
    t1 = t0 + tpl_len_s
    templates = {}
    for (net, sta, loc, cha) in stations:
        try:
            st = fetch_stream(net, sta, loc, cha, t0, t1, attach_resp=True)
            st = remove_instrument(st, PRE_FILT, WATER_LEVEL, output="VEL")
            st = basic_preprocess(st, fmin, fmax, target_fs)
            x = st[0].data
            if time_norm_mode:
                x = time_normalize(x, target_fs, mode=time_norm_mode, ram_win_sec=ram_win_sec)
            if apply_whiten:
                x = spectral_whiten(x, target_fs, fmin, fmax)
            # final normalization to unit energy (L2)
            e = np.sqrt(np.sum(x**2)) + EPS
            templates[(net, sta, loc, cha)] = x / e
        except Exception as e:
            print(f"[TPL WARN] {net}.{sta}.{loc}.{cha}: {e}")
    if not templates:
        raise RuntimeError("No templates were built. Check TEMPLATE_START/LENGTH and station list.")
    # Force all templates to same length (truncate to min)
    L = min(len(v) for v in templates.values())
    for k in list(templates.keys()):
        templates[k] = templates[k][:L]
    return templates

# ============================================================
# 6) Fast normalized cross-correlation (NCC) via FFT
# ============================================================
def ncc_fft(data, template):
    """
    Normalized cross-correlation of 'data' with 'template' (time-reversed inside).
    Returns array of length len(data) - len(template) + 1
    """
    N = len(template)
    if len(data) < N:
        return np.array([])
    # numerator: correlation
    num = fftconvolve(data, template[::-1], mode="valid")
    # denominator: sqrt( sum(data^2 over sliding window) * sum(template^2) )
    denom_data = fftconvolve(data**2, np.ones(N, dtype=np.float64), mode="valid")
    denom = np.sqrt(np.maximum(denom_data, 0.0))  # template already unit energy
    denom = np.maximum(denom, EPS)
    return num / denom

# ============================================================
# 7) Scan one day for all stations & stack
# ============================================================
def scan_one_day(day_start, day_len_s, stations, templates, fmin, fmax, target_fs,
                 time_norm_mode, ram_win_sec, apply_whiten=True):
    """
    Returns (times_center, station_ncc_dict, stacked_score)
    times_center are UTCDateTime array for the center time of each NCC sample.
    """
    t0 = UTCDateTime(day_start)
    t1 = t0 + day_len_s
    station_ncc = {}
    # First pass: compute per-station NCC arrays
    for (net, sta, loc, cha) in stations:
        key = (net, sta, loc, cha)
        if key not in templates:
            continue
        tpl = templates[key]
        try:
            st = fetch_stream(net, sta, loc, cha, t0, t1, attach_resp=True)
            st = remove_instrument(st, PRE_FILT, WATER_LEVEL, output="VEL")
            st = basic_preprocess(st, fmin, fmax, target_fs)
            x = st[0].data
            if time_norm_mode:
                x = time_normalize(x, target_fs, mode=time_norm_mode, ram_win_sec=ram_win_sec)
            if apply_whiten:
                x = spectral_whiten(x, target_fs, fmin, fmax)
            ncc = ncc_fft(x, tpl)
            if ncc.size:
                station_ncc[key] = ncc
        except Exception as e:
            print(f"[SCAN WARN] {net}.{sta}.{loc}.{cha} ({t0}): {e}")
            continue

    if not station_ncc:
        return np.array([]), {}, np.array([])

    # Align lengths and build time vector (center time of template window)
    L = min(len(arr) for arr in station_ncc.values())
    for k in list(station_ncc.keys()):
        station_ncc[k] = station_ncc[k][:L]
    dt = 1.0 / target_fs
    # NCC index i represents the window x[i : i+N-1]; use center time of the window
    tpl_len = len(next(iter(templates.values())))
    center_offset = (tpl_len / 2.0) * dt
    times = np.array([t0 + (i * dt) + center_offset for i in range(L)], dtype=object)

    # Stack across stations (robust mean)
    M = np.vstack([v for v in station_ncc.values()])
    stacked = np.nanmean(M, axis=0)
    return times, station_ncc, stacked

# ============================================================
# 8) Robust thresholding & detection
# ============================================================
def mad(x):
    med = np.median(x)
    return med, 1.4826 * np.median(np.abs(x - med))

def detect_events(times, station_ncc, stacked, net_sigma=8.0, sta_sigma=5.0,
                  coinc_min=3, peak_sep=0.5, tpl_len_s=10.0, fs=20.0):
    if stacked.size == 0:
        return []
    # Network threshold
    med_s, mad_s = mad(stacked)
    net_thr = med_s + net_sigma * (mad_s if mad_s > 0 else 1.0)
    # Per-station thresholds
    sta_thr = {}
    for k, v in station_ncc.items():
        m, d = mad(v)
        sta_thr[k] = m + sta_sigma * (d if d > 0 else 1.0)

    # Coincidence count over stations at each sample
    keys = list(station_ncc.keys())
    M = np.vstack([station_ncc[k] for k in keys])
    coinc = np.sum(M >= np.array([sta_thr[k] for k in keys])[:, None], axis=0)

    # Candidate mask: both network score and coincidence condition
    mask = (stacked >= net_thr) & (coinc >= coinc_min)

    # Peak picking on network score with minimum separation
    min_dist = int(max(1, peak_sep * tpl_len_s * fs))
    idx, _ = find_peaks(stacked, distance=min_dist)
    idx = [i for i in idx if mask[i]]

    detections = []
    for i in idx:
        det_time = times[i]
        det_score = stacked[i]
        # simple metadata
        detections.append({
            "time": det_time,
            "score": float(det_score),
            "coinc": int(coinc[i]),
            "net_thr": float(net_thr)
        })
    return detections

# ============================================================
# 9) Plotting helpers
# ============================================================
def plot_stacked(times, stacked, net_thr, detections, title="HV network matched filter score"):
    if stacked.size == 0:
        print("[PLOT] Nothing to plot.")
        return
    t0 = times[0]
    x = np.array([(ti - t0) for ti in times], dtype=float)  # seconds from day start
    plt.figure(figsize=(12, 4))
    plt.plot(x, stacked, linewidth=1)
    plt.axhline(net_thr, linestyle="--")
    for det in detections:
        xs = (det["time"] - t0)
        plt.axvline(xs, alpha=0.5)
    plt.xlabel(f"Seconds since {t0}")
    plt.ylabel("Stacked NCC")
    plt.title(title)
    plt.tight_layout()
    plt.show()

def summary_print(dets, label="detections"):
    if not dets:
        print(f"[INFO] No {label}.")
        return
    print(f"[INFO] {len(dets)} {label}:")
    for d in dets[:20]:
        print(f"  - {d['time']}  score={d['score']:.2f}  coinc={d['coinc']}")
    if len(dets) > 20:
        print(f"  ... and {len(dets)-20} more")

# ============================================================
# 10) Run the HV matched filter
# ============================================================
DAY_SEC = 24 * 3600

def run_hv_matched_filter():
    # 1) Stations
    stations = find_hv_stations(CENTER_LAT, CENTER_LON, STARTTIME, ENDTIME,
                                maxradius_deg=MAXRADIUS_DEG, n_stations=N_STATIONS)
    print("[HV] Stations:", ", ".join([f"{n}.{s}.{l}.{c}" for (n,s,l,c) in stations]))

    # 2) Templates
    templates = build_templates(stations, TEMPLATE_START, TEMPLATE_LENGTHS,
                                FREQMIN, FREQMAX, TARGET_FS,
                                TIME_NORMALIZATION, RAM_WINDOW_SEC, apply_whiten=APPLY_WHITENING)
    tpl_len_s = len(next(iter(templates.values()))) / TARGET_FS
    print(f"[TPL] Built {len(templates)} templates; template length = {tpl_len_s:.2f} s")

    # 3) Scan over days
    tscan0, tscan1 = UTCDateTime(STARTTIME), UTCDateTime(ENDTIME)
    all_detections = []
    day = tscan0
    while day < tscan1:
        window = min(DAY_SEC, (tscan1 - day))
        print(f"[SCAN] {day} + {window/3600:.1f} h")
        times, sta_ncc, stacked = scan_one_day(day, window, stations, templates,
                                               FREQMIN, FREQMAX, TARGET_FS,
                                               TIME_NORMALIZATION, RAM_WINDOW_SEC,
                                               apply_whiten=APPLY_WHITENING)
        if stacked.size:
            # thresholds for this day
            med_s, mad_s = mad(stacked)
            net_thr = med_s + NET_SIGMA * (mad_s if mad_s > 0 else 1.0)
            dets = detect_events(times, sta_ncc, stacked, NET_SIGMA, STA_SIGMA,
                                 COINC_MIN, PEAK_SEP, tpl_len_s, TARGET_FS)
            plot_stacked(times, stacked, net_thr, dets,
                         title=f"HV matched filter score ({day.date})")
            summary_print(dets, label=f"detections on {day.date}")
            all_detections.extend(dets)
        else:
            print("[SCAN] No usable data that day.")
        day += DAY_SEC

    # 4) Summary
    print("\n=== SUMMARY ===")
    summary_print(all_detections, label="total detections")
    return {"stations": stations, "templates": templates, "detections": all_detections}

if __name__ == "__main__":
    _ = run_hv_matched_filter()


# =========================
# 0) HV: Hawaiian Volcano Observatory Network (HV) parameters
# =========================
from obspy import UTCDateTime

# Example  (Hawaii HV, 2018):
CENTER_LAT, CENTER_LON = 19.3, -155.5
STARTTIME = "2018-05-01"
ENDTIME   = "2018-08-01"

# Search radius (degrees) for HV stations around the center:
MAXRADIUS_DEG = 10.0

# Prefer inter-station distances in this range (km) for microseism band:
MIN_DIST_KM, MAX_DIST_KM = 100.0, 400.0

# Band of interest (volcano):
FREQMIN = 1   # Hz  (20 s)
FREQMAX = 5   # Hz  (5 s)

# Windowing and correlation
WINDOW_LENGTH = 24 * 3600     # seconds (one-day windows)
DEFAULT_MAX_LAG = 300         # seconds; will auto-bump based on pair distance
TIME_NORMALIZATION = "onebit" # "onebit", "ram", or None
RAM_WINDOW_SEC = 10

# Response removal (to velocity)
PRE_FILT = (0.005, 0.006, 45.0, 50.0)
WATER_LEVEL = 60

# =========================
# 1) Imports and utilities
# =========================
import numpy as np
from obspy.clients.fdsn import Client
from obspy import Stream
from obspy.signal.util import next_pow_2
from obspy.geodetics import gps2dist_azimuth
from scipy.signal import fftconvolve
import matplotlib.pyplot as plt
from tqdm import tqdm

# =========================
# 2) HV station discovery
# =========================
def find_ta_pair(center_lat, center_lon, start, end,
                 maxradius_deg=6.0, min_km=100.0, max_km=400.0,
                 preferred_channels=("BHZ", "HHZ")):
    """
    Query EarthScope FDSN for HV stations active in [start, end] within the search radius.
    Return two tuples: (net, sta, loc, cha, lat, lon) for a pair whose distance falls in [min_km, max_km].
    Preference is given to vertical BH/HH channels with location codes available.
    """
    client = Client("IRIS")  # EarthScope Data Services FDSN endpoint
    t0, t1 = UTCDateTime(start), UTCDateTime(end)

    # Get stations with vertical broadband channels; include channel level for loc/code
    inv = client.get_stations(network="HV",
                              starttime=t0, endtime=t1,
                              latitude=center_lat, longitude=center_lon,
                              maxradius=maxradius_deg,
                              channel="BHZ,HHZ",
                              level="channel")
    

    # Build candidate list with explicit (net, sta, loc, cha, lat, lon)
    cand = []
    for net in inv:
        for sta in net:
            sta_lat, sta_lon = sta.latitude, sta.longitude
            # Sort channels to prefer BHZ then HHZ, and valid location codes
            chans = sorted([ch for ch in sta.channels if ch.code in preferred_channels],
                           key=lambda ch: preferred_channels.index(ch.code))
            if not chans:
                continue
            ch = chans[0]
            loc = ch.location_code or ""
            cand.append((net.code, sta.code, loc, ch.code, sta_lat, sta_lon))

    if len(cand) < 2:
        raise RuntimeError("Not enough HV stations found—try moving the center, expanding MAXRADIUS_DEG, or changing dates.")

    # Evaluate all pairs; keep those within distance bounds; pick mid-range distance
    best = None
    target = 0.5 * (min_km + max_km)
    for i in range(len(cand)):
        for j in range(i + 1, len(cand)):
            a, b = cand[i], cand[j]
            dist_m, _, _ = gps2dist_azimuth(a[4], a[5], b[4], b[5])
            dist_km = dist_m / 1000.0
            if min_km <= dist_km <= max_km:
                score = abs(dist_km - target)
                if best is None or score < best[0]:
                    best = (score, dist_km, a, b)

    if best is None:
        # Fall back: choose the closest pair and warn
        min_pair = None
        for i in range(len(cand)):
            for j in range(i + 1, len(cand)):
                a, b = cand[i], cand[j]
                dist_m, _, _ = gps2dist_azimuth(a[4], a[5], b[4], b[5])
                dist_km = dist_m / 1000.0
                if min_pair is None or dist_km < min_pair[0]:
                    min_pair = (dist_km, a, b)
        dist_km, a, b = min_pair
        print(f"[INFO] No pair in [{min_km:.0f},{max_km:.0f}] km; using closest pair at {dist_km:.1f} km.")
        chosen = (dist_km, a, b)
    else:
        _, dist_km, a, b = best
        chosen = (dist_km, a, b)

    dist_km, A, B = chosen
    # Return: (net, sta, loc, cha), distance_km
    return (A[0], A[1], A[2], A[3]), (B[0], B[1], B[2], B[3]), dist_km

# =========================
# 3) Download + preprocessing
# =========================
def fetch_stream(net, sta, loc, cha, t0, t1, client_name="IRIS", attach_resp=True):
    client = Client(client_name)
    st = client.get_waveforms(net, sta, loc, cha, t0, t1, attach_response = attach_resp)
    st.merge(fill_value="interpolate")
    return st

def basic_preprocess(st, freqmin, freqmax):
    st = st.copy()
    st.detrend("linear")
    st.detrend("demean")
    st.taper(max(10.0 / st[0].stats.sampling_rate, 0.05))
    st.filter("bandpass", freqmin=freqmin, freqmax=freqmax, corners=4, zerophase=True)
    return st

def remove_instrument(st, pre_filt, water_level, output="VEL"):
    st = st.copy()
    for tr in st:
        tr.remove_response(pre_filt=pre_filt, output=output, water_level=water_level)
    return st

# =========================
# 4) Amplitude stabilization & spectral whitening
# =========================
def time_normalize(data, fs, mode="onebit", ram_win_sec=10):
    x = data.astype(np.float64).copy()
    if mode is None:
        return x
    if mode.lower() == "onebit":
        return np.sign(x)
    if mode.lower() == "ram":
        n = max(1, int(ram_win_sec * fs))
        kern = np.ones(n, dtype=np.float64) / n
        denom = np.convolve(np.abs(x), kern, mode="same")
        denom[denom == 0] = 1.0
        return x / denom
    raise ValueError("Unknown time normalization mode")

def spectral_whiten(x, fs, fmin, fmax, eps=1e-10, smooth_bins=5):
    n = len(x)
    nfft = next_pow_2(2 * n)
    X = np.fft.rfft(x, n=nfft)
    freqs = np.fft.rfftfreq(nfft, d=1.0 / fs)

    amp = np.abs(X)
    if smooth_bins > 1:
        w = np.ones(smooth_bins) / smooth_bins
        amp = np.convolve(amp, w, mode="same")

    band = (freqs >= fmin) & (freqs <= fmax)
    Xw = np.zeros_like(X)
    Xw[band] = X[band] / (amp[band] + eps)

    xw = np.fft.irfft(Xw, n=nfft)[:n]
    return xw

# =========================
# 5) Cross-correlation & trimming
# =========================
def xcorr_full(x, y):
    # r_xy[k] = sum_n x[n] * y[n + k]
    return fftconvolve(x, y[::-1], mode="full")

def trim_to_maxlag(ccf, fs, max_lag):
    n = len(ccf)
    zero_idx = (n - 1) // 2  # for equal-length inputs
    half = int(max_lag * fs)
    left = max(0, zero_idx - half)
    right = min(n, zero_idx + half + 1)
    ccf_trim = ccf[left:right]
    lags_sec = np.arange(left - zero_idx, right - zero_idx) / fs
    return ccf_trim, lags_sec

# =========================
# 6) Per-window pipeline
# =========================
def process_one_window(t0, window_length, pair, freqmin, freqmax, pre_filt, water_level,
                       time_norm_mode, ram_win_sec, max_lag):
    t1 = t0 + window_length
    (n1, s1, l1, c1), (n2, s2, l2, c2) = pair
    try:
        st1 = fetch_stream(n1, s1, l1, c1, t0, t1, "IRIS", True)
        st2 = fetch_stream(n2, s2, l2, c2, t0, t1, "IRIS", True)
    except Exception as e:
        print(f"[WARN] Fetch failed {t0}—{t1}: {e}")
        return None, None

    if len(st1) == 0 or len(st2) == 0:
        return None, None

    tr1, tr2 = st1[0], st2[0]
    fs1, fs2 = tr1.stats.sampling_rate, tr2.stats.sampling_rate
    if not np.isclose(fs1, fs2):
        fs = min(fs1, fs2)
        st1.resample(fs)
        st2.resample(fs)
        tr1, tr2 = st1[0], st2[0]
    fs = tr1.stats.sampling_rate

    st1 = remove_instrument(st1, PRE_FILT, WATER_LEVEL, output="VEL")
    st2 = remove_instrument(st2, PRE_FILT, WATER_LEVEL, output="VEL")
    st1 = basic_preprocess(st1, freqmin, freqmax)
    st2 = basic_preprocess(st2, freqmin, freqmax)

    x = st1[0].data.astype(np.float64)
    y = st2[0].data.astype(np.float64)

    if time_norm_mode:
        x = time_normalize(x, fs, mode=time_norm_mode, ram_win_sec=ram_win_sec)
        y = time_normalize(y, fs, mode=time_norm_mode, ram_win_sec=ram_win_sec)

    x = spectral_whiten(x, fs, freqmin, freqmax)
    y = spectral_whiten(y, fs, freqmin, freqmax)

    ccf_full = xcorr_full(x, y)
    ccf_trim, lags_sec = trim_to_maxlag(ccf_full, fs, max_lag)

    m = np.max(np.abs(ccf_trim)) or 1.0
    ccf_trim = ccf_trim / m
    return ccf_trim, lags_sec

# =========================
# 7) Loop & stack
# =========================
def daterange(starttime, endtime, step_seconds):
    t0 = UTCDateTime(starttime)
    t1 = UTCDateTime(endtime)
    while t0 + step_seconds <= t1:
        yield t0
        t0 += step_seconds

def stack_ccfs(ccf_list):
    arr = np.vstack(ccf_list)
    return np.nanmean(arr, axis=0)

# =========================
# 8) Plot helpers
# =========================
def plot_ccf(lags, ccf, title="Stacked CCF"):
    plt.figure(figsize=(10, 4))
    plt.plot(lags, ccf, linewidth=1.5)
    plt.axvline(0, linestyle="--")
    plt.xlabel("Lag (s)")
    plt.ylabel("Amplitude (arb.)")
    plt.title(title)
    plt.tight_layout()
    plt.show()

def plot_record_section(lags, ccfs, title="Per-window CCFs (record section)"):
    M = np.vstack(ccfs)
    extent = [lags[0], lags[-1], 0, M.shape[0]]
    plt.figure(figsize=(10, 6))
    plt.imshow(M, aspect="auto", extent=extent, origin="lower")
    plt.axvline(0, linestyle="--")
    plt.xlabel("Lag (s)")
    plt.ylabel("Window index (time order)")
    plt.title(title)
    plt.colorbar(label="CCF amplitude")
    plt.tight_layout()
    plt.show()

# =========================
# 9) Run the HV pipeline
# =========================
def run_ta_demo():
    # Pick a HV station pair automatically
    A, B, dist_km = find_ta_pair(CENTER_LAT, CENTER_LON, STARTTIME, ENDTIME,
                                 maxradius_deg=MAXRADIUS_DEG,
                                 min_km=MIN_DIST_KM, max_km=MAX_DIST_KM)
    print(f"[HV] Using stations: {A[0]}.{A[1]}.{A[2]}.{A[3]} ↔ {B[0]}.{B[1]}.{B[2]}.{B[3]}  (~{dist_km:.1f} km apart)")

    # Heuristic lag: enough to cover surface-wave travel times (slower ~2 km/s)
    max_lag = max(DEFAULT_MAX_LAG, int(dist_km / 2.0))  # seconds
    if max_lag > DEFAULT_MAX_LAG:
        print(f"[INFO] MAX_LAG set to {max_lag} s to cover ~{dist_km:.0f} km at ~2 km/s")

    ccfs = []
    lags = None
    for t0 in tqdm(list(daterange(STARTTIME, ENDTIME, WINDOW_LENGTH))):
        ccf, lags_sec = process_one_window(
            t0, WINDOW_LENGTH, (A, B), FREQMIN, FREQMAX, PRE_FILT, WATER_LEVEL,
            TIME_NORMALIZATION, RAM_WINDOW_SEC, max_lag
        )
        if ccf is not None:
            ccfs.append(ccf)
            if lags is None:
                lags = lags_sec

    if not ccfs:
        raise RuntimeError("No windows succeeded—try a different date range or expand the search radius.")

    ccf_stack = stack_ccfs(ccfs)
    plot_ccf(lags, ccf_stack, title=f"HV stacked CCF: {A[1]} ↔ {B[1]}  ({dist_km:.0f} km)")
    plot_record_section(lags, ccfs, title="HV per-window CCFs (each day)")
    return {"pair_A": A, "pair_B": B, "distance_km": dist_km, "lags": lags, "ccf_stack": ccf_stack, "ccfs": ccfs}

if __name__ == "__main__":
    _ = run_ta_demo()
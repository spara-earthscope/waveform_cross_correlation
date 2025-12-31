# Waveform Cross-Correlation: A Step-by-Step Tutorial (Undergrad Level)

## What you’ll learn

* The intuition and math behind cross-correlation
* How to build a practical processing workflow (preprocess → correlate → stack → interpret)
* What the causal vs. acausal sides mean and when results approximate an empirical Green’s function (EGF)
* Common pitfalls and example use cases (tomography, detection, nodal arrays)

---

## 1) Big picture & intuition

**Cross-correlation ≈ sliding dot product.** You slide one waveform past another and measure similarity at each lag; a strong peak means the shapes align at that time shift. In frequency space, correlation brings in a **complex conjugate**. 

**Important property:** Unlike convolution, **cross-correlation is not commutative**—swapping the order flips the result. 

**(Optional refresher)** Moving between time and frequency uses the Fourier transform; in the correlation derivation you substitute variables and obtain the frequency-domain form. 

---

## 2) When and why we use it

* **Matched filtering / detection:** Correlate a known template with noisy data; a spike indicates a detection.  
* **Ambient-noise interferometry (EGFs):** Correlating long records at two stations can approximate the Green’s function sensitive to structure between them; stacking boosts weak signals. 

**Limitation reminders:** Sources aren’t everywhere (often dominated by ocean microseisms), and some processing choices affect how directly your result relates to the true Green’s function. 

---

## 3) The standard processing workflow

> This is the backbone you can adapt to your dataset.

### Step 3.1 — Gather and segment data

* Acquire continuous waveform data for a station pair (or array).
* Cut into windows (e.g., **day-long windows** are common; shorter or longer can work depending on compute/storage and the periods of interest). 

### Step 3.2 — Preprocess each window

1. **Remove instrument response** (convert to ground motion).
2. **Detrend / de-mean.**
3. **Band-pass filter** to the periods you care about (e.g., microseism bands). 

### Step 3.3 — Stabilize amplitudes (to avoid loud transients dominating)

* **Time-domain normalization** (e.g., running mean/one-bit): suppress earthquakes and other bursts. 
* **Spectral whitening**: flatten the spectrum; conceptually similar to the time-domain step but in frequency. 

### Step 3.4 — Compute cross-correlation for each window

* Correlate the two preprocessed traces over a max lag that covers the travel times of interest. (Remember: **order matters** for correlation.) 

### Step 3.5 — Stack many windows

* **Stack (average)** the per-window correlations to enhance stable features and suppress incoherent noise. This is what reveals surface-wave energy in ambient-noise studies. 
* The final stacked correlation is your **CCF** (cross-correlation function).

### Step 3.6 — Interpret the CCF (EGF perspective)

* **Acausal (negative lag)** vs. **causal (positive lag)** sides correspond to “who is the virtual source vs. receiver.” Averaging sides is one way studies connect CCFs to EGFs. 
* Expect strong **Rayleigh-wave** energy (retrograde particle motion) in surface-wave bands. 

---

## 4) Minimal pseudocode (Python ecosystem)

Several Python tools can run this workflow end-to-end: **ObsPy**, **MSNoise**, **NoisePy**—each with pros/cons; cite them if you use them. 

```text
for each station pair (A,B):
  read continuous data for days D1..DN
  for each day Di:
    x, y = load(A, Di), load(B, Di)
    x = remove_instrument(x); y = remove_instrument(y)
    x = detrend(x); y = detrend(y)
    x = bandpass(x, fmin, fmax); y = bandpass(y, fmin, fmax)
    x = time_normalize(x); y = time_normalize(y)
    X = spectral_whiten(fft(x)); Y = spectral_whiten(fft(y))
    ccf_i = ifft( X * conj(Y) )   # correlation in frequency domain
    save(ccf_i)
  CCF = stack(all ccf_i)
  analyze(CCF)  # pick peaks, measure group/phase times, etc.
```

*(The conj() appears because correlation uses a complex conjugate in frequency space.)* 

---

## 5) Reading the result (what to look for)

* **Clear surface-wave moveout** with distance when you plot many station pairs: record sections often show coherent Rayleigh energy around your band (e.g., ~16 s).  
* **Retrograde motion** animation/diagnostic plots help verify you’re seeing Rayleigh waves. 

---

## 6) Practical tips & pitfalls

* **Window length:** Long enough to capture sources traversing the array but balanced against compute/storage. 
* **Source distribution:** Don’t expect isotropy—microseisms often dominate; this affects how “EGF-like” your CCF is. 
* **Processing choices matter:** Different normalization/whitening choices can change how directly your CCF approximates a Green’s function; be explicit and compare variants. 
* **Order awareness:** Correlation is not commutative; keep your station order consistent in filenames/metadata. 

---

## 7) Short hands-on lab (30–60 min)

1. **Pick a band** (e.g., 0.05–0.1 Hz for ocean microseisms).
2. **Grab one day of data** for two nearby broadband stations.
3. **Run Steps 3.2–3.4** for that day; plot the per-day CCF.
4. **Repeat for 7–30 days**; **stack** CCFs and re-plot; identify the first Rayleigh-wave arrival. 
5. **Flip order** of inputs and confirm the acausal/causal flip. 

**Stretch goal:** Use picked lags vs. distance to estimate Rayleigh phase velocity (intro eikonal tomography idea). 

---

## 8) Showcase applications (for inspiration)

* **Crustal/sedimentary tomography** with USArray and nodal networks (e.g., Long Beach).  
* **Event detection & subspace methods** via matched filtering. 
* **Relocations:** Differential times from waveform correlation tighten earthquake sequences (e.g., Magna, Utah). 
* **Volcano & hydrothermal systems:** Imaging geyser plumbing and magma pressurization with cross-correlation. 
* **Structural health:** Building damage detection using correlation changes. 

---

## 9) Quick self-check

* Can you explain why correlation is **not** commutative? (What flips?) 
* What do **causal** and **acausal** sides represent physically? 
* Why do we **stack** many windows? 
* Name two **normalization** steps and why they’re used. 

---

## 10) Where to go next

* Try a turnkey package (**MSNoise** or **NoisePy**) for automated pipelines; start with **ObsPy** to learn the pieces. Remember to **cite the software** you use. 

---

*This tutorial is distilled from the provided transcript/lecture on waveform cross-correlation and ambient-noise applications.* 

```
# --- USER SETTINGS (edit me) ---
NETWORK1, STATION1, LOCATION1, CHANNEL1 = "IU", "ANMO", "00", "BHZ"   # Albuquerque, NM (example)
NETWORK2, STATION2, LOCATION2, CHANNEL2 = "IU", "CCM", "00", "BHZ"    # Cathedral Caves, MO (example)

# UTC start/end for a small demo (try 7–30 days for nicer stacks)
STARTTIME = "2020-01-01"
ENDTIME   = "2020-01-08"

# Band of interest (microseism example)
FREQMIN = 0.05   # Hz
FREQMAX = 0.20   # Hz

# Per-window settings
WINDOW_LENGTH = 24 * 3600  # seconds (1 day windows)
MAX_LAG = 600              # seconds of correlation lag to keep on each side (±MAX_LAG)
TIME_NORMALIZATION = "onebit"   # "onebit", "ram", or None
RAM_WINDOW_SEC = 10             # for "ram" (running absolute-mean) in seconds

# Remove instrument to velocity; pre-filter for deconvolution
PRE_FILT = (0.005, 0.006, 45.0, 50.0)  # (f1, f2, f3, f4) in Hz
WATER_LEVEL = 60                       # response water level
```

## Imports

```
import numpy as np
from obspy import UTCDateTime, Stream
from obspy.clients.fdsn import Client
from obspy.signal.util import next_pow_2
from scipy.signal import fftconvolve
import matplotlib.pyplot as plt
from tqdm import tqdm
```

## Downloading and PreProcessing

```
def fetch_stream(net, sta, loc, cha, t0, t1, client_name="IRIS", attach_resp=True):
    """
    Download waveforms and (optionally) attach response for later removal.
    """
    client = Client(client_name)
    st = client.get_waveforms(net, sta, loc, cha, t0, t1, attach_response=attach_resp)
    st.merge(fill_value="interpolate")
    return st

def basic_preprocess(st, freqmin, freqmax):
    """
    Detrend, demean, taper, band-pass. Operates in-place on a copy.
    """
    st = st.copy()
    st.detrend("linear")
    st.detrend("demean")
    st.taper(max(10.0 / st[0].stats.sampling_rate, 0.05))
    st.filter("bandpass", freqmin=freqmin, freqmax=freqmax, corners=4, zerophase=True)
    return st

def remove_instrument(st, pre_filt, water_level, output="VEL"):
    """
    Remove instrument response to velocity (default). Requires attached response.
    """
    st = st.copy()
    for tr in st:
        tr.remove_response(pre_filt=pre_filt, output=output, water_level=water_level)
    return st
```

## Amplitude stabilization (time-domain) and spectral whitening

```
def time_normalize(data, fs, mode="onebit", ram_win_sec=10):
    """
    Returns a copy of data with time-domain normalization applied.
    - onebit: sign(data)
    - ram: running absolute-mean normalization
    """
    x = data.astype(np.float64).copy()
    if mode is None:
        return x
    if mode.lower() == "onebit":
        return np.sign(x)
    if mode.lower() == "ram":
        n = max(1, int(ram_win_sec * fs))
        # Running absolute mean via convolution
        kern = np.ones(n, dtype=np.float64) / n
        denom = np.convolve(np.abs(x), kern, mode="same")
        denom[denom == 0] = 1.0
        return x / denom
    raise ValueError("Unknown time normalization mode")

def spectral_whiten(x, fs, fmin, fmax, eps=1e-10, smooth_bins=5):
    """
    Simple spectral whitening:
    - flatten amplitude within [fmin, fmax]
    - zero outside band (with tiny tapers to reduce ringing).
    """
    n = len(x)
    nfft = next_pow_2(2 * n)
    X = np.fft.rfft(x, n=nfft)
    freqs = np.fft.rfftfreq(nfft, d=1.0/fs)

    amp = np.abs(X)
    # Simple moving-average smoothing to avoid spikes
    if smooth_bins > 1:
        w = np.ones(smooth_bins) / smooth_bins
        amp = np.convolve(amp, w, mode="same")

    band = (freqs >= fmin) & (freqs <= fmax)
    Xw = np.zeros_like(X)
    Xw[band] = X[band] / (amp[band] + eps)

    # Inverse FFT and truncate back to original length
    xw = np.fft.irfft(Xw, n=nfft)[:n]
    return xw
```

## Cross-correlation (FFT-based, full lags → trimmed to ±MAX_LAG)

```
def xcorr_full(x, y):
    """
    Full linear cross-correlation using FFT.
    Returns ccf (length 2N-1) with lags from -(N-1) .. +(N-1):
       r_xy[k] = sum_n x[n] * y[n + k]
    """
    # Using convolution with reversed y is equivalent to correlation.
    return fftconvolve(x, y[::-1], mode="full")

def trim_to_maxlag(ccf, fs, max_lag):
    """
    Center the full CCF at zero lag and return only the ±max_lag window.
    """
    n = len(ccf)
    lags = np.arange(- (n - 1)//2, (n - 1)//2 + 1)  # placeholder; we’ll compute properly below

    # For full-mode length L=2N-1: zero lag index = N-1
    zero_idx = (n - 1) // 2
    half = int(max_lag * fs)
    left = max(0, zero_idx - half)
    right = min(n, zero_idx + half + 1)
    ccf_trim = ccf[left:right]
    # Build lag-time vector in seconds
    lags_sec = np.arange(left - zero_idx, right - zero_idx) / fs
    return ccf_trim, lags_sec
```

## The per-day pipeline (station pair → one window → one CCF)

```
def process_one_window(t0, window_length, pair, freqmin, freqmax, pre_filt, water_level,
                       time_norm_mode, ram_win_sec, max_lag):
    """
    Fetch, remove instrument, preprocess, normalize/whiten, correlate.
    pair: ((net, sta, loc, cha), (net, sta, loc, cha))
    Returns (ccf_trim, lags_sec) or (None, None) if fetch fails.
    """
    t1 = t0 + window_length
    (n1, s1, l1, c1), (n2, s2, l2, c2) = pair
    try:
        st1 = fetch_stream(n1, s1, l1, c1, t0, t1, attach_resp=True)
        st2 = fetch_stream(n2, s2, l2, c2, t0, t1, attach_resp=True)
    except Exception as e:
        print(f"[WARN] Fetch failed {t0}—{t1}: {e}")
        return None, None

    # Ensure single, matching channels
    if len(st1) == 0 or len(st2) == 0:
        return None, None
    tr1, tr2 = st1[0], st2[0]
    fs1, fs2 = tr1.stats.sampling_rate, tr2.stats.sampling_rate
    if not np.isclose(fs1, fs2):
        # resample to min(fs1, fs2)
        fs = min(fs1, fs2)
        st1.resample(fs)
        st2.resample(fs)
        tr1, tr2 = st1[0], st2[0]
    fs = tr1.stats.sampling_rate

    # Remove instrument, band-pass, etc.
    st1 = remove_instrument(st1, PRE_FILT, WATER_LEVEL, output="VEL")
    st2 = remove_instrument(st2, PRE_FILT, WATER_LEVEL, output="VEL")
    st1 = basic_preprocess(st1, freqmin, freqmax)
    st2 = basic_preprocess(st2, freqmin, freqmax)

    x = tr1 = st1[0].data.astype(np.float64)
    y = tr2 = st2[0].data.astype(np.float64)

    # Time-domain normalization (optional)
    if time_norm_mode:
        x = time_normalize(x, fs, mode=time_norm_mode, ram_win_sec=ram_win_sec)
        y = time_normalize(y, fs, mode=time_norm_mode, ram_win_sec=ram_win_sec)

    # Spectral whitening (optional but common in ambient noise)
    x = spectral_whiten(x, fs, freqmin, freqmax)
    y = spectral_whiten(y, fs, freqmin, freqmax)

    # Correlate (full) and trim to ±MAX_LAG
    ccf_full = xcorr_full(x, y)
    ccf_trim, lags_sec = trim_to_maxlag(ccf_full, fs, max_lag)

    # Normalize CCF by max abs to stabilize stacking (optional)
    m = np.max(np.abs(ccf_trim)) or 1.0
    ccf_trim = ccf_trim / m
    return ccf_trim, lags_sec
```

## Multi-day loop + stacking

```
def daterange(starttime, endtime, step_seconds):
    t0 = UTCDateTime(starttime)
    t1 = UTCDateTime(endtime)
    while t0 + step_seconds <= t1:
        yield t0
        t0 += step_seconds

def stack_ccfs(ccf_list):
    arr = np.vstack(ccf_list)  # shape: (nwin, nlag)
    return np.nanmean(arr, axis=0)

def run_pair_pipeline():
    pair = ((NETWORK1, STATION1, LOCATION1, CHANNEL1),
            (NETWORK2, STATION2, LOCATION2, CHANNEL2))
    ccfs = []
    lags = None

    for t0 in tqdm(list(daterange(STARTTIME, ENDTIME, WINDOW_LENGTH))):
        ccf, lags_sec = process_one_window(
            t0, WINDOW_LENGTH, pair, FREQMIN, FREQMAX, PRE_FILT, WATER_LEVEL,
            TIME_NORMALIZATION, RAM_WINDOW_SEC, MAX_LAG
        )
        if ccf is not None:
            ccfs.append(ccf)
            if lags is None:
                lags = lags_sec

    if not ccfs:
        raise RuntimeError("No windows succeeded. Try different dates/stations/channels.")

    ccf_stack = stack_ccfs(ccfs)
    return lags, ccf_stack, ccfs
```

## Plotting & quick interpretation

```
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
    """
    Quick-look: imshow each per-window CCF as a row (not distance-sorted).
    """
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
```

## Run

```
if __name__ == "__main__":
    lags, ccf_stack, ccfs = run_pair_pipeline()
    plot_ccf(lags, ccf_stack, title=f"Stacked CCF: {NETWORK1}.{STATION1} ↔ {NETWORK2}.{STATION2}")
    plot_record_section(lags, ccfs, title="Per-window CCFs (each day)")
```

## Flip order to check causal vs. acausal

```
def run_flipped():
    global NETWORK1, STATION1, LOCATION1, CHANNEL1
    global NETWORK2, STATION2, LOCATION2, CHANNEL2
    # Swap station definitions
    (NETWORK1, STATION1, LOCATION1, CHANNEL1), (NETWORK2, STATION2, LOCATION2, CHANNEL2) = \
        (NETWORK2, STATION2, LOCATION2, CHANNEL2), (NETWORK1, STATION1, LOCATION1, CHANNEL1)
    lags, ccf_stack, _ = run_pair_pipeline()
    plot_ccf(lags, ccf_stack, title="Stacked CCF (flipped order)")
```

## Estimate a simple Rayleigh “apparent velocity” from the first peak (very rough)

```
def first_peak_time(ccf, lags, search_window=(5, 200)):
    """Find the max on the causal side within a search window (s)."""
    i0 = np.searchsorted(lags, search_window[0])
    i1 = np.searchsorted(lags, search_window[1])
    idx = i0 + np.argmax(ccf[i0:i1])
    return lags[idx], ccf[idx]

# Example usage (you’ll need an inter-station distance in km)
# dist_km =  ...   # e.g., from station metadata (obspy.geodetics.locations2degrees + kilometers2degrees)
# t_peak, amp = first_peak_time(ccf_stack, lags, search_window=(5, 300))
# c_app = dist_km * 1000.0 / t_peak   # m/s
# print(f"First-peak apparent velocity ~ {c_app:.1f} m/s")
```

## TA demo: end-to-end Python

```
# =========================
# 0) Transportable Array (TA) parameters
# =========================
from obspy import UTCDateTime

# Choose a place/time when TA operated nearby.
# Example A (Lower 48, 2011 near the Rockies):
CENTER_LAT, CENTER_LON = 40.0, -105.3   # Front Range (CO/WY)
STARTTIME = "2011-06-01"
ENDTIME   = "2011-06-08"

# Example B (Alaska TA, 2018 around Anchorage) — uncomment to try:
# CENTER_LAT, CENTER_LON = 61.2, -149.9
# STARTTIME = "2018-07-01"
# ENDTIME   = "2018-07-08"

# Search radius (degrees) for TA stations around the center:
MAXRADIUS_DEG = 6.0

# Prefer inter-station distances in this range (km) for microseism band:
MIN_DIST_KM, MAX_DIST_KM = 100.0, 400.0

# Band of interest (microseisms):
FREQMIN = 0.05   # Hz  (20 s)
FREQMAX = 0.20   # Hz  (5 s)

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
# 2) TA station discovery
# =========================
def find_ta_pair(center_lat, center_lon, start, end,
                 maxradius_deg=6.0, min_km=100.0, max_km=400.0,
                 preferred_channels=("BHZ", "HHZ")):
    """
    Query EarthScope FDSN for TA stations active in [start, end] within the search radius.
    Return two tuples: (net, sta, loc, cha, lat, lon) for a pair whose distance falls in [min_km, max_km].
    Preference is given to vertical BH/HH channels with location codes available.
    """
    client = Client("IRIS")  # EarthScope Data Services FDSN endpoint
    t0, t1 = UTCDateTime(start), UTCDateTime(end)

    # Get stations with vertical broadband channels; include channel level for loc/code
    inv = client.get_stations(network="TA",
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
        raise RuntimeError("Not enough TA stations found—try moving the center, expanding MAXRADIUS_DEG, or changing dates.")

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
    st = client.get_waveforms(net, sta, loc, cha, t0, t1, attach_response=attach_resp)
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
        st1 = fetch_stream(n1, s1, l1, c1, t0, t1, attach_response=True)
        st2 = fetch_stream(n2, s2, l2, c2, t0, t1, attach_response=True)
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
# 9) Run the TA pipeline
# =========================
def run_ta_demo():
    # Pick a TA station pair automatically
    A, B, dist_km = find_ta_pair(CENTER_LAT, CENTER_LON, STARTTIME, ENDTIME,
                                 maxradius_deg=MAXRADIUS_DEG,
                                 min_km=MIN_DIST_KM, max_km=MAX_DIST_KM)
    print(f"[TA] Using stations: {A[0]}.{A[1]}.{A[2]}.{A[3]} ↔ {B[0]}.{B[1]}.{B[2]}.{B[3]}  (~{dist_km:.1f} km apart)")

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
    plot_ccf(lags, ccf_stack, title=f"TA stacked CCF: {A[1]} ↔ {B[1]}  ({dist_km:.0f} km)")
    plot_record_section(lags, ccfs, title="TA per-window CCFs (each day)")
    return {"pair_A": A, "pair_B": B, "distance_km": dist_km, "lags": lags, "ccf_stack": ccf_stack, "ccfs": ccfs}

if __name__ == "__main__":
    _ = run_ta_demo()

```

## What STA/LTA does

It’s a simple onset detector: compare the recent signal energy (STA) to the background energy (LTA).

- During quiet time, both are similar → ratio ≈ 1.
- When a phase (P, S, LP, blast, etc.) arrives, short-term energy jumps while long-term background hasn’t caught up → ratio ≫ 1 → trigger “ON”.
- After the transient passes, the ratio falls → trigger “OFF”.

### Minimal workflow

1. Preprocess: remove mean/trend, band-pass to target band, optional instrument removal if needed.
2. Characteristic function: use x^2 or envelope.
3. Compute STA, LTA, ratio.
4. Hysteresis thresholds + coincidence across stations (optional).
5. Refine picks (e.g., maximum slope, polarity, travel-time consistency).

### Common pitfalls (and fixes)

- Wrong band → swamped by microseisms or cultural noise. Fix: pick a band where your phase is distinctive.
- LTA contamination by big events → elevated baseline. Fix: longer LTA, freeze LTA briefly after trigger, or robust baselines (median/MAD).
- Chatter around threshold. Fix: use Ton > Toff, apply a short refractory time after ON.
- Gain/clock issues across stations → spurious coincidence. Fix: check metadata, use per-station thresholds.

### Obspy one-liners

```
from obspy.signal.trigger import classic_sta_lta, recursive_sta_lta, trigger_onset

# Assume x_f is bandpassed, float64
sta_s, lta_s = 0.8, 12.0
nsta, nlta = int(sta_s*fs), int(lta_s*fs)

# Classic boxcar
cft = classic_sta_lta(x_f, nsta, nlta)

# Or recursive (Allen)
cft_rec = recursive_sta_lta(x_f, nsta, nlta)

# Get ON/OFF using hysteresis
on_off = trigger_onset(cft, 4.0, 2.0)   # Ton=4.0, Toff=2.0
# on_off rows are sample indices [on, off)
```


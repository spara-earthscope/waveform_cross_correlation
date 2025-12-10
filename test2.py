import numpy as np
import matplotlib.pyplot as plt
from obspy import read, UTCDateTime
from obspy.signal.cross_correlation import correlate, xcorr_max
from obspy.clients.fdsn import Client

# --- Step 1: Download or load data ---
# Replace with your own data or use IRIS client
client = Client("IRIS")
t = UTCDateTime("2020-01-01T00:00:00")
st1 = client.get_waveforms("IU", "ANMO", "00", "BHZ", t, t + 3600)  # 1 hour of data
st2 = client.get_waveforms("IU", "CCM", "00", "BHZ", t, t + 3600)

# --- Step 2: Preprocess data ---
def preprocess(st):
    st.detrend("linear")
    st.taper(max_percentage=0.05)
    st.filter("bandpass", freqmin=0.5, freqmax=2.0)  # 10-20 s period
    st.resample(1.0)  # Resample to 1 Hz
    return st

st1 = preprocess(st1)
st2 = preprocess(st2)

# --- Step 3: Cross-correlate ---
tr1 = st1[0]
tr2 = st2[0]
corr = correlate(tr1, tr2, shift=100, normalize="naive")

# --- Step 4: Plot the cross-correlation ---
plt.figure(figsize=(10, 4))
plt.plot(corr, 'k-', label="Cross-correlation")
plt.xlabel("Lag (s)")
plt.ylabel("Amplitude")
plt.title("Cross-correlation: ANMO-CCM")
plt.grid()
plt.legend()
plt.show()

# --- Step 5: Find the maximum correlation ---
shift, value = xcorr_max(corr, abs_max=True)
print(f"Max correlation at lag: {shift} samples, value: {value:.2f}")

# --- Step 6: Stack multiple windows (example with 3 windows) ---
# In practice, use many more windows (e.g., 1 day = 24 windows)
corr_stack = corr.copy()
for i in range(1, 3):
    st1_window = client.get_waveforms("IU", "ANMO", "00", "BHZ", t + i*3600, t + (i+1)*3600)
    st2_window = client.get_waveforms("IU", "CCM", "00", "BHZ", t + i*3600, t + (i+1)*3600)
    st1_window = preprocess(st1_window)
    st2_window = preprocess(st2_window)
    corr_window = correlate(st1_window[0], st2_window[0], shift=100, normalize="naive")
    corr_stack += corr_window

corr_stack /= 3  # Average

# --- Step 7: Plot the stacked cross-correlation ---
plt.figure(figsize=(10, 4))
plt.plot(corr_stack, 'r-', label="Stacked cross-correlation")
plt.xlabel("Lag (s)")
plt.ylabel("Amplitude")
plt.title("Stacked Cross-correlation: ANMO-CCM (3 windows)")
plt.grid()
plt.legend()
plt.show()

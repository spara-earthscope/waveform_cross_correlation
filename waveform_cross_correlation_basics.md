# Seismic Matched Filtering with Waveform Cross Correlation

Goal: Find volcanic tremors by cross-correlating it to a known long period event.

## Basics

**Theory**: Matched filtering is a detection technique where you slide a known waveform (a “template”) across continuous data and compute a similarity score at each time shift. If the data contains a waveform shaped like your template (even if it’s much weaker than the noise), the score spikes.

Mathematically, the matched filter is just the cross-correlation of the data with the time-reversed template. In white Gaussian noise, it’s provably optimal for maximizing signal-to-noise ratio (SNR).

- Noise is (mostly) random; it doesn’t reinforce at a specific lag.
- A repeated signal has the same shape every time. When aligned with its template, the inner product (correlation) adds up coherently, producing a large peak.
- Normalizing the correlation makes the detector insensitive to amplitude differences (helpful when events have similar shapes but different sizes).

**Method**: Correlate a known template with noisy data. A spike indicates a detection.

**Process**:

1. 

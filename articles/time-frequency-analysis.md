# Time-Frequency Analysis

## Introduction

Time-frequency analysis is essential for understanding how oscillatory
brain activity changes over time in response to stimuli or during
cognitive tasks. Unlike classical ERP analysis that focuses on
time-domain averages, time-frequency methods reveal dynamic spectral
changes in specific frequency bands (delta, theta, alpha, beta, gamma)
that are often masked in traditional averaging.

Common neuroscience applications include:

- **Event-Related Spectral Perturbation (ERSP)**: Power changes relative
  to baseline in specific frequency bands, revealing induced oscillatory
  activity that may not be time-locked
- **Intertrial Coherence (ITC)**: Phase consistency across trials,
  capturing evoked oscillatory responses
- **Theta-gamma coupling**: Investigating cross-frequency interactions
  during memory encoding
- **Mu/beta suppression**: Motor imagery and action observation studies
- **Alpha desynchronization**: Attentional modulation and visual
  processing

PhysioEEG provides three complementary time-frequency decomposition
methods, each with distinct advantages depending on your research
question, signal characteristics, and required time-frequency
resolution.

## Setup

``` r

library(PhysioEEG)
library(SummarizedExperiment)

# The single-signal transforms (Morlet, STFT, multitaper) operate on 2D
# continuous data (time x channels); ERSP and ITC operate on 3D epoched data
# (time x channels x trials) and compute the wavelet internally.

# 2D continuous data (19-channel 10-20 layout, 4 s at 500 Hz)
eeg <- make_eeg(n_channels = 19, n_time = 2000, sr = 500)

# 3D epoched ERP data (40 trials of 1 s) for ERSP and ITC
eeg_epochs <- make_eeg_erp(
  n_channels = 19,
  n_epochs = 40,
  sr = 500,
  epoch_sec = 1.0
)

# Inspect dimensions
dim(assay(eeg, "raw"))         # timepoints x channels (continuous)
dim(assay(eeg_epochs, "raw"))  # timepoints x channels x trials (epoched)
samplingRate(eeg)              # 500 Hz
```

## Morlet Wavelet Transform

The Morlet wavelet transform provides excellent time-frequency
resolution trade-off control via the number of cycles parameter. It
convolves the signal with complex Morlet wavelets at specified
frequencies, yielding both power and phase information.

**Parameter selection:**

- **Frequencies**: Define your frequency range based on your research
  question. For general cognitive studies, 4-40 Hz captures theta
  through low gamma. For detailed alpha analysis, use finer spacing
  (e.g., 8-13 Hz in 0.5 Hz steps).
- **n_cycles**: Controls time vs frequency resolution trade-off. Lower
  values (3-5 cycles) provide better temporal precision but broader
  frequency resolution; higher values (7-12 cycles) provide sharper
  frequency resolution but smear events in time. A common compromise is
  to use frequency-dependent cycles: `n_cycles = f / 3` (e.g., 4 Hz →
  1.3 cycles, 40 Hz → 13 cycles).

``` r

# Detailed alpha-band analysis with higher frequency resolution
eeg <- eegMorletWavelet(
  eeg,
  frequencies = seq(8, 13, by = 0.5),  # Fine-grained alpha range
  n_cycles = 7,                         # More cycles for sharper frequency precision
  assay_name = "raw",
  output_assay = "alpha_power"
)

# Broad-band Morlet transform (theta-gamma). eegMorletWavelet() uses a fixed
# (scalar) number of cycles across frequencies. Stored under the default name
# "wavelet_power", this is the product visualized in the spectrogram section.
eeg <- eegMorletWavelet(
  eeg,
  frequencies = seq(4, 40, by = 2),    # 4-40 Hz in 2 Hz steps
  n_cycles = 5,                         # good time-frequency compromise
  assay_name = "raw",
  output_assay = "wavelet_power"
)

# Output (stored in metadata): timepoints x frequencies x channels
dim(metadata(eeg)$wavelet_power)
```

**When to use Morlet wavelets:**

- Default choice for most EEG time-frequency analysis
- Excellent for event-related designs where temporal precision matters
- Ideal when you need consistent time-frequency resolution across the
  spectrum
- Preferred for cross-frequency coupling analysis

## Short-Time Fourier Transform

STFT divides the signal into short overlapping windows and computes the
FFT for each segment. It provides uniform time-frequency resolution
across all frequencies, determined by the window size.

**Window size trade-offs:**

- **Window size**: Larger windows (e.g., 1000 ms) provide better
  frequency resolution but poor temporal precision. Smaller windows
  (e.g., 100 ms) capture rapid changes but have coarse frequency
  resolution. A 200-500 ms window is typical for cognitive EEG.
- **Overlap**: Higher overlap (75-90%) provides smoother time-frequency
  representations but increases computation. 50% overlap is a reasonable
  default.
- **Window type**: Hanning window reduces spectral leakage and is the
  standard choice. Hamming window has slightly better frequency
  resolution but more leakage.

``` r

# Standard STFT with 250 ms windows (good balance for cognitive EEG)
eeg <- eegSTFT(
  eeg,
  window_sec = 0.25,      # 250 ms window
  overlap = 0.75,         # 75% overlap for smooth representation
  window_type = "hanning",
  assay_name = "raw",
  output_assay = "stft_power"
)

# Output (stored in metadata): timepoints x frequencies x channels
dim(metadata(eeg)$stft_power)

# For high temporal precision (e.g., rapid stimulus changes)
eeg <- eegSTFT(
  eeg,
  window_sec = 0.1,       # 100 ms windows
  overlap = 0.50,         # 50% overlap
  window_type = "hanning",
  assay_name = "raw",
  output_assay = "stft_fast"
)

# For high frequency resolution (e.g., resting-state narrow-band analysis)
eeg <- eegSTFT(
  eeg,
  window_sec = 1.0,       # 1000 ms windows
  overlap = 0.90,         # High overlap to maintain some temporal detail
  window_type = "hamming",
  assay_name = "raw",
  output_assay = "stft_fine"
)
```

**When to use STFT:**

- When uniform time-frequency resolution is desired
- Continuous EEG or resting-state data
- When computational speed is a priority (faster than wavelets for many
  frequencies)
- When you need direct control over time-frequency resolution in
  physical units

## Multitaper Spectral Analysis

The Thomson multitaper method uses multiple orthogonal tapers (windowing
functions) to reduce variance in spectral estimates while maintaining
good frequency resolution. It provides robust spectral estimates,
especially for stationary or quasi-stationary signals.

**Parameter selection:**

- **Bandwidth**: Determines frequency smoothing and number of available
  tapers. Bandwidth of 2-4 Hz is typical for EEG. Higher bandwidth
  provides more tapers (lower variance) but coarser frequency
  resolution.
- **n_tapers**: Number of tapers to use. More tapers reduce variance but
  increase smoothing. Typically set to `floor(2 * bandwidth * T - 1)`
  where T is the time window duration in seconds.

``` r

# Multitaper analysis with 3 Hz bandwidth
# Good for reducing noise in oscillatory power estimates
eeg <- eegMultitaper(
  eeg,
  bandwidth = 3,          # 3 Hz smoothing bandwidth
  n_tapers = 5,           # Use 5 Slepian tapers
  assay_name = "raw",
  output_assay = "multitaper_power"
)

# Output (stored in metadata): timepoints x frequencies x channels
dim(metadata(eeg)$multitaper_power)

# For high-frequency resolution (narrow-band analysis)
eeg <- eegMultitaper(
  eeg,
  bandwidth = 1.5,        # Narrow bandwidth for fine frequency resolution
  n_tapers = 3,
  assay_name = "raw",
  output_assay = "multitaper_narrow"
)
```

**When to use multitaper:**

- Resting-state or continuous EEG with quasi-stationary spectral content
- When you need robust spectral estimates with minimal variance
- Narrow-band oscillation detection (e.g., individual alpha frequency)
- When signal-to-noise ratio is low
- Spectral analysis of relatively long, stable epochs

## Event-Related Spectral Perturbation (ERSP)

ERSP quantifies time-frequency power changes relative to a baseline
period, typically pre-stimulus. It reveals both event-related
synchronization (ERS, power increase) and desynchronization (ERD, power
decrease), capturing induced oscillatory activity that may not be
strictly phase-locked to the stimulus.

**Baseline normalization methods:**

- **Decibel (dB)**: `10 * log10(power / baseline)` - symmetric,
  interpretable as percent change on log scale
- **Percent change**: `100 * (power - baseline) / baseline` - intuitive
  but asymmetric
- **Z-score**: `(power - mean_baseline) / sd_baseline` - accounts for
  variance, good for statistics

``` r

# eegERSP() takes the 3D epoched data and computes the Morlet wavelet and
# baseline normalization internally. The baseline window is given in samples.
freqs <- seq(4, 40, by = 2)
eeg_epochs <- eegERSP(
  eeg_epochs,
  baseline = c(1, 50),   # baseline = first 50 samples (100 ms) of each epoch
  frequencies = freqs,
  n_cycles = 5,
  assay_name = "raw"
)

# ERSP stored in metadata (already averaged across trials)
ersp_data <- metadata(eeg_epochs)$ersp_data
dim(ersp_data)  # timepoints x frequencies x channels

# Extract single-channel ERSP for visualization
pz_idx <- match("Pz", colData(eeg_epochs)$label)
pz_ersp <- ersp_data[, , pz_idx]
dim(pz_ersp)  # timepoints x frequencies
```

**Interpreting ERSP:**

- **Positive values (ERS)**: Power increase relative to baseline.
  Example: frontal theta increase during memory encoding, gamma increase
  during sensory processing.
- **Negative values (ERD)**: Power decrease relative to baseline.
  Example: central mu/beta suppression during motor preparation,
  occipital alpha decrease during visual attention.
- **Time course**: ERSP reveals when spectral changes occur,
  distinguishing early sensory from late cognitive processes.
- **Spatial distribution**: Topographic patterns indicate functional
  networks (e.g., frontal theta vs posterior alpha).

## Intertrial Coherence (ITC)

ITC measures phase consistency across trials at each time-frequency
point, ranging from 0 (random phase) to 1 (perfect phase alignment).
Unlike ERSP which captures power changes, ITC specifically identifies
phase-locked (evoked) oscillatory activity.

**Key concepts:**

- **ITC = 1**: All trials have identical phase at that time-frequency
  point (perfect evoked response)
- **ITC = 0**: Random phase across trials (purely induced activity)
- **ITC vs ERP**: High ITC at low frequencies corresponds to traditional
  ERP components. However, ITC can also detect phase-locked
  higher-frequency oscillations invisible in the averaged ERP.

``` r

# eegITC() computes the complex Morlet transform per trial internally and
# measures phase consistency across trials from the 3D epoched data.
eeg_epochs <- eegITC(
  eeg_epochs,
  frequencies = freqs,
  n_cycles = 5,
  assay_name = "raw"
)

# ITC stored in metadata
itc_data <- metadata(eeg_epochs)$itc_data
dim(itc_data)  # timepoints x frequencies x channels

# Extract channel-specific ITC (Cz is a central fronto-central lead)
cz_idx <- match("Cz", colData(eeg_epochs)$label)
cz_itc <- itc_data[, , cz_idx]
dim(cz_itc)  # timepoints x frequencies

# High ITC at low frequencies captures ERP phase-locking.
# Theta-band (4-8 Hz) ITC time course:
theta_cols <- which(freqs >= 4 & freqs <= 8)
theta_itc <- rowMeans(cz_itc[, theta_cols, drop = FALSE])
plot(theta_itc, type = "l", xlab = "Time (samples)", ylab = "ITC")
```

**Interpreting ITC:**

- **High ITC + high ERSP**: Evoked power increase (both phase-locked and
  power increase)
- **High ITC + low ERSP**: Phase-locked activity without net power
  change (phase reset)
- **Low ITC + high ERSP**: Induced power increase without phase-locking
  (e.g., stimulus-triggered but temporally jittered oscillations)
- **Low-frequency ITC peaks**: Correspond to traditional ERP components
  (N1, P3, etc.)
- **High-frequency ITC**: May indicate early sensory-evoked gamma or
  beta responses

## Visualization with eegPlotSpectrogram

Time-frequency heatmaps (spectrograms) are the standard visualization
for ERSP, ITC, and raw power. PhysioEEG provides customizable
spectrogram plotting with sensible defaults.

``` r

# Plot the Morlet power spectrogram for channel Pz
eegPlotSpectrogram(
  eeg,
  channel = "Pz",
  freq_range = c(4, 40),      # Focus on theta-gamma range
  time_range = NULL,          # full time range (time axis is in seconds)
  log_power = TRUE,           # convert power to dB scale
  palette = "viridis",
  assay_name = "wavelet_power"
)

# Plot the spectrogram for a central lead (Cz)
eegPlotSpectrogram(
  eeg,
  channel = "Cz",
  freq_range = c(4, 30),
  time_range = NULL,
  log_power = TRUE,
  palette = "viridis",
  assay_name = "wavelet_power"
)

# Zoom into the alpha band for an occipital channel
eegPlotSpectrogram(
  eeg,
  channel = "O1",
  freq_range = c(8, 13),      # Alpha band detail
  time_range = NULL,
  log_power = TRUE,           # Convert to dB scale
  palette = "magma",
  assay_name = "wavelet_power"
)
```

**Visualization tips:**

- **ERSP**: Use diverging colormaps (RdBu, RdYlBu) with center at zero
  to distinguish ERS from ERD
- **ITC**: Use sequential colormaps (viridis, magma, plasma) since
  values range from 0 to 1
- **Raw power**: Apply log scaling to compress dynamic range and make
  oscillations visible
- **Frequency range**: Narrow ranges (e.g., 8-13 Hz for alpha) reveal
  fine structure; broad ranges (4-80 Hz) show multi-band dynamics
- **Time range**: Include pre-stimulus baseline to contextualize
  post-stimulus changes

## Method Comparison Guide

Choosing the right time-frequency method depends on your data
characteristics, research question, and computational resources.

| **Method** | **Time Resolution** | **Frequency Resolution** | **Best For** | **Computational Cost** |
|----|----|----|----|----|
| **Morlet Wavelet** | Variable (tunable via n_cycles) | Variable (tunable via n_cycles) | Event-related designs, flexible resolution trade-off, cross-frequency coupling | Moderate-High (depends on \# frequencies) |
| **STFT** | Uniform (set by window) | Uniform (set by window) | Continuous EEG, resting-state, fast computation, uniform resolution | Low-Moderate |
| **Multitaper** | Moderate | High (good for narrow bands) | Resting-state, low SNR, narrow-band oscillations, robust spectral estimates | Moderate |

**Decision tree:**

1.  **Epoched event-related data with precise timing requirements**: Use
    Morlet wavelets with frequency-dependent n_cycles for optimal
    time-frequency resolution matching signal characteristics.

2.  **Continuous or resting-state EEG**: Use STFT for speed or
    multitaper for robust spectral estimates with low variance.

3.  **Narrow-band oscillation detection (e.g., individual alpha
    frequency)**: Use multitaper with low bandwidth for superior
    frequency resolution.

4.  **Broad frequency range with varying dynamics**: Use Morlet wavelets
    with frequency-dependent parameters (low n_cycles for low
    frequencies, high n_cycles for high frequencies).

5.  **Low signal-to-noise ratio**: Use multitaper to reduce spectral
    variance, or average multiple Morlet wavelet trials.

6.  **Computational constraints**: Use STFT, which is faster than
    wavelets when analyzing many frequencies.

**Frequency resolution comparison:**

- **Wavelet**: Resolution = f / n_cycles. Example: 10 Hz with 5 cycles →
  2 Hz resolution.
- **STFT**: Resolution = sampling_rate / window_size. Example: 500 Hz
  with 250-sample window → 2 Hz resolution.
- **Multitaper**: Resolution ≈ bandwidth. Example: 3 Hz bandwidth → ~3
  Hz resolution.

**Time resolution comparison:**

- **Wavelet**: Resolution = n_cycles / f. Example: 10 Hz with 5 cycles →
  0.5 s resolution.
- **STFT**: Resolution = window_size / sampling_rate. Example: 250
  samples at 500 Hz → 0.5 s resolution.
- **Multitaper**: Resolution depends on window size and taper
  properties, generally moderate.

## Session Info

``` r

sessionInfo()
```

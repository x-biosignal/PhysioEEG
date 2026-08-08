# Compute cross-spectral density between two signals using Welch's method

Segments the signals, applies a Hanning window, computes FFTs, and
averages the cross-spectrum and auto-spectra across windows.

## Usage

``` r
.welch_cross_spectrum(x, y, sr, window_sec, overlap)
```

## Arguments

- x:

  Numeric vector, first signal.

- y:

  Numeric vector, second signal.

- sr:

  Sampling rate in Hz.

- window_sec:

  Window length in seconds.

- overlap:

  Overlap fraction (0 to 1).

## Value

A list with components:

- Sxy:

  Complex vector of cross-spectral density (positive freqs).

- Sxx:

  Numeric vector of auto-spectral density of x.

- Syy:

  Numeric vector of auto-spectral density of y.

- freqs:

  Numeric vector of frequencies in Hz.

- n_windows:

  Integer number of windows averaged.

- Sxy_segments:

  List of per-window cross-spectra (complex vectors).

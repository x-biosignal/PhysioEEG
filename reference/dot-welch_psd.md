# Welch power spectral density (internal)

Averaged periodogram (Hann window, 50 percent overlap) giving a stable
one-sided PSD estimate, used by the clinical band-power routines.

## Usage

``` r
.welch_psd(signal, sr, nperseg = NULL, noverlap = NULL)
```

## Arguments

- signal:

  Numeric signal.

- sr:

  Sampling rate in Hz.

- nperseg:

  Segment length in samples (default: about a one-second power-of-two
  segment, capped at the signal length).

- noverlap:

  Overlap in samples (default: half of `nperseg`).

## Value

A list with `freq` (Hz) and `psd` (one-sided).

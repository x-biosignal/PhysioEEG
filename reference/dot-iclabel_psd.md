# Welch power spectral density (internal)

Welch power spectral density (internal)

## Usage

``` r
.iclabel_psd(s, sr, nperseg = 512)
```

## Arguments

- s:

  Numeric signal.

- sr:

  Sampling rate in Hz.

- nperseg:

  Segment length (default: 512).

## Value

List with `f` (frequencies) and `p` (one-sided power).

# IIR Butterworth filter wrapper

Applies a Butterworth IIR filter using
[`signal::butter`](https://rdrr.io/pkg/signal/man/butter.html) and
[`signal::filtfilt`](https://rdrr.io/pkg/signal/man/filtfilt.html) for
zero-phase filtering. Falls back to FIR if the signal package is not
available.

## Usage

``` r
.iir_butterworth(data_matrix, sr, lowcut = NULL, highcut = NULL, order = 4)
```

## Arguments

- data_matrix:

  Numeric matrix (time x channels).

- sr:

  Sampling rate in Hz.

- lowcut:

  Low cutoff frequency in Hz (NULL for lowpass).

- highcut:

  High cutoff frequency in Hz (NULL for highpass).

- order:

  Filter order (default: 4).

## Value

Filtered numeric matrix (same dimensions).

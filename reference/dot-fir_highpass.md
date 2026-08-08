# Windowed-sinc FIR highpass filter

Applies a windowed-sinc FIR highpass filter using spectral inversion of
a lowpass filter. The lowpass kernel at the cutoff frequency is computed
and then inverted to create a highpass response.

## Usage

``` r
.fir_highpass(signal, sr, cutoff, order = NULL)
```

## Arguments

- signal:

  Numeric vector of the input signal.

- sr:

  Sampling rate in Hz.

- cutoff:

  Cutoff frequency in Hz.

- order:

  Filter order. If `NULL`, auto-selected.

## Value

Numeric vector of the filtered signal (same length as input).

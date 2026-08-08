# Compute instantaneous phase via Hilbert transform

Constructs the analytic signal by zeroing negative frequency components
of the FFT and doubling positive frequencies, then returns the
instantaneous phase (argument of the analytic signal).

## Usage

``` r
.hilbert_phase(signal)
```

## Arguments

- signal:

  Numeric vector of the input signal.

## Value

Numeric vector of instantaneous phase in radians (-pi, pi\].

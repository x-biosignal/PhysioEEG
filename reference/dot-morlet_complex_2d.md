# Compute Morlet wavelet complex coefficients for 2D data (internal)

Computes complex wavelet coefficients for a 2D matrix (time x channels)
at specified frequencies using Morlet wavelets with FFT convolution.
Used by `eegITC` to extract phase information.

## Usage

``` r
.morlet_complex_2d(data, sr, frequencies, n_cycles)
```

## Arguments

- data:

  Numeric matrix (time x channels).

- sr:

  Sampling rate in Hz.

- frequencies:

  Numeric vector of frequencies.

- n_cycles:

  Number of wavelet cycles.

## Value

3D complex array (time x frequencies x channels).

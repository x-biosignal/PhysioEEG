# Compute Morlet wavelet power for 2D data (internal)

Computes wavelet power (\|W\|^2) for a 2D matrix (time x channels) at
specified frequencies using Morlet wavelets with FFT convolution.

## Usage

``` r
.morlet_power_2d(data, sr, frequencies, n_cycles)
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

3D array (time x frequencies x channels) of power values.

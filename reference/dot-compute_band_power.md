# Compute band power for a single signal (Welch PSD)

Integrates a Welch (averaged-periodogram) power spectral density over
the band, replacing the earlier single-window FFT for a lower-variance
estimate.

## Usage

``` r
.compute_band_power(signal, sr, low, high, method = "welch", nperseg = NULL)
```

## Arguments

- signal:

  Numeric vector of the input signal.

- sr:

  Sampling rate in Hz.

- low:

  Low frequency bound in Hz.

- high:

  High frequency bound in Hz.

- method:

  Spectral estimator; currently `"welch"`.

- nperseg:

  Optional Welch segment length in samples.

## Value

Numeric scalar of the integrated spectral power in the band.

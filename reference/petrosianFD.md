# Petrosian Fractal Dimension

Computes the Petrosian fractal dimension of a one-dimensional signal
(Petrosian, 1995) – a fast waveform-complexity measure based on the
number of sign changes in the signal's derivative (i.e. the density of
local extrema). A smooth, oscillatory signal has few derivative sign
changes and a Petrosian FD near 1; a complex or noise-like signal has
many and a higher FD. \\P = \log\_{10} N / (\log\_{10} N + \log\_{10}(N
/ (N + 0.4\\N\_\delta)))\\, where \\N\\ is the signal length and
\\N\_\delta\\ the number of sign changes in the first difference.

## Usage

``` r
petrosianFD(x)
```

## Arguments

- x:

  A numeric vector (the time series).

## Value

A single numeric value, the Petrosian fractal dimension (typically in
\\\[1, \sim1.1\]\\ for physiological signals).

## References

Petrosian, A. (1995). "Kolmogorov complexity of finite sequences and
recognition of different preictal EEG patterns." In *Proceedings of the
Eighth IEEE Symposium on Computer-Based Medical Systems*, 212–217.
[doi:10.1109/CBMS.1995.465426](https://doi.org/10.1109/CBMS.1995.465426)

## See also

[`svdEntropy`](https://x-biosignal.github.io/PhysioEEG/reference/svdEntropy.md)
for the singular-value complexity view,
[`eegComplexity`](https://x-biosignal.github.io/PhysioEEG/reference/eegComplexity.md)
for the multi-measure wrapper.

## Examples

``` r
petrosianFD(sin(seq(0, 20 * pi, length.out = 500)))  # low (smooth oscillation)
#> [1] 1.002561
set.seed(1); petrosianFD(rnorm(500))                 # higher (noise-like)
#> [1] 1.039285
```

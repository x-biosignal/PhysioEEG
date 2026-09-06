# Slope Entropy

Computes the slope entropy (SlopEn) of a one-dimensional signal
(Cuesta-Frau, 2019) – a symbolic complexity measure defined on the
*slopes* between successive samples. Each first difference is turned
into the *angle* \\\theta = \arctan(\Delta x)\\ (in degrees) and mapped
to one of five symbols by two angle thresholds: steep-up (\\\>
\gamma\\), gentle-up (\\(\delta, \gamma\]\\), flat (\\\[-\delta,
\delta\]\\), gentle-down and steep-down; SlopEn is the Shannon entropy
of the resulting length-(m-1) slope patterns. By combining a coarse
(steep/gentle) and fine (flat) angular resolution it captures the shape
of local trends, a mechanism distinct from the amplitude, NCDF-pattern,
template-matching and increment-magnitude entropies.

## Usage

``` r
slopeEntropy(x, dimension = 3L, thresholds = c(0.1, 45), delay = 1L)
```

## Arguments

- x:

  A numeric vector (the time series).

- dimension:

  Embedding dimension \\m\\ (default 3); slope patterns have length
  `dimension - 1`.

- thresholds:

  Increasing positive angle thresholds (degrees) that split the slope
  into symbol classes (default `c(0.1, 45)`: a flat band within +/-0.1
  deg, a gentle band up to +/-45 deg, and a steep class beyond).

- delay:

  Time delay (lag) in samples for the difference (default 1).

## Value

A single numeric value, the slope entropy (bits).

## Details

This implementation reproduces `NeuroKit2`'s `entropy_slope`
bit-for-bit.

## References

Cuesta-Frau, D. (2019). "Slope Entropy: A New Time Series Complexity
Estimator Based on Both Symbolic Patterns and Amplitude Information."
*Entropy*, 21(12), 1167.
[doi:10.3390/e21121167](https://doi.org/10.3390/e21121167)

## See also

[`incrementEntropy`](https://x-biosignal.github.io/PhysioEEG/reference/incrementEntropy.md),
[`dispersionEntropy`](https://x-biosignal.github.io/PhysioEEG/reference/dispersionEntropy.md),
[`fuzzyEntropy`](https://x-biosignal.github.io/PhysioEEG/reference/fuzzyEntropy.md).

## Examples

``` r
set.seed(1)
slopeEntropy(rnorm(1000), dimension = 3)
#> [1] 3.865113
```

# Fuzzy Entropy

Computes the fuzzy entropy (FuzzyEn) of a one-dimensional signal (Chen
et al., 2007) – the *fuzzy* generalization of the sample entropy. Sample
entropy counts template matches with a crisp threshold (a Heaviside
step: two vectors either match within tolerance \\r\\ or do not); fuzzy
entropy replaces that hard step with a smooth exponential membership
\\\mu = \exp(-(d^{n})/r)\\, so near-matches contribute partially. This
makes FuzzyEn markedly more robust and continuous than sample entropy,
especially on short or noisy segments where the crisp count is unstable.
Each embedding vector has its own baseline (mean) removed before the
Chebyshev distance is taken, following Chen et al.

## Usage

``` r
fuzzyEntropy(x, dimension = 2L, r = 0.2, delay = 1L, n = 1)
```

## Arguments

- x:

  A numeric vector (the time series).

- dimension:

  Embedding dimension \\m\\ (default 2).

- r:

  Tolerance as a fraction of the signal's standard deviation; the
  distance tolerance is `r * sd(x)` (default 0.2, matching NeuroKit2's
  `tolerance = "sd"`).

- delay:

  Time delay (lag) in samples (default 1).

- n:

  Fuzzy power (exponent) of the membership function (default 1).

## Value

A single numeric value, the fuzzy entropy.

## Details

This implementation reproduces `NeuroKit2`'s `entropy_fuzzy` bit-for-bit
(the default `tolerance = "sd"` there corresponds to `r = 0.2`).

## References

Chen, W., Zhuang, J., Yu, W., & Wang, Z. (2009). "Measuring complexity
using FuzzyEn, ApEn, and SampEn." *Medical Engineering & Physics*,
31(1), 61–68.
[doi:10.1016/j.medengphy.2008.04.005](https://doi.org/10.1016/j.medengphy.2008.04.005)

## See also

[`dispersionEntropy`](https://x-biosignal.github.io/PhysioEEG/reference/dispersionEntropy.md),
[`svdEntropy`](https://x-biosignal.github.io/PhysioEEG/reference/svdEntropy.md),
[`renyiEntropy`](https://x-biosignal.github.io/PhysioEEG/reference/renyiEntropy.md),
[`tsallisEntropy`](https://x-biosignal.github.io/PhysioEEG/reference/tsallisEntropy.md).

## Examples

``` r
set.seed(1)
fuzzyEntropy(rnorm(300), dimension = 2)   # fuzzy generalization of sample entropy
#> [1] 1.647206
```

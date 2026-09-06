# Increment Entropy

Computes the increment entropy (IncrEn) of a one-dimensional signal (Liu
et al., 2016) – a symbolic complexity measure defined on the signal's
*increments* (successive differences) rather than its amplitudes. Each
increment in an embedding vector is mapped to a word letter combining
its *sign* (\\-, 0, +\\) with a *magnitude class* in \\\\0, 1, \dots,
q\\\\ obtained by quantizing \\\|{\rm increment}\|\\ relative to the
vector's own standard deviation; IncrEn is the Shannon entropy of the
resulting word distribution, normalized by \\m - 1\\. Because it encodes
both the direction and the graded size of changes, it is sensitive to
dynamics that amplitude-based entropies miss, and is a distinct
mechanism from the amplitude-histogram
([`renyiEntropy`](https://x-biosignal.github.io/PhysioEEG/reference/renyiEntropy.md),
[`tsallisEntropy`](https://x-biosignal.github.io/PhysioEEG/reference/tsallisEntropy.md)),
NCDF-pattern
([`dispersionEntropy`](https://x-biosignal.github.io/PhysioEEG/reference/dispersionEntropy.md))
and template-matching
([`fuzzyEntropy`](https://x-biosignal.github.io/PhysioEEG/reference/fuzzyEntropy.md))
entropies.

## Usage

``` r
incrementEntropy(x, dimension = 2L, q = 4L)
```

## Arguments

- x:

  A numeric vector (the time series).

- dimension:

  Embedding dimension \\m\\ of the increment words (default 2).

- q:

  Number of magnitude quantization levels (default 4); each increment is
  binned into one of `q + 1` magnitude classes by its size relative to
  the embedding vector's standard deviation.

## Value

A single numeric value, the increment entropy.

## Details

This implementation reproduces `NeuroKit2`'s `entropy_increment`
bit-for-bit.

## References

Liu, X., Jiang, A., Xu, N., & Xue, J. (2016). "Increment Entropy as a
Measure of Complexity for Time Series." *Entropy*, 18(1), 22.
[doi:10.3390/e18010022](https://doi.org/10.3390/e18010022)

## See also

[`dispersionEntropy`](https://x-biosignal.github.io/PhysioEEG/reference/dispersionEntropy.md),
[`fuzzyEntropy`](https://x-biosignal.github.io/PhysioEEG/reference/fuzzyEntropy.md),
[`renyiEntropy`](https://x-biosignal.github.io/PhysioEEG/reference/renyiEntropy.md),
[`tsallisEntropy`](https://x-biosignal.github.io/PhysioEEG/reference/tsallisEntropy.md).

## Examples

``` r
set.seed(1)
incrementEntropy(rnorm(1000), dimension = 2, q = 4)
#> [1] 4.588031
```

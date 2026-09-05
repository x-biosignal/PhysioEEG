# Singular Value Decomposition (SVD) Entropy

Computes the SVD entropy of a one-dimensional signal (Roberts et al.,
1999) – a measure of the signal's dimensionality / the number of
eigenvectors needed to explain it. The signal is time-delay embedded
into an `order`-dimensional matrix; the normalized singular values of
that matrix are treated as a probability distribution, and their Shannon
entropy is the SVD entropy. Low SVD entropy indicates a low-dimensional
(structured, e.g. oscillatory) signal; high SVD entropy indicates a
high-dimensional (complex or noise-like) signal.

## Usage

``` r
svdEntropy(x, order = 3L, delay = 1L, normalize = FALSE)
```

## Arguments

- x:

  A numeric vector (the time series).

- order:

  Embedding dimension (default 3); the number of singular values.

- delay:

  Embedding delay / lag in samples (default 1).

- normalize:

  If `TRUE`, divide by `log2(order)` to bound the result in \\\[0,
  1\]\\; otherwise return the entropy in bits (default `FALSE`).

## Value

A single numeric value, the SVD entropy (in bits, or normalized).

## References

Roberts, S. J., Penny, W., & Rezek, I. (1999). "Temporal and spatial
complexity measures for electroencephalogram based brain-computer
interfacing." *Medical & Biological Engineering & Computing*, 37(1),
93–98. [doi:10.1007/BF02513272](https://doi.org/10.1007/BF02513272)

## See also

[`eegComplexity`](https://x-biosignal.github.io/PhysioEEG/reference/eegComplexity.md)
for the multi-measure complexity wrapper (permutation entropy,
Lempel-Ziv, Hjorth, DFA, ...).

## Examples

``` r
set.seed(1)
svdEntropy(sin(seq(0, 20 * pi, length.out = 500)))  # low (structured)
#> [1] 0.44732
svdEntropy(rnorm(500))                               # high (noise-like)
#> [1] 1.584797
```

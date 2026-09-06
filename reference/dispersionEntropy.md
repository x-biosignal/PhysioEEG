# Dispersion Entropy

Computes the dispersion entropy (DispEn) of a one-dimensional signal
(Rostaghi & Azami, 2016) – a fast, robust symbolic-dynamics complexity
measure. Unlike the amplitude-histogram entropies
([`renyiEntropy`](https://x-biosignal.github.io/PhysioEEG/reference/renyiEntropy.md),
[`tsallisEntropy`](https://x-biosignal.github.io/PhysioEEG/reference/tsallisEntropy.md))
and the ordinal permutation entropy, DispEn maps the signal to `c`
amplitude classes through the normal cumulative distribution function
(NCDF), forms embedding vectors of length `dimension` (the *dispersion
patterns*), and takes the Shannon entropy of the pattern distribution,
normalized by \\\ln(c^{m})\\. The NCDF mapping makes it insensitive to
outliers and much cheaper than sample entropy, with no
threshold-crossing boundary effects. Set `reverse = TRUE` to return the
reverse dispersion entropy (RDEn), the mean-square deviation of the
pattern distribution from uniformity (larger = more regular / less
complex).

## Usage

``` r
dispersionEntropy(x, c = 6L, dimension = 3L, delay = 1L, reverse = FALSE)
```

## Arguments

- x:

  A numeric vector (the time series).

- c:

  Number of amplitude classes / symbols (default 6; Rostaghi & Azami
  recommend 4–8).

- dimension:

  Embedding dimension \\m\\ of the dispersion patterns (default 3).

- delay:

  Time delay (lag) in samples (default 1).

- reverse:

  If `TRUE`, return the reverse dispersion entropy (RDEn) instead of
  DispEn (default `FALSE`).

## Value

A single numeric value: the dispersion entropy (or, if `reverse = TRUE`,
the reverse dispersion entropy).

## Details

This implementation reproduces `NeuroKit2`'s `entropy_dispersion`
bit-for-bit, including its normalization convention (a base-2 Shannon
term divided by the natural-log \\\ln(c^{m})\\, so DispEn ranges in
\\\[0, 1/\ln 2\]\\).

## References

Rostaghi, M., & Azami, H. (2016). "Dispersion entropy: A measure for
time-series analysis." *IEEE Signal Processing Letters*, 23(5), 610–614.
[doi:10.1109/LSP.2016.2542881](https://doi.org/10.1109/LSP.2016.2542881)

## See also

[`renyiEntropy`](https://x-biosignal.github.io/PhysioEEG/reference/renyiEntropy.md),
[`tsallisEntropy`](https://x-biosignal.github.io/PhysioEEG/reference/tsallisEntropy.md),
[`svdEntropy`](https://x-biosignal.github.io/PhysioEEG/reference/svdEntropy.md),
[`petrosianFD`](https://x-biosignal.github.io/PhysioEEG/reference/petrosianFD.md).

## Examples

``` r
set.seed(1)
dispersionEntropy(rnorm(1000))                 # DispEn (c = 6, m = 3)
#> [1] 1.414845
dispersionEntropy(rnorm(1000), reverse = TRUE) # reverse dispersion entropy
#> [1] 0.0007654149
```

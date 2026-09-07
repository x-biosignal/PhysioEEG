# Multiscale Entropy

Computes the multiscale entropy (MSE; Costa, Goldberger & Peng, 2002) of
a one-dimensional signal: the sample entropy of the signal after
coarse-graining it at a series of time scales. At scale \\s\\ the signal
is partitioned into non-overlapping windows of length \\s\\, each window
is replaced by its mean, and the sample entropy (embedding `m`,
tolerance `r` times the SD of the ORIGINAL signal, held fixed across
scales) of that coarse-grained series is the MSE at scale \\s\\. A
signal whose MSE stays high or rises at coarse scales carries structure
across time scales; white noise, by contrast, has an MSE that falls
monotonically from scale 1.

## Usage

``` r
multiscaleEntropy(x, m = 2L, r = 0.2, scales = 1:8)
```

## Arguments

- x:

  A numeric vector (the time series).

- m:

  Embedding dimension for the sample entropy (default 2).

- r:

  Tolerance as a fraction of the ORIGINAL signal's SD (default 0.2),
  held fixed across all scales (Costa's convention).

- scales:

  Integer vector of time scales (default `1:8`).

## Value

A named numeric vector of sample entropies, one per scale (`scale1`,
`scale2`, ...). `NA` for a scale whose coarse-grained series is too
short or yields no template matches.

## References

Costa, M., Goldberger, A. L., & Peng, C.-K. (2002). "Multiscale entropy
analysis of complex physiologic time series." *Physical Review Letters*,
89(6), 068102.
[doi:10.1103/PhysRevLett.89.068102](https://doi.org/10.1103/PhysRevLett.89.068102)

## See also

[`eegComplexity`](https://x-biosignal.github.io/PhysioEEG/reference/eegComplexity.md)
for the multi-measure wrapper; PhysioMoCap's `sampleEntropy` for the
single-scale measure.

## Examples

``` r
set.seed(1); multiscaleEntropy(rnorm(1000), scales = 1:5)  # noise: falls with scale
#>   scale1   scale2   scale3   scale4   scale5 
#> 2.157936 1.832660 1.524710 1.409100 1.163269 
```

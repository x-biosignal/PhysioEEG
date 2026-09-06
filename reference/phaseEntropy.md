# Phase Entropy

Computes the phase entropy (PhasEn) of a one-dimensional signal (Rohila
& Sharma, 2019) from its Second-Order Difference Plot (SODP). Each point
plots the first difference \\x(n+\tau)-x(n)\\ against the second
difference \\x(n+2\tau)-x(n+\tau)\\; the polar angle of each point is
binned into \\k\\ equal sectors of \\\[0, 2\pi)\\, and PhasEn is the
(normalized) Shannon entropy of the angle-weighted sector distribution.
It summarizes how the signal's local rate-of-change and its
change-of-rate co-vary – a phase-plane geometry distinct from the
amplitude, template-matching and symbolic (dispersion/increment/slope)
entropies.

## Usage

``` r
phaseEntropy(x, delay = 1L, k = 4L)
```

## Arguments

- x:

  A numeric vector (the time series).

- delay:

  Time delay (lag) in samples for the differences (default 1).

- k:

  Number of angular sectors of the SODP (default 4).

## Value

A single numeric value, the phase entropy.

## Details

This implementation reproduces `NeuroKit2`'s `entropy_phase`
bit-for-bit, including its normalization (a base-2 Shannon term divided
by the natural-log \\\ln k\\, so PhasEn ranges in \\\[0, 1/\ln 2\]\\).

## References

Rohila, A., & Sharma, A. (2019). "Phase entropy: a new complexity
measure for heart rate variability." *Physiological Measurement*,
40(10), 105006.
[doi:10.1088/1361-6579/ab499e](https://doi.org/10.1088/1361-6579/ab499e)

## See also

[`slopeEntropy`](https://x-biosignal.github.io/PhysioEEG/reference/slopeEntropy.md),
[`incrementEntropy`](https://x-biosignal.github.io/PhysioEEG/reference/incrementEntropy.md),
[`dispersionEntropy`](https://x-biosignal.github.io/PhysioEEG/reference/dispersionEntropy.md).

## Examples

``` r
set.seed(1)
phaseEntropy(rnorm(1000), k = 4)
#> [1] 1.154952
```

# Katz Fractal Dimension

Computes Katz's fractal dimension (KFD) of a one-dimensional signal
(Katz, 1988) – a waveform-complexity measure derived from the total
length of the signal curve relative to its extent. Let \\L\\ be the sum
of the absolute successive differences (the curve length), \\a\\ their
mean, and \\d\\ the maximum absolute distance from the first sample;
then \\\mathrm{KFD} = \log\_{10}(L/a) / \log\_{10}(d/a)\\. A smoother,
more regular signal has a lower KFD; a more convoluted (noise-like)
waveform a higher one.

## Usage

``` r
katzFD(x)
```

## Arguments

- x:

  A numeric vector (the time series).

## Value

A single numeric value, Katz's fractal dimension.

## Details

This uses the amplitude-only convention adopted by the de-facto
reference implementations (Python `antropy`, `NeuroKit2`, `pyeeg`,
`mne-features`), in which both the successive distances and the extent
are measured in the signal's amplitude units alone; it reproduces
`antropy.katz_fd` and `NeuroKit2`'s `fractal_katz` bit-for-bit. (Katz's
original definition uses the Euclidean distance in the
sample-index/amplitude plane, which mixes the two axes' units, is
scale-dependent, and is not comparable across tools; the amplitude-only
form is what the field uses and what makes KFD reproducible.)

## References

Katz, M. J. (1988). "Fractals and the analysis of waveforms." *Computers
in Biology and Medicine*, 18(3), 145–156.
[doi:10.1016/0010-4825(88)90041-8](https://doi.org/10.1016/0010-4825%2888%2990041-8)

## See also

[`petrosianFD`](https://x-biosignal.github.io/PhysioEEG/reference/petrosianFD.md)
for the derivative-sign-change fractal view,
[`svdEntropy`](https://x-biosignal.github.io/PhysioEEG/reference/svdEntropy.md),
[`eegComplexity`](https://x-biosignal.github.io/PhysioEEG/reference/eegComplexity.md)
for the multi-measure wrapper.

## Examples

``` r
katzFD(sin(seq(0, 20 * pi, length.out = 500)))  # low (smooth oscillation)
#> [1] 2.461042
set.seed(1); katzFD(rnorm(500))                 # higher (noise-like)
#> [1] 4.653128
```

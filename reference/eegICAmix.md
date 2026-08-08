# Access ICA Metadata

Returns the ICA results stored in `metadata(x)$ica`, including the
mixing matrix, unmixing matrix, channel means, and whitening matrix.

## Usage

``` r
eegICAmix(x)
```

## Arguments

- x:

  A PhysioExperiment object with ICA results.

## Value

A list containing: `mixing` (mixing matrix A, n_channels x
n_components), `unmixing` (unmixing matrix, n_components x n_channels),
`mean` (channel mean vector), `whiten` (whitening matrix), and `method`
(character, algorithm used).

## References

Hyvarinen, A., & Oja, E. (2000). Independent component analysis:
algorithms and applications. Neural Networks, 13(4-5), 411-430.

Bell, A. J., & Sejnowski, T. J. (1995). An information-maximization
approach to blind separation and blind deconvolution. Neural
Computation, 7(6), 1129-1159.

## See also

[`eegICA()`](https://x-biosignal.github.io/PhysioEEG/reference/eegICA.md),
[`eegICAremove()`](https://x-biosignal.github.io/PhysioEEG/reference/eegICAremove.md),
[`eegICAdetect()`](https://x-biosignal.github.io/PhysioEEG/reference/eegICAdetect.md)

## Examples

``` r
if (FALSE) { # \dontrun{
pe <- make_eeg(n_time = 5000, sr = 500)
pe <- eegICA(pe, n_components = 4, method = "fastica")
ica_info <- eegICAmix(pe)
} # }
```

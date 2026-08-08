# Detect Artifact ICA Components

Automatically identifies artifact components using one of four methods:
correlation with frontal channels, kurtosis, spatial weight pattern, or
the ICLabel-style seven-class classifier
([`eegICLabel()`](https://x-biosignal.github.io/PhysioEEG/reference/eegICLabel.md)).
The `"iclabel"` method delegates to
[`eegICLabel()`](https://x-biosignal.github.io/PhysioEEG/reference/eegICLabel.md)
and collapses its result to the two-class contract (any non-brain argmax
is an `"artifact"`), with the score being one minus the brain
probability.

## Usage

``` r
eegICAdetect(
  x,
  method = c("correlation", "kurtosis", "spatial", "iclabel"),
  threshold = 0.3,
  ica_assay = "ica_components"
)
```

## Arguments

- x:

  A PhysioExperiment object with ICA results (from `eegICA`).

- method:

  Detection method: `"correlation"` (frontal channel correlation),
  `"kurtosis"` (excess kurtosis), `"spatial"` (spatial weight pattern),
  or `"iclabel"` (seven-class classifier).

- threshold:

  Threshold for artifact detection. For `"correlation"`, absolute
  correlation \> threshold marks artifact (default: 0.3).

- ica_assay:

  Assay name containing ICA activations (default: `"ica"`).

## Value

A data.frame with columns: `component` (integer index), `type`
(`"artifact"` or `"neural"`), `method` (detection method used), and
`score` (numeric detection score).

## References

Hyvarinen, A., & Oja, E. (2000). Independent component analysis:
algorithms and applications. Neural Networks, 13(4-5), 411-430.

Bell, A. J., & Sejnowski, T. J. (1995). An information-maximization
approach to blind separation and blind deconvolution. Neural
Computation, 7(6), 1129-1159.

## See also

[`eegICA()`](https://x-biosignal.github.io/PhysioEEG/reference/eegICA.md),
[`eegICAremove()`](https://x-biosignal.github.io/PhysioEEG/reference/eegICAremove.md),
[`eegICAmix()`](https://x-biosignal.github.io/PhysioEEG/reference/eegICAmix.md)

## Examples

``` r
if (FALSE) { # \dontrun{
pe <- make_eeg(n_time = 5000, sr = 500)
pe <- eegICA(pe, n_components = 4, method = "fastica")
artifacts <- eegICAdetect(pe, method = "kurtosis")
} # }
```

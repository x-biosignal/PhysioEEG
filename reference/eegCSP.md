# Common Spatial Pattern (CSP) Analysis

Computes Common Spatial Pattern filters for two-class EEG
discrimination. CSP maximizes the variance ratio between two conditions,
making it a standard spatial filtering technique for motor imagery BCI.

## Usage

``` r
eegCSP(x, labels, n_filters = 3, assay_name = NULL, output_assay = "csp")
```

## Arguments

- x:

  A PhysioExperiment object with epoched (3D) EEG data (time x channels
  x trials).

- labels:

  Character or factor vector of class labels, one per trial. Must
  contain exactly two unique classes.

- n_filters:

  Number of CSP filter pairs to retain (default: 3). The total number of
  spatial filters will be `2 * n_filters`.

- assay_name:

  Input assay name (default: first assay).

- output_assay:

  Output assay name for CSP features (default: `"csp"`).

## Value

Modified PhysioExperiment with CSP log-variance features in
`output_assay` (a matrix of trials x `2 * n_filters`) and CSP filter
information stored in `metadata(x)$csp` as a list containing `filters`
(spatial filter matrix), `eigenvalues` (selected eigenvalues), and
`classes` (unique class labels).

## References

Blankertz, B., et al. (2008). Optimizing spatial filters for robust EEG
single-trial analysis. IEEE Signal Processing Magazine, 25(1), 41-56.

## See also

[`eegBCIfeatures()`](https://x-biosignal.github.io/PhysioEEG/reference/eegBCIfeatures.md),
[`eegBCIclassify()`](https://x-biosignal.github.io/PhysioEEG/reference/eegBCIclassify.md),
[`eegMotorImagery()`](https://x-biosignal.github.io/PhysioEEG/reference/eegMotorImagery.md)

## Examples

``` r
if (FALSE) { # \dontrun{
pe <- make_eeg_bci(n_trials = 30, n_channels = 8, sr = 256)
labels <- metadata(pe)$labels
result <- eegCSP(pe, labels = labels, n_filters = 3)
csp_features <- SummarizedExperiment::assay(result, "csp")
} # }
```

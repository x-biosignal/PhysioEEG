# BCI Feature Extraction

Extracts features from epoched EEG data for Brain-Computer Interface
classification. Supports band power, CSP, and Riemannian geometry
methods.

## Usage

``` r
eegBCIfeatures(
  x,
  method = c("bandpower", "csp", "riemannian"),
  labels = NULL,
  bands = NULL,
  assay_name = NULL
)
```

## Arguments

- x:

  A PhysioExperiment object with epoched (3D) EEG data (time x channels
  x trials).

- method:

  Feature extraction method: `"bandpower"` (log band power), `"csp"`
  (Common Spatial Pattern log-variance), or `"riemannian"` (tangent
  space projection of covariance matrices).

- labels:

  Character or factor vector of class labels (required for `"csp"`
  method). One label per trial.

- bands:

  Named list of frequency bands for `"bandpower"` method. Default:
  `list(mu = c(8, 13), beta = c(13, 30))`.

- assay_name:

  Input assay name (default: first assay).

## Value

A numeric matrix with `n_trials` rows and feature columns. Number of
columns depends on method:

- `"bandpower"`: `n_channels * n_bands`

- `"csp"`: `2 * n_filters` (default: 6)

- `"riemannian"`: `n_channels * (n_channels + 1) / 2`

## References

Blankertz, B., et al. (2008). Optimizing spatial filters for robust EEG
single-trial analysis. IEEE Signal Processing Magazine, 25(1), 41-56.

## See also

[`eegCSP()`](https://x-biosignal.github.io/PhysioEEG/reference/eegCSP.md),
[`eegMotorImagery()`](https://x-biosignal.github.io/PhysioEEG/reference/eegMotorImagery.md),
[`eegBCIclassify()`](https://x-biosignal.github.io/PhysioEEG/reference/eegBCIclassify.md)

## Examples

``` r
if (FALSE) { # \dontrun{
pe <- make_eeg_bci(n_trials = 20, n_channels = 8, sr = 256)
features <- eegBCIfeatures(pe, method = "bandpower")
} # }
```

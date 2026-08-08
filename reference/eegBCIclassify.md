# BCI Classification

Classifies BCI features using Linear Discriminant Analysis (LDA) or
shrinkage LDA. Implements Fisher's LDA with optional Ledoit-Wolf
shrinkage regularization for robust classification with high-dimensional
or small-sample data. Optionally performs k-fold cross-validation.

## Usage

``` r
eegBCIclassify(
  x,
  features = NULL,
  labels,
  method = c("lda", "shrinkage_lda"),
  cv_folds = NULL,
  assay_name = NULL
)
```

## Arguments

- x:

  A PhysioExperiment object with epoched (3D) EEG data.

- features:

  Optional pre-computed feature matrix (n_trials x n_features). If
  `NULL`, features are extracted using
  [`eegBCIfeatures`](https://x-biosignal.github.io/PhysioEEG/reference/eegBCIfeatures.md)
  with the `"bandpower"` method.

- labels:

  Character or factor vector of class labels, one per trial. Must
  contain exactly two unique classes.

- method:

  Classification method: `"lda"` (Fisher's LDA) or `"shrinkage_lda"`
  (LDA with Ledoit-Wolf shrinkage).

- cv_folds:

  Number of cross-validation folds (default: `NULL`, meaning no
  cross-validation). If set (e.g., 5), performs k-fold CV and reports
  out-of-fold predictions. The CV accuracy is stored as
  `attr(result, "cv_accuracy")`.

- assay_name:

  Input assay name used when extracting features (default: first assay).

## Value

A data.frame with columns: `trial`, `predicted_class`, `confidence`, and
`true_class`. When `cv_folds` is not `NULL`, predictions are out-of-fold
and `attr(result, "cv_accuracy")` contains the cross-validated accuracy.
The trained LDA model (on all data) is stored in
`metadata(x)$bci_model`, containing `weights`, `threshold`, `classes`,
`method`, and `class_means`.

## References

Blankertz, B., et al. (2008). Optimizing spatial filters for robust EEG
single-trial analysis. IEEE Signal Processing Magazine, 25(1), 41-56.

## See also

[`eegCSP()`](https://x-biosignal.github.io/PhysioEEG/reference/eegCSP.md),
[`eegBCIfeatures()`](https://x-biosignal.github.io/PhysioEEG/reference/eegBCIfeatures.md),
[`eegMotorImagery()`](https://x-biosignal.github.io/PhysioEEG/reference/eegMotorImagery.md)

## Examples

``` r
if (FALSE) { # \dontrun{
pe <- make_eeg_bci(n_trials = 20, n_channels = 8, sr = 256)
labels <- metadata(pe)$labels
features <- eegBCIfeatures(pe, method = "bandpower")
result <- eegBCIclassify(pe, features = features, labels = labels, method = "lda")

# With 5-fold cross-validation
result_cv <- eegBCIclassify(pe, features = features, labels = labels,
                            method = "lda", cv_folds = 5)
attr(result_cv, "cv_accuracy")
} # }
```

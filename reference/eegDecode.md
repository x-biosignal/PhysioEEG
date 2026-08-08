# Leakage-free MVPA Decoding of Epoched EEG

Decodes trial class labels from epoched EEG with a strictly fold-safe
pipeline: the spatial filter (CSP), feature extraction, feature scaling
and classifier are all estimated on the training fold only and applied,
frozen, to the test fold (Varoquaux et al., 2017). This avoids the
label/data leakage of fitting CSP or scaling on the full dataset before
cross-validation. Folds may respect a run/block grouping to avoid
temporally-contiguous-epoch leakage.

## Usage

``` r
eegDecode(
  x,
  labels = NULL,
  pipeline = c("csp+lda", "bandpower+lda", "riemannian"),
  cv = c("stratified-kfold", "leave-one-run-out"),
  n_folds = 5,
  groups = NULL,
  n_permutations = 200,
  inner_cv = NULL,
  assay_name = NULL,
  n_csp = 3,
  bands = NULL,
  seed = NULL
)
```

## Arguments

- x:

  A PhysioExperiment with an epoched (3D: time x channels x trials)
  assay.

- labels:

  Class labels, one per trial. If `NULL`, taken from
  `metadata(x)$labels`.

- pipeline:

  Decoding pipeline: `"csp+lda"` (common spatial patterns then LDA),
  `"bandpower+lda"` (per-channel band-power then LDA), or `"riemannian"`
  (tangent-space covariance features then LDA).

- cv:

  Cross-validation scheme: `"stratified-kfold"` or `"leave-one-run-out"`
  (the latter requires `groups`).

- n_folds:

  Number of folds for k-fold CV (default: 5).

- groups:

  Optional run/block grouping (one per trial); folds keep each group
  intact so contiguous epochs are never split across train/test.

- n_permutations:

  Number of label permutations for the null distribution and p-value
  (default: 200; 0 to skip).

- inner_cv:

  Optional number of inner CV folds for nested hyperparameter selection
  (shrinkage LDA vs plain LDA). `NULL` disables nesting.

- assay_name:

  Input assay (default: `defaultAssay(x)`).

- n_csp:

  Number of CSP filter pairs for `"csp+lda"` (default: 3).

- bands:

  Named list of frequency bands for `"bandpower+lda"` (default:
  `list(mu = c(8, 13), beta = c(13, 30))`).

- seed:

  Optional RNG seed for reproducible fold assignment/permutations.

## Value

A list of class `"eeg_decode"` with:

- accuracy:

  pooled out-of-fold accuracy.

- accuracy_ci:

  95 percent Wilson confidence interval for the accuracy.

- fold_accuracy:

  per-fold accuracy.

- auc:

  binary area under the ROC curve.

- confusion:

  confusion matrix (true x predicted).

- permutation:

  list with `p_value`, `null_accuracy`, `n_permutations`.

- predictions:

  a `data.frame` (trial, fold, true, predicted, score).

- csp_filters:

  per-fold CSP filters (`"csp+lda"` only).

## References

Varoquaux, G., et al. (2017). "Assessing and tuning brain decoders:
cross-validation, caveats, and guidelines." *NeuroImage*, 145, 166-179.
[doi:10.1016/j.neuroimage.2016.10.038](https://doi.org/10.1016/j.neuroimage.2016.10.038)

Lemm, S., et al. (2011). "Introduction to machine learning for brain
imaging." *NeuroImage*, 56(2), 387-399.
[doi:10.1016/j.neuroimage.2010.11.004](https://doi.org/10.1016/j.neuroimage.2010.11.004)

## See also

[`eegCSP`](https://x-biosignal.github.io/PhysioEEG/reference/eegCSP.md),
[`eegBCIfeatures`](https://x-biosignal.github.io/PhysioEEG/reference/eegBCIfeatures.md).
The older
[`eegBCIclassify`](https://x-biosignal.github.io/PhysioEEG/reference/eegBCIclassify.md)
cross-validation path is leakage-prone and deprecated in favour of
`eegDecode`.

## Examples

``` r
pe <- make_eeg_bci(n_trials = 20, n_channels = 8, sr = 128, trial_sec = 3)
res <- eegDecode(pe, pipeline = "csp+lda", n_permutations = 0, seed = 1)
res$accuracy
#> [1] 1
```

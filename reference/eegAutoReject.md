# Automated epoch rejection and interpolation (autoreject)

Implements the autoreject algorithm (Jas et al. 2017) for epoched EEG.
Per-channel peak-to-peak rejection thresholds are chosen by
cross-validation, then each epoch is kept, repaired by interpolating its
worst bad channels, or dropped, according to two hyperparameters:
`consensus` (the fraction of bad channels above which an epoch is
dropped) and `n_interpolate` (the number of worst bad channels to
interpolate in a repaired epoch). Either hyperparameter is chosen by
cross-validation when not supplied. Bad channels are repaired with
spherical-spline interpolation (the same weights used by
[`eegInterpolate()`](https://x-biosignal.github.io/PhysioEEG/reference/eegInterpolate.md)),
so electrode positions must be present (run
[`eegMontage()`](https://x-biosignal.github.io/PhysioEEG/reference/eegMontage.md)
first).

## Usage

``` r
eegAutoReject(
  x,
  consensus = NULL,
  n_interpolate = NULL,
  cv_folds = 5L,
  thresh_range = NULL,
  n_thresh = 20L,
  assay_name = "epoched",
  output_assay = "autoreject_clean"
)
```

## Arguments

- x:

  A PhysioExperiment with 3D epoched data (time x channels x epochs) in
  `metadata(x)[[assay_name]]` (from
  [`eegEpoch()`](https://x-biosignal.github.io/PhysioEEG/reference/eegEpoch.md))
  or in a 3D assay of that name (from
  [`make_eeg_erp()`](https://x-biosignal.github.io/PhysioEEG/reference/make_eeg_erp.md)).

- consensus:

  Fraction of bad channels (from 0 to 1) above which an epoch is
  dropped. `NULL` selects it by cross-validation.

- n_interpolate:

  Number of worst bad channels to interpolate per repaired epoch. `NULL`
  selects it by cross-validation.

- cv_folds:

  Number of cross-validation folds (default: 5).

- thresh_range:

  Optional candidate peak-to-peak threshold grid: either a length-2
  `c(min, max)` span (expanded to `n_thresh` values) or an explicit
  vector of candidate thresholds. `NULL` derives the span from the data.

- n_thresh:

  Number of grid thresholds when `thresh_range` is `NULL` or a span
  (default: 20).

- assay_name:

  Name of the 3D epoched data in metadata or assays (default:
  `"epoched"`).

- output_assay:

  Name for the cleaned 3D data stored in metadata (default:
  `"autoreject_clean"`).

## Value

The PhysioExperiment with cleaned 3D epoched data in
`metadata(x)[[output_assay]]` (dropped epochs removed, repaired epochs
interpolated) and a rejection log in `metadata(x)$autoreject`: the
per-channel `thresholds`, chosen `consensus` and `n_interpolate`, the
epoch-by-channel `bad_matrix`, the `interpolated` channel-epochs, and
the `dropped_epochs`. The step is recorded in provenance.

## Details

The cross-validation objective, at every stage, is the root-mean-square
error between the mean of the good/cleaned training epochs and the
median of the held-out test epochs - the criterion introduced by
autoreject.

## References

Jas, M., Engemann, D. A., Bekhti, Y., Raimondo, F., & Gramfort, A.
(2017). Autoreject: Automated artifact rejection for MEG and EEG data.
NeuroImage, 159, 417-429.

## See also

[`eegArtifactReject()`](https://x-biosignal.github.io/PhysioEEG/reference/eegArtifactReject.md),
[`eegInterpolate()`](https://x-biosignal.github.io/PhysioEEG/reference/eegInterpolate.md),
[`eegEpoch()`](https://x-biosignal.github.io/PhysioEEG/reference/eegEpoch.md),
[`eegMontage()`](https://x-biosignal.github.io/PhysioEEG/reference/eegMontage.md)

## Examples

``` r
if (FALSE) { # \dontrun{
pe <- make_eeg_erp(n_epochs = 40, n_channels = 19, sr = 250)
pe <- eegMontage(pe, system = "10-20")
pe <- eegAutoReject(pe, assay_name = "raw")
dim(S4Vectors::metadata(pe)$autoreject$bad_matrix)
} # }
```

# Fit an equivalent current dipole (ECD) to an EEG topography

Localizes the single equivalent current dipole that best explains a
scalp topography, by nonlinearly optimizing the dipole position (a grid
search inside the head sphere refined with Nelder-Mead) while the dipole
moment is the closed-form least-squares fit at each candidate position.
Uses the same forward physics as
[`eegForwardModel()`](https://x-biosignal.github.io/PhysioEEG/reference/eegForwardModel.md)
and the electrode positions from `colData` / the montage. This is the
classic focal-source complement to the distributed inverses in
[`eegSourceEstimate()`](https://x-biosignal.github.io/PhysioEEG/reference/eegSourceEstimate.md).

## Usage

``` r
eegDipoleFit(
  x,
  time = NULL,
  method = c("spherical", "bem_simplified"),
  n_grid = 8L,
  assay_name = NULL
)
```

## Arguments

- x:

  A `PhysioExperiment`.

- time:

  Sample index to fit (a single index), a length-2 range to average
  over, or `NULL` (default) to fit the peak global-field-power sample.

- method:

  Forward physics: `"spherical"` (single-sphere) or `"bem_simplified"`
  (Berg-Scherg 3-shell).

- n_grid:

  Grid resolution per axis for the initial search (default 8).

- assay_name:

  Assay to use (default: the object's default assay).

## Value

An `eeg_dipole_fit` object: `position` (x, y, z in normalized head
radius), `moment`, `orientation` (unit), `amplitude`, `gof` (goodness of
fit, variance explained), `residual`, and the fitted/observed
topographies.

## References

Scherg (1990); Mosher, Lewis & Leahy (1992).

## See also

[`eegSourceEstimate()`](https://x-biosignal.github.io/PhysioEEG/reference/eegSourceEstimate.md),
[`eegForwardModel()`](https://x-biosignal.github.io/PhysioEEG/reference/eegForwardModel.md)

## Examples

``` r
pe <- make_eeg(n_time = 100, n_channels = 19, sr = 100)
fit <- eegDipoleFit(pe)
fit$gof
#> [1] 0.3307643
```

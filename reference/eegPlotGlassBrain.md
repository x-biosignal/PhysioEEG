# Plot Source Localization as Multi-View Glass-Brain Projections

Creates deterministic sagittal, axial, and coronal maximum-intensity
projections from source amplitudes with real three-dimensional
coordinates. The bundled outline is schematic and is not an MRI, patient
anatomy, or a validation of the inverse solution.

## Usage

``` r
eegPlotGlassBrain(
  x,
  source_data = NULL,
  views = c("sagittal", "axial", "coronal"),
  threshold_pct = 90
)
```

## Arguments

- x:

  A PhysioExperiment object.

- source_data:

  Optional explicit source data. Supply a data frame with finite `x`,
  `y`, `z`, and `amplitude` columns, or a structured list with
  `positions` and declared amplitude/matrix reduction fields. If `NULL`,
  uses coordinate-bearing metadata from a new
  [`eegSourceEstimate()`](https://x-biosignal.github.io/PhysioEEG/reference/eegSourceEstimate.md)
  or
  [`eegBeamformer()`](https://x-biosignal.github.io/PhysioEEG/reference/eegBeamformer.md)
  result.

- views:

  Non-empty unique character vector drawn exactly from `"sagittal"`,
  `"axial"`, and `"coronal"`. Request order is preserved.

- threshold_pct:

  Finite percentile in `[0, 100]`. The threshold is the type-8 quantile
  of absolute amplitude over the complete source set; ties are retained.

## Value

A patchwork object for normal projections. When fewer than three rows
survive thresholding, returns the first requested view as an
[`eegPlotSource()`](https://x-biosignal.github.io/PhysioEEG/reference/eegPlotSource.md)
ggplot fallback without applying a second threshold. Immutable
`"glassbrain_data"` attributes retain the resolved sources, panel data,
display transform, and threshold settings.

## Details

Projection axes are sagittal `y-z` with depth `x`, axial `x-y` with
depth `z`, and coronal `x-z` with depth `y`. Free-orientation source
estimates are reduced by root-sum-square over orientations at each time
and RMS over time. This non-negative summary is not a signed
instantaneous current. A fixed 80 by 80 display grid retains the
greatest absolute amplitude in each projected bin, with exact ties
resolved by the smallest source ID. The affine display transform comes
from the complete unthresholded coordinate cloud and does not modify
amplitudes. The bundled outline is schematic: the plot is not MRI
registration, patient anatomy, or validation of localization accuracy.
Bare amplitude vectors and two-dimensional coordinates are rejected.

## See also

[`eegForwardModel()`](https://x-biosignal.github.io/PhysioEEG/reference/eegForwardModel.md),
[`eegSourceEstimate()`](https://x-biosignal.github.io/PhysioEEG/reference/eegSourceEstimate.md),
[`eegBeamformer()`](https://x-biosignal.github.io/PhysioEEG/reference/eegBeamformer.md),
[`eegPlotSource()`](https://x-biosignal.github.io/PhysioEEG/reference/eegPlotSource.md)

## Examples

``` r
if (FALSE) { # \dontrun{
pe <- make_eeg(n_time = 500, n_channels = 19, sr = 250)
fm <- eegForwardModel(pe, method = "spherical", n_sources = 50)
localized <- eegSourceEstimate(pe, fm, method = "mne")
eegPlotGlassBrain(localized)
} # }
```

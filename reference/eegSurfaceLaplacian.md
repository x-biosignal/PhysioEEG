# Surface Laplacian (current source density) of EEG

Applies the spherical-spline surface Laplacian (current source density,
CSD) of Perrin et al. (1989). CSD is **reference-free** and acts as a
spatial high-pass filter, deblurring volume conduction so that each
channel reflects the radial current under it rather than the whole-head
average — sharpening focal topographies and attenuating broadly
distributed activity. It needs electrode positions (from `colData`
columns `pos_x/pos_y/pos_z`, else matched to a standard 10-10/10-20
montage by channel label).

## Usage

``` r
eegSurfaceLaplacian(
  x,
  order = 4L,
  lambda = 1e-05,
  n_terms = 50L,
  assay_name = NULL,
  output_assay = "csd"
)
```

## Arguments

- x:

  A `PhysioExperiment`.

- order:

  Spline flexibility `m` (default 4; higher = stiffer).

- lambda:

  Smoothing / regularization added to the spline system (default
  `1e-5`).

- n_terms:

  Number of Legendre terms in the g/h series (default 50).

- assay_name:

  Assay to transform (default: the object's default assay).

- output_assay:

  Name of the assay to store the CSD in (default `"csd"`).

## Value

The `PhysioExperiment` with the CSD stored in `output_assay` (same shape
as the input assay; units are µV/m² up to a scale/head-radius constant).

## References

Perrin, Pernier, Bertrand & Echallier (1989), Electroenceph Clin
Neurophysiol; Kayser & Tenke (2006), Clin Neurophysiol.

## See also

[`eegRereference()`](https://x-biosignal.github.io/PhysioEEG/reference/eegRereference.md),
[`eegPlotTopomap()`](https://x-biosignal.github.io/PhysioEEG/reference/eegPlotTopomap.md),
[`eegConnectivityMatrix()`](https://x-biosignal.github.io/PhysioEEG/reference/eegConnectivityMatrix.md)

## Examples

``` r
pe <- make_eeg(n_time = 500, n_channels = 19, sr = 250)
pe <- eegSurfaceLaplacian(pe)
"csd" %in% SummarizedExperiment::assayNames(pe)
#> [1] TRUE
```

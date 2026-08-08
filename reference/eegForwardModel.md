# Construct EEG Forward Model (Leadfield Matrix)

Constructs a leadfield matrix that maps brain source activity to scalp
electrode potentials. Uses electrode positions from `colData(x)` if
available, otherwise falls back to standard 10-20 system positions on a
unit sphere.

## Usage

``` r
eegForwardModel(
  x,
  method = c("spherical", "bem_simplified", "bem", "sphere_analytic", "nyhead",
    "fsaverage"),
  n_sources = 500,
  assay_name = NULL
)
```

## Arguments

- x:

  A PhysioExperiment object with EEG data.

- method:

  Forward model method: `"spherical"` (current dipole in infinite
  homogeneous medium with conductivity 0.33 S/m), `"bem_simplified"`
  (3-shell Berg and Scherg 1994 approximation), `"sphere_analytic"` (the
  analytic single-sphere leadfield over a structured spherical cortical
  source space), `"bem"` (a 3-shell boundary element method over the
  same source space), or `"nyhead"` / `"fsaverage"` (published head
  models loaded on demand). The last four require the PhysioHeadModels
  package; the `"nyhead"` / `"fsaverage"` methods additionally require
  the corresponding dataset to have been downloaded (see
  [`PhysioHeadModels::fetchNYHead`](https://rdrr.io/pkg/PhysioHeadModels/man/fetchNYHead.html)).

- n_sources:

  Number of dipole sources to distribute inside the head model (default:
  500).

- assay_name:

  Name of the assay to reference for channel count. If `NULL`, the
  default assay is used.

## Value

A list with components:

- leadfield:

  Numeric matrix of dimensions n_electrodes x (n_sources \* 3). Each
  source has 3 orientation columns (x, y, z).

- source_positions:

  Data frame with columns `x`, `y`, `z` for each source location.

- electrode_positions:

  Data frame with columns `label`, `x`, `y`, `z` for each electrode.

- n_sources:

  Integer number of source dipoles.

- source_normals:

  For the realistic methods, an n_sources x 3 matrix of cortical surface
  normals (absent for the legacy methods).

- method:

  The resolved forward-model method.

- coordinate_provenance:

  Source-space origin, units, and anatomical fidelity used by
  coordinate-aware plotting functions.

## References

Pascual-Marqui, R. D. (2002). Standardized low-resolution brain
electromagnetic tomography (sLORETA). Methods and Findings in
Experimental and Clinical Pharmacology, 24(Suppl D), 5-12.

## See also

[`eegSourceEstimate()`](https://x-biosignal.github.io/PhysioEEG/reference/eegSourceEstimate.md),
[`eegBeamformer()`](https://x-biosignal.github.io/PhysioEEG/reference/eegBeamformer.md),
[`eegSourcePower()`](https://x-biosignal.github.io/PhysioEEG/reference/eegSourcePower.md)

## Examples

``` r
if (FALSE) { # \dontrun{
pe <- make_eeg(n_time = 1000, n_channels = 19, sr = 250)
fm <- eegForwardModel(pe, method = "spherical", n_sources = 100)
dim(fm$leadfield)
} # }
```

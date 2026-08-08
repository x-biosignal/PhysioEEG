# EEG Source Estimation

Estimates brain source activity from scalp EEG data using distributed
source imaging methods. Requires a forward model from
[`eegForwardModel`](https://x-biosignal.github.io/PhysioEEG/reference/eegForwardModel.md).

## Usage

``` r
eegSourceEstimate(
  x,
  forward_model,
  method = c("sloreta", "eloreta", "mne"),
  lambda = 0.05,
  assay_name = NULL,
  output_assay = "source"
)
```

## Arguments

- x:

  A PhysioExperiment object with EEG data.

- forward_model:

  A forward model list as returned by
  [`eegForwardModel`](https://x-biosignal.github.io/PhysioEEG/reference/eegForwardModel.md).

- method:

  Source estimation method: `"sloreta"` (standardized low-resolution
  tomography), `"eloreta"` (exact low-resolution tomography), or `"mne"`
  (minimum norm estimate).

- lambda:

  Regularization parameter (default: 0.05). Higher values produce
  smoother solutions.

- assay_name:

  Name of the input assay. If `NULL`, the default assay is used.

- output_assay:

  Name for the output assay containing source estimates (default:
  `"source"`). The metadata names `"source_estimate"`,
  `"beamformer_info"`, and `"source_plot_default"` are reserved.

## Value

Modified PhysioExperiment with source estimates stored in
`output_assay`. The assay is a matrix of dimensions n_time x (n_sources
\* 3). Sets `metadata(x)$source_estimate` with a list containing:
`method`, `lambda`, `n_sources`, `n_source_cols`, `source_positions`,
`output_assay`, `orientation_count`, and coordinate provenance.

## References

Pascual-Marqui, R. D. (2002). Standardized low-resolution brain
electromagnetic tomography (sLORETA). Methods and Findings in
Experimental and Clinical Pharmacology, 24(Suppl D), 5-12.

## See also

[`eegForwardModel()`](https://x-biosignal.github.io/PhysioEEG/reference/eegForwardModel.md),
[`eegBeamformer()`](https://x-biosignal.github.io/PhysioEEG/reference/eegBeamformer.md),
[`eegSourcePower()`](https://x-biosignal.github.io/PhysioEEG/reference/eegSourcePower.md),
[`eegPlotSource()`](https://x-biosignal.github.io/PhysioEEG/reference/eegPlotSource.md),
[`eegPlotGlassBrain()`](https://x-biosignal.github.io/PhysioEEG/reference/eegPlotGlassBrain.md)

## Examples

``` r
if (FALSE) { # \dontrun{
pe <- make_eeg(n_time = 500, n_channels = 19, sr = 250)
fm <- eegForwardModel(pe, method = "spherical", n_sources = 50)
pe <- eegSourceEstimate(pe, fm, method = "mne")
} # }
```

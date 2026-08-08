# EEG Beamformer Source Localization

Applies spatial filtering (beamforming) to localize neural source power.
Linearly Constrained Minimum Variance (LCMV) beamformer operates in the
time domain. Dynamic Imaging of Coherent Sources (DICS) operates in the
frequency domain.

## Usage

``` r
eegBeamformer(
  x,
  forward_model,
  method = c("lcmv", "dics"),
  freq_range = NULL,
  assay_name = NULL,
  output_assay = "beamformer"
)
```

## Arguments

- x:

  A PhysioExperiment object with EEG data.

- forward_model:

  A forward model list as returned by
  [`eegForwardModel`](https://x-biosignal.github.io/PhysioEEG/reference/eegForwardModel.md).

- method:

  Beamformer method: `"lcmv"` (Linearly Constrained Minimum Variance) or
  `"dics"` (Dynamic Imaging of Coherent Sources).

- freq_range:

  Numeric vector of length 2 specifying frequency range in Hz for DICS
  method (e.g., `c(8, 13)` for alpha band). Ignored for LCMV.

- assay_name:

  Name of the input assay. If `NULL`, the default assay is used.

- output_assay:

  Name for the output assay containing beamformer results (default:
  `"beamformer"`). The metadata names `"source_estimate"`,
  `"beamformer_info"`, and `"source_plot_default"` are reserved.

## Value

Modified PhysioExperiment with source power stored in `output_assay` as
a matrix with one column named `"power"` (n_sources rows). Each value
represents the estimated source power at the corresponding dipole
location. Plotting metadata, including source positions, the output
name, method, and coordinate provenance, is stored in
`metadata(x)$beamformer_info`.

## References

Pascual-Marqui, R. D. (2002). Standardized low-resolution brain
electromagnetic tomography (sLORETA). Methods and Findings in
Experimental and Clinical Pharmacology, 24(Suppl D), 5-12.

Van Veen, B. D., et al. (1997). Localization of brain electrical
activity via linearly constrained minimum variance spatial filtering.
IEEE Transactions on Biomedical Engineering, 44(9), 867-880.

## See also

[`eegForwardModel()`](https://x-biosignal.github.io/PhysioEEG/reference/eegForwardModel.md),
[`eegSourceEstimate()`](https://x-biosignal.github.io/PhysioEEG/reference/eegSourceEstimate.md),
[`eegSourcePower()`](https://x-biosignal.github.io/PhysioEEG/reference/eegSourcePower.md),
[`eegPlotGlassBrain()`](https://x-biosignal.github.io/PhysioEEG/reference/eegPlotGlassBrain.md)

## Examples

``` r
if (FALSE) { # \dontrun{
pe <- make_eeg(n_time = 1000, n_channels = 19, sr = 250)
fm <- eegForwardModel(pe, method = "spherical", n_sources = 50)
pe <- eegBeamformer(pe, fm, method = "lcmv")
} # }
```

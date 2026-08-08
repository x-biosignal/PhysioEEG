# SSVEP Frequency Detection

Detects Steady-State Visual Evoked Potentials using canonical
correlation analysis (CCA) or filter-bank CCA (FBCCA). For each
candidate stimulus frequency, a CCA is computed between the EEG data and
sinusoidal reference signals at the frequency and its harmonics.

## Usage

``` r
eegSSVEP(
  x,
  frequencies,
  n_harmonics = 3,
  method = c("cca", "fbcca"),
  assay_name = NULL
)
```

## Arguments

- x:

  A PhysioExperiment object with EEG data. If 3D (time x channels x
  trials), data is averaged across trials before analysis.

- frequencies:

  Numeric vector of target stimulus frequencies in Hz.

- n_harmonics:

  Number of harmonics to include in reference signals (default: 3).

- method:

  Detection method: `"cca"` (canonical correlation analysis) or
  `"fbcca"` (filter-bank CCA).

- assay_name:

  Input assay name (default: first assay).

## Value

A data.frame with columns: `frequency` (numeric target frequency in Hz),
`correlation` (numeric CCA correlation), `snr` (numeric signal-to-noise
ratio), and `predicted_class` (numeric predicted stimulus frequency).

## References

Blankertz, B., et al. (2008). Optimizing spatial filters for robust EEG
single-trial analysis. IEEE Signal Processing Magazine, 25(1), 41-56.

Norcia, A. M., et al. (2015). The steady-state visual evoked potential
in vision research: a review. Journal of Vision, 15(6), 4.

## See also

[`eegCSP()`](https://x-biosignal.github.io/PhysioEEG/reference/eegCSP.md),
[`eegBCIfeatures()`](https://x-biosignal.github.io/PhysioEEG/reference/eegBCIfeatures.md),
[`eegBCIclassify()`](https://x-biosignal.github.io/PhysioEEG/reference/eegBCIclassify.md)

## Examples

``` r
if (FALSE) { # \dontrun{
pe <- make_eeg_bci(n_trials = 10, n_channels = 8, sr = 256)
result <- eegSSVEP(pe, frequencies = c(10, 12, 15), method = "cca")
} # }
```

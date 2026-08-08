# Motor Imagery ERD/ERS Computation

Computes event-related desynchronization (ERD) and event-related
synchronization (ERS) for motor imagery BCI analysis. ERD/ERS is
expressed as a percentage change from baseline power in specified
frequency bands.

## Usage

``` r
eegMotorImagery(
  x,
  bands = NULL,
  baseline_fraction = 0.25,
  assay_name = NULL,
  output_assay = "erd_ers"
)
```

## Arguments

- x:

  A PhysioExperiment object with epoched (3D) EEG data (time x channels
  x trials).

- bands:

  Named list of frequency bands, each a numeric vector `c(low, high)`.
  Default: `list(mu = c(8, 13), beta = c(13, 30))`.

- baseline_fraction:

  Fraction of each trial to use as baseline (default: 0.25, i.e., the
  first 25 percent of the trial).

- assay_name:

  Input assay name (default: first assay).

- output_assay:

  Output assay name (default: `"erd_ers"`).

## Value

Modified PhysioExperiment with ERD/ERS percentage values in
`output_assay`. The output is a matrix of `n_trials` rows x
`n_channels * n_bands` columns. Sets `metadata(x)$erd_ers_bands` with
the band definitions used.

## References

Blankertz, B., et al. (2008). Optimizing spatial filters for robust EEG
single-trial analysis. IEEE Signal Processing Magazine, 25(1), 41-56.

## See also

[`eegCSP()`](https://x-biosignal.github.io/PhysioEEG/reference/eegCSP.md),
[`eegBCIfeatures()`](https://x-biosignal.github.io/PhysioEEG/reference/eegBCIfeatures.md),
[`eegBCIclassify()`](https://x-biosignal.github.io/PhysioEEG/reference/eegBCIclassify.md)

## Examples

``` r
if (FALSE) { # \dontrun{
pe <- make_eeg_bci(n_trials = 20, n_channels = 8, sr = 256)
result <- eegMotorImagery(pe)
erd_data <- SummarizedExperiment::assay(result, "erd_ers")
} # }
```

# Baseline Correct ERP Data

Subtracts the mean amplitude of a pre-stimulus baseline period from each
epoch. This is an essential preprocessing step before ERP measurement.

## Usage

``` r
eegERPbaseline(
  x,
  baseline = c(-200, 0),
  epoch_start = 0,
  assay_name = NULL,
  output_assay = "baseline_corrected"
)
```

## Arguments

- x:

  A PhysioExperiment object with epoched (3D) EEG data.

- baseline:

  Baseline period in milliseconds as `c(start, end)`. For example,
  `c(-200, 0)` for 200ms pre-stimulus baseline. Default is `c(-200, 0)`.

- epoch_start:

  Start time of the epoch in milliseconds relative to stimulus onset.
  Default is 0 (epoch starts at stimulus).

- assay_name:

  Input assay name (default: first assay).

- output_assay:

  Output assay name (default: `"baseline_corrected"`).

## Value

Modified PhysioExperiment with baseline-corrected data in
`output_assay`. Creates a new 3D assay (time x channels x epochs) where
each epoch has had the mean of its baseline period subtracted.

## References

Luck, S. J. (2014). An Introduction to the Event-Related Potential
Technique (2nd ed.). MIT Press.

## See also

[`eegERPdetect()`](https://x-biosignal.github.io/PhysioEEG/reference/eegERPdetect.md),
[`eegERPmeasure()`](https://x-biosignal.github.io/PhysioEEG/reference/eegERPmeasure.md),
[`eegEpoch()`](https://x-biosignal.github.io/PhysioEEG/reference/eegEpoch.md),
[`eegFilter()`](https://x-biosignal.github.io/PhysioEEG/reference/eegFilter.md)

## Examples

``` r
if (FALSE) { # \dontrun{
pe <- make_eeg_erp(n_epochs = 40, sr = 250)
result <- eegERPbaseline(pe, baseline = c(-200, 0), epoch_start = -200)
} # }
```

# ERP Difference Waveform

Computes the difference waveform between two PhysioExperiment objects by
subtracting the assay data of `y` from `x`.

## Usage

``` r
eegERPdifference(x, y, assay_name = NULL, output_assay = "difference")
```

## Arguments

- x:

  A PhysioExperiment object (minuend).

- y:

  A PhysioExperiment object (subtrahend). Must have the same dimensions
  as `x`.

- assay_name:

  Input assay name (default: first assay).

- output_assay:

  Output assay name (default: `"difference"`).

## Value

Modified `x` with difference waveform stored in `output_assay`. The
difference assay has the same dimensions as the input assay.

## References

Luck, S. J. (2014). An Introduction to the Event-Related Potential
Technique (2nd ed.). MIT Press.

## See also

[`eegERPdetect()`](https://x-biosignal.github.io/PhysioEEG/reference/eegERPdetect.md),
[`eegERPmeasure()`](https://x-biosignal.github.io/PhysioEEG/reference/eegERPmeasure.md),
[`eegERPtest()`](https://x-biosignal.github.io/PhysioEEG/reference/eegERPtest.md),
[`eegERPgrandAverage()`](https://x-biosignal.github.io/PhysioEEG/reference/eegERPgrandAverage.md)

## Examples

``` r
if (FALSE) { # \dontrun{
pe1 <- make_eeg_erp(n_epochs = 20, sr = 250)
pe2 <- make_eeg_erp(n_epochs = 20, sr = 250)
result <- eegERPdifference(pe1, pe2)
} # }
```

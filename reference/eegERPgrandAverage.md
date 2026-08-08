# Compute Grand Average ERP

Averages ERP waveforms across multiple PhysioExperiment objects
(participants). Each input should already be averaged across trials.

## Usage

``` r
eegERPgrandAverage(..., assay_name = NULL, output_assay = "grand_average")
```

## Arguments

- ...:

  PhysioExperiment objects to average, or a list of them.

- assay_name:

  Input assay name (default: first assay).

- output_assay:

  Output assay name (default: `"grand_average"`).

## Value

A PhysioExperiment (the first input object) with the grand average
waveform stored in `output_assay`. The grand average assay has the same
dimensions as the input assays.

## References

Luck, S. J. (2014). An Introduction to the Event-Related Potential
Technique (2nd ed.). MIT Press.

## See also

[`eegERPdetect()`](https://x-biosignal.github.io/PhysioEEG/reference/eegERPdetect.md),
[`eegERPmeasure()`](https://x-biosignal.github.io/PhysioEEG/reference/eegERPmeasure.md),
[`eegERPtest()`](https://x-biosignal.github.io/PhysioEEG/reference/eegERPtest.md),
[`eegERPdifference()`](https://x-biosignal.github.io/PhysioEEG/reference/eegERPdifference.md)

## Examples

``` r
if (FALSE) { # \dontrun{
pe1 <- make_eeg_erp(n_epochs = 20, sr = 250)
pe2 <- make_eeg_erp(n_epochs = 20, sr = 250)
result <- eegERPgrandAverage(pe1, pe2)
} # }
```

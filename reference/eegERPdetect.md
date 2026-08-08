# Detect ERP Components

Detects known event-related potential (ERP) components in epoched (3D)
EEG data. Averages across epochs and finds peaks within predefined time
windows.

## Usage

``` r
eegERPdetect(
  x,
  component = c("N100", "P300", "N400", "P600", "MMN", "LPP"),
  channels = NULL,
  epoch_start = 0,
  assay_name = NULL
)
```

## Arguments

- x:

  A PhysioExperiment object with epoched (3D) EEG data.

- component:

  ERP component to detect: `"N100"`, `"P300"`, `"N400"`, `"P600"`,
  `"MMN"`, or `"LPP"`.

- channels:

  Character vector of channel labels to analyze. If `NULL`, all channels
  are used.

- epoch_start:

  Start time of the epoch in milliseconds relative to stimulus onset.
  Default is 0 (epoch starts at stimulus).

- assay_name:

  Input assay name (default: first assay).

## Value

A data.frame with columns: `channel` (character label), `component`
(character name), `latency_ms` (numeric peak latency in milliseconds),
and `amplitude` (numeric peak amplitude).

## References

Luck, S. J. (2014). An Introduction to the Event-Related Potential
Technique (2nd ed.). MIT Press.

## See also

[`eegERPmeasure()`](https://x-biosignal.github.io/PhysioEEG/reference/eegERPmeasure.md),
[`eegERPlatency()`](https://x-biosignal.github.io/PhysioEEG/reference/eegERPlatency.md),
[`eegERPbaseline()`](https://x-biosignal.github.io/PhysioEEG/reference/eegERPbaseline.md),
[`eegEpoch()`](https://x-biosignal.github.io/PhysioEEG/reference/eegEpoch.md),
[`eegPlotERP()`](https://x-biosignal.github.io/PhysioEEG/reference/eegPlotERP.md)

## Examples

``` r
if (FALSE) { # \dontrun{
pe <- make_eeg_erp(n_epochs = 40, sr = 250)
result <- eegERPdetect(pe, component = "P300")
} # }
```

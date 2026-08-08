# Measure ERP Amplitude

Measures ERP amplitude in a specified time window using peak, mean, or
adaptive mean methods. For epoched (3D) data, averages across epochs
first.

## Usage

``` r
eegERPmeasure(
  x,
  window,
  method = c("peak", "mean", "adaptive_mean"),
  polarity = c("positive", "negative"),
  epoch_start = 0,
  assay_name = NULL
)
```

## Arguments

- x:

  A PhysioExperiment object with EEG data.

- window:

  Numeric vector of length 2: `c(start_ms, end_ms)`.

- method:

  Measurement method: `"peak"` (peak amplitude), `"mean"` (mean
  amplitude), or `"adaptive_mean"` (mean in +/-25ms around peak).

- polarity:

  Expected polarity: `"positive"` or `"negative"`.

- epoch_start:

  Start time of the epoch in milliseconds relative to stimulus onset.
  Default is 0 (epoch starts at stimulus).

- assay_name:

  Input assay name (default: first assay).

## Value

A data.frame with columns: `channel` (character label), `amplitude`
(numeric measured amplitude), `latency_ms` (numeric latency in
milliseconds), and `method` (character measurement method used).

## References

Luck, S. J. (2014). An Introduction to the Event-Related Potential
Technique (2nd ed.). MIT Press.

## See also

[`eegERPdetect()`](https://x-biosignal.github.io/PhysioEEG/reference/eegERPdetect.md),
[`eegERPlatency()`](https://x-biosignal.github.io/PhysioEEG/reference/eegERPlatency.md),
[`eegERPbaseline()`](https://x-biosignal.github.io/PhysioEEG/reference/eegERPbaseline.md),
[`eegPlotERP()`](https://x-biosignal.github.io/PhysioEEG/reference/eegPlotERP.md)

## Examples

``` r
if (FALSE) { # \dontrun{
pe <- make_eeg_erp(n_epochs = 40, sr = 250)
result <- eegERPmeasure(pe, window = c(250, 500), method = "peak",
                        polarity = "positive")
} # }
```

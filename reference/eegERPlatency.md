# Fractional Area Latency

Computes the fractional area latency of an ERP component. This is the
time point at which a specified fraction of the total area under the
curve (in the given window) has accumulated.

## Usage

``` r
eegERPlatency(
  x,
  window,
  fraction = 0.5,
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

- fraction:

  Fraction of total area (default: 0.5 for median latency).

- polarity:

  Expected polarity: `"positive"` or `"negative"`.

- epoch_start:

  Start time of the epoch in milliseconds relative to stimulus onset.
  Default is 0 (epoch starts at stimulus).

- assay_name:

  Input assay name (default: first assay).

## Value

A data.frame with columns: `channel` (character label), `latency_ms`
(numeric fractional area latency in milliseconds), and `fraction`
(numeric fraction used).

## References

Luck, S. J. (2014). An Introduction to the Event-Related Potential
Technique (2nd ed.). MIT Press.

## See also

[`eegERPdetect()`](https://x-biosignal.github.io/PhysioEEG/reference/eegERPdetect.md),
[`eegERPmeasure()`](https://x-biosignal.github.io/PhysioEEG/reference/eegERPmeasure.md),
[`eegERPbaseline()`](https://x-biosignal.github.io/PhysioEEG/reference/eegERPbaseline.md),
[`eegPlotERP()`](https://x-biosignal.github.io/PhysioEEG/reference/eegPlotERP.md)

## Examples

``` r
if (FALSE) { # \dontrun{
pe <- make_eeg_erp(n_epochs = 40, sr = 250)
result <- eegERPlatency(pe, window = c(250, 500), fraction = 0.5,
                        polarity = "positive")
} # }
```

# Automatic Sleep Staging

Classifies EEG epochs into sleep stages (Wake, N1, N2, N3, REM) using
spectral power analysis based on simplified AASM criteria. Each epoch is
scored by computing spectral band powers (delta, theta, alpha, sigma,
beta) and applying rule-based classification.

## Usage

``` r
eegSleepStage(
  x,
  method = c("spectral", "rule_based"),
  epoch_sec = 30,
  assay_name = NULL
)
```

## Arguments

- x:

  A PhysioExperiment object with EEG data.

- method:

  Staging method: `"spectral"` (FFT power-based) or `"rule_based"`
  (synonym, same algorithm).

- epoch_sec:

  Epoch duration in seconds (default: 30, AASM standard).

- assay_name:

  Input assay name. If `NULL`, uses the default assay.

## Value

A data.frame with columns:

- epoch:

  Integer epoch number.

- stage:

  Character sleep stage: "W", "N1", "N2", "N3", or "REM".

- start_sample:

  Integer start sample of the epoch.

- end_sample:

  Integer end sample of the epoch.

- delta_power:

  Numeric delta band power (0.5-4 Hz).

- theta_power:

  Numeric theta band power (4-8 Hz).

- alpha_power:

  Numeric alpha band power (8-13 Hz).

- sigma_power:

  Numeric sigma band power (12-16 Hz).

- beta_power:

  Numeric beta band power (16-30 Hz).

The result is also stored in `metadata(x)$sleep_stages`.

## References

Berry, R. B., et al. (2017). AASM Scoring Manual Updates for 2017.
Journal of Clinical Sleep Medicine, 13(5), 665-666.

## See also

[`eegSpindleDetect()`](https://x-biosignal.github.io/PhysioEEG/reference/eegSpindleDetect.md),
[`eegKcomplexDetect()`](https://x-biosignal.github.io/PhysioEEG/reference/eegKcomplexDetect.md),
[`eegSlowWaveDetect()`](https://x-biosignal.github.io/PhysioEEG/reference/eegSlowWaveDetect.md),
[`eegSleepMetrics()`](https://x-biosignal.github.io/PhysioEEG/reference/eegSleepMetrics.md),
[`eegPlotHypnogram()`](https://x-biosignal.github.io/PhysioEEG/reference/eegPlotHypnogram.md)

## Examples

``` r
if (FALSE) { # \dontrun{
pe <- make_eeg_sleep(n_time = 90000, n_channels = 2, sr = 500)
stages <- eegSleepStage(pe, epoch_sec = 30)
head(stages)
} # }
```

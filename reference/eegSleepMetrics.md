# Compute Sleep Metrics

Calculates summary sleep architecture metrics from staged sleep data.
Requires that
[`eegSleepStage()`](https://x-biosignal.github.io/PhysioEEG/reference/eegSleepStage.md)
has been run and results are stored in `metadata(x)$sleep_stages`.

## Usage

``` r
eegSleepMetrics(x)
```

## Arguments

- x:

  A PhysioExperiment object with sleep staging results in
  `metadata(x)$sleep_stages`.

## Value

A data.frame with columns:

- metric:

  Character name of the sleep metric.

- value:

  Numeric value of the metric.

- unit:

  Character unit of measurement.

Metrics include:

- `total_sleep_time`: Time in N1+N2+N3+REM (minutes)

- `sleep_efficiency`: TST / total recording time \* 100 (percent)

- `waso`: Wake after sleep onset (minutes)

- `sleep_latency`: Time to first non-Wake epoch (minutes)

- `rem_latency`: Time from first sleep to first REM (minutes)

- `pct_N1`: Percentage of TST in N1

- `pct_N2`: Percentage of TST in N2

- `pct_N3`: Percentage of TST in N3

- `pct_REM`: Percentage of TST in REM

- `pct_W`: Percentage of total time in Wake

## References

Berry, R. B., et al. (2017). AASM Scoring Manual Updates for 2017.
Journal of Clinical Sleep Medicine, 13(5), 665-666.

## See also

[`eegSleepStage()`](https://x-biosignal.github.io/PhysioEEG/reference/eegSleepStage.md),
[`eegSpindleDetect()`](https://x-biosignal.github.io/PhysioEEG/reference/eegSpindleDetect.md),
[`eegKcomplexDetect()`](https://x-biosignal.github.io/PhysioEEG/reference/eegKcomplexDetect.md),
[`eegSlowWaveDetect()`](https://x-biosignal.github.io/PhysioEEG/reference/eegSlowWaveDetect.md),
[`eegPlotHypnogram()`](https://x-biosignal.github.io/PhysioEEG/reference/eegPlotHypnogram.md)

## Examples

``` r
if (FALSE) { # \dontrun{
pe <- make_eeg_sleep(n_time = 150000, n_channels = 2, sr = 500)
stages <- eegSleepStage(pe, epoch_sec = 30)
metadata(pe)$sleep_stages <- stages
metrics <- eegSleepMetrics(pe)
print(metrics)
} # }
```

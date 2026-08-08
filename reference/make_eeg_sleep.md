# Create EEG with Sleep Stage Characteristics

Generates continuous EEG with segments having distinct frequency content
matching AASM sleep stages: Wake (alpha), N1 (theta), N2 (spindles + K),
N3 (delta/SWA), REM (mixed low-voltage). Each 30-second epoch is
assigned a stage in a cyclic pattern.

## Usage

``` r
make_eeg_sleep(n_time = 150000, n_channels = 2, sr = 500)
```

## Arguments

- n_time:

  Number of time points (default: 150000 = 5min at 500Hz).

- n_channels:

  Number of channels (default: 2, C3 and C4).

- sr:

  Sampling rate in Hz (default: 500).

## Value

A PhysioExperiment with simulated sleep EEG in the `"raw"` assay. Sets
`metadata(x)$sleep_stages` with a data.frame containing columns:
`epoch`, `stage` (ground truth), `start_sample`, and `end_sample`.

## See also

[`make_eeg()`](https://x-biosignal.github.io/PhysioEEG/reference/make_eeg.md),
[`eegSleepStage()`](https://x-biosignal.github.io/PhysioEEG/reference/eegSleepStage.md),
[`eegSpindleDetect()`](https://x-biosignal.github.io/PhysioEEG/reference/eegSpindleDetect.md),
[`eegKcomplexDetect()`](https://x-biosignal.github.io/PhysioEEG/reference/eegKcomplexDetect.md),
[`eegSleepMetrics()`](https://x-biosignal.github.io/PhysioEEG/reference/eegSleepMetrics.md)

## Examples

``` r
pe <- make_eeg_sleep(n_time = 30000, n_channels = 2, sr = 500)
dim(SummarizedExperiment::assay(pe, "raw"))  # 30000 x 2
#> [1] 30000     2
```

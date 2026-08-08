# Detect Slow Waves

Identifies slow oscillations in EEG data by bandpass filtering in the
slow wave frequency range, finding zero crossings, and measuring
negative half-wave amplitudes and slopes.

## Usage

``` r
eegSlowWaveDetect(
  x,
  min_amplitude = 75,
  freq_range = c(0.5, 2),
  min_duration_ms = 250,
  max_duration_ms = 1000,
  assay_name = NULL
)
```

## Arguments

- x:

  A PhysioExperiment object with EEG data.

- min_amplitude:

  Minimum absolute negative peak amplitude in microvolts (default: 75).

- freq_range:

  Numeric vector of length 2 specifying the slow wave frequency range in
  Hz (default: `c(0.5, 2)`).

- min_duration_ms:

  Minimum half-wave duration in milliseconds (default: 250).

- max_duration_ms:

  Maximum half-wave duration in milliseconds (default: 1000).

- assay_name:

  Input assay name. If `NULL`, uses the default assay.

## Value

A data.frame with columns:

- channel:

  Integer channel index.

- start_sample:

  Integer sample at first zero crossing.

- end_sample:

  Integer sample at second zero crossing.

- negative_peak:

  Numeric negative peak amplitude.

- positive_peak:

  Numeric positive peak amplitude.

- duration_ms:

  Numeric half-wave duration in milliseconds.

- slope:

  Numeric slope from negative to positive peak (microvolts per
  millisecond).

## References

Berry, R. B., et al. (2017). AASM Scoring Manual Updates for 2017.
Journal of Clinical Sleep Medicine, 13(5), 665-666.

## See also

[`eegSleepStage()`](https://x-biosignal.github.io/PhysioEEG/reference/eegSleepStage.md),
[`eegSpindleDetect()`](https://x-biosignal.github.io/PhysioEEG/reference/eegSpindleDetect.md),
[`eegKcomplexDetect()`](https://x-biosignal.github.io/PhysioEEG/reference/eegKcomplexDetect.md),
[`eegSleepMetrics()`](https://x-biosignal.github.io/PhysioEEG/reference/eegSleepMetrics.md)

## Examples

``` r
if (FALSE) { # \dontrun{
pe <- make_eeg_sleep(n_time = 150000, n_channels = 2, sr = 500)
slow_waves <- eegSlowWaveDetect(pe)
head(slow_waves)
} # }
```

# Detect K-Complexes

Identifies K-complexes in EEG data by lowpass filtering at 4 Hz, finding
negative peaks exceeding a threshold amplitude, and verifying the
characteristic negative-positive waveform morphology.

## Usage

``` r
eegKcomplexDetect(
  x,
  min_neg_amplitude = 75,
  min_duration_ms = 500,
  max_duration_ms = 1500,
  assay_name = NULL
)
```

## Arguments

- x:

  A PhysioExperiment object with EEG data.

- min_neg_amplitude:

  Minimum absolute negative peak amplitude in microvolts (default: 75).
  Peaks must be more negative than `-min_neg_amplitude`.

- min_duration_ms:

  Minimum K-complex duration in milliseconds (default: 500).

- max_duration_ms:

  Maximum K-complex duration in milliseconds (default: 1500).

- assay_name:

  Input assay name. If `NULL`, uses the default assay.

## Value

A data.frame with columns:

- channel:

  Integer channel index.

- negative_peak_sample:

  Integer sample of the negative peak.

- positive_peak_sample:

  Integer sample of the positive peak.

- negative_amplitude:

  Numeric amplitude at the negative peak.

- positive_amplitude:

  Numeric amplitude at the positive peak.

- duration_ms:

  Numeric total duration in milliseconds.

## References

Berry, R. B., et al. (2017). AASM Scoring Manual Updates for 2017.
Journal of Clinical Sleep Medicine, 13(5), 665-666.

## See also

[`eegSleepStage()`](https://x-biosignal.github.io/PhysioEEG/reference/eegSleepStage.md),
[`eegSpindleDetect()`](https://x-biosignal.github.io/PhysioEEG/reference/eegSpindleDetect.md),
[`eegSlowWaveDetect()`](https://x-biosignal.github.io/PhysioEEG/reference/eegSlowWaveDetect.md),
[`eegSleepMetrics()`](https://x-biosignal.github.io/PhysioEEG/reference/eegSleepMetrics.md)

## Examples

``` r
if (FALSE) { # \dontrun{
pe <- make_eeg_sleep(n_time = 150000, n_channels = 2, sr = 500)
kcomplexes <- eegKcomplexDetect(pe)
head(kcomplexes)
} # }
```

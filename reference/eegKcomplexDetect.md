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
pe <- make_eeg_sleep(n_time = 150000, n_channels = 2, sr = 500)
kcomplexes <- eegKcomplexDetect(pe)
head(kcomplexes)
#>   channel negative_peak_sample positive_peak_sample negative_amplitude
#> 1       1               135267               135990          -79.77818
#> 2       1               136267               136990          -79.44237
#> 3       1               137266               137991          -79.58419
#> 4       1               138267               138992          -79.96084
#> 5       1               139267               139990          -79.41695
#> 6       1               140266               140992          -79.85121
#>   positive_amplitude duration_ms
#> 1           56.01905        1446
#> 2           56.40773        1446
#> 3           56.18535        1450
#> 4           56.37444        1450
#> 5           55.69484        1446
#> 6           56.07785        1452
```

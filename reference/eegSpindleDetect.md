# Detect Sleep Spindles

Identifies sleep spindles in EEG data by bandpass filtering in the sigma
frequency range, computing an RMS envelope, and detecting continuous
segments exceeding a threshold. Spindle frequency is estimated from zero
crossings in the bandpassed signal.

## Usage

``` r
eegSpindleDetect(
  x,
  method = c("sigma", "wavelet"),
  freq_range = c(11, 16),
  min_duration_ms = 500,
  max_duration_ms = 2000,
  threshold_sd = 1.5,
  assay_name = NULL
)
```

## Arguments

- x:

  A PhysioExperiment object with EEG data.

- method:

  Detection method: `"sigma"` (sigma-band filtering) or `"wavelet"`
  (synonym, same algorithm).

- freq_range:

  Numeric vector of length 2 specifying the spindle frequency range in
  Hz (default: `c(11, 16)`).

- min_duration_ms:

  Minimum spindle duration in milliseconds (default: 500).

- max_duration_ms:

  Maximum spindle duration in milliseconds (default: 2000).

- threshold_sd:

  Number of standard deviations above the mean for detection threshold
  (default: 1.5).

- assay_name:

  Input assay name. If `NULL`, uses the default assay.

## Value

A data.frame with columns:

- channel:

  Integer channel index.

- start_sample:

  Integer start sample of the spindle.

- end_sample:

  Integer end sample of the spindle.

- duration_ms:

  Numeric spindle duration in milliseconds.

- peak_sample:

  Integer sample index of peak amplitude.

- peak_amplitude:

  Numeric peak RMS amplitude.

- frequency_hz:

  Numeric estimated spindle frequency in Hz.

## References

Berry, R. B., et al. (2017). AASM Scoring Manual Updates for 2017.
Journal of Clinical Sleep Medicine, 13(5), 665-666.

Molle, M., et al. (2011). Fast and slow spindles during the sleep slow
oscillation. Sleep, 34(10), 1411-1421.

## See also

[`eegSleepStage()`](https://x-biosignal.github.io/PhysioEEG/reference/eegSleepStage.md),
[`eegKcomplexDetect()`](https://x-biosignal.github.io/PhysioEEG/reference/eegKcomplexDetect.md),
[`eegSlowWaveDetect()`](https://x-biosignal.github.io/PhysioEEG/reference/eegSlowWaveDetect.md),
[`eegSleepMetrics()`](https://x-biosignal.github.io/PhysioEEG/reference/eegSleepMetrics.md)

## Examples

``` r
if (FALSE) { # \dontrun{
pe <- make_eeg_sleep(n_time = 150000, n_channels = 2, sr = 500)
spindles <- eegSpindleDetect(pe)
head(spindles)
} # }
```

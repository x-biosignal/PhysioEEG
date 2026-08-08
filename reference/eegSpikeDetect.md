# Detect Epileptic Spikes in EEG Data

Identifies epileptiform spike discharges in multi-channel EEG using
either morphology-based derivative analysis or template matching via
cross-correlation. The morphology method detects sharp transients by
thresholding the first derivative, while the template method
cross-correlates a canonical spike waveform with each channel.

## Usage

``` r
eegSpikeDetect(
  x,
  method = c("morphology", "template"),
  threshold_sd = 4,
  min_duration_ms = 20,
  max_duration_ms = 200,
  min_amplitude = 50,
  assay_name = NULL
)
```

## Arguments

- x:

  A PhysioExperiment object with EEG data.

- method:

  Detection method: `"morphology"` (derivative-based) or `"template"`
  (cross-correlation with canonical spike shape).

- threshold_sd:

  Number of standard deviations above the mean derivative magnitude for
  detection (default: 4).

- min_duration_ms:

  Minimum spike duration in milliseconds (default: 20).

- max_duration_ms:

  Maximum spike duration in milliseconds (default: 200).

- min_amplitude:

  Minimum peak amplitude in microvolts for a valid spike detection
  (default: 50).

- assay_name:

  Name of the assay to use. If `NULL`, the default assay is used.

## Value

A data.frame with columns:

- channel:

  Integer channel index.

- sample:

  Integer sample index of the spike peak.

- time_sec:

  Time of the spike in seconds.

- amplitude:

  Peak amplitude at the spike location.

- duration_ms:

  Estimated spike duration in milliseconds.

- confidence:

  Confidence score for the detection.

## References

Nuwer, M. R., et al. (1999). IFCN standards for digital recording of
clinical EEG. Electroencephalography and Clinical Neurophysiology,
106(3), 259-261.

## See also

[`eegQEEG()`](https://x-biosignal.github.io/PhysioEEG/reference/eegQEEG.md),
[`eegAsymmetry()`](https://x-biosignal.github.io/PhysioEEG/reference/eegAsymmetry.md),
[`eegSlowing()`](https://x-biosignal.github.io/PhysioEEG/reference/eegSlowing.md)

## Examples

``` r
if (FALSE) { # \dontrun{
pe <- make_eeg_spikes(n_time = 30000, n_channels = 19, sr = 500, n_spikes = 15)
spikes <- eegSpikeDetect(pe, method = "morphology")
head(spikes)
} # }
```

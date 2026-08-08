# Burst-Suppression Detection

Identifies periods of burst and suppression in EEG data by computing the
root mean square (RMS) amplitude in sliding windows. Suppression
segments are defined as consecutive windows where RMS falls below the
specified threshold for at least `min_duration_ms`.

## Usage

``` r
eegSuppression(
  x,
  threshold = 10,
  min_duration_ms = 500,
  window_ms = 500,
  assay_name = NULL
)
```

## Arguments

- x:

  A PhysioExperiment object with EEG data.

- threshold:

  RMS amplitude threshold below which signal is considered suppressed
  (default: 10).

- min_duration_ms:

  Minimum duration in milliseconds for a segment to be classified as
  suppression (default: 500).

- window_ms:

  Window length in milliseconds for RMS computation (default: 500).

- assay_name:

  Name of the input assay. If `NULL`, the default assay is used.

## Value

A data.frame with columns:

- type:

  Character: `"burst"` or `"suppression"`.

- start_sample:

  Integer start sample of the segment.

- end_sample:

  Integer end sample of the segment.

- duration_ms:

  Numeric duration of the segment in milliseconds.

An attribute `"bsr"` (Burst-Suppression Ratio) is attached, representing
the percentage of total time in suppression.

## References

Nuwer, M. R., et al. (1999). IFCN standards for digital recording of
clinical EEG. Electroencephalography and Clinical Neurophysiology,
106(3), 259-261.

## See also

[`eegSpikeDetect()`](https://x-biosignal.github.io/PhysioEEG/reference/eegSpikeDetect.md),
[`eegQEEG()`](https://x-biosignal.github.io/PhysioEEG/reference/eegQEEG.md),
[`eegSlowing()`](https://x-biosignal.github.io/PhysioEEG/reference/eegSlowing.md)

## Examples

``` r
if (FALSE) { # \dontrun{
pe <- make_eeg(n_time = 10000, n_channels = 4, sr = 500)
result <- eegSuppression(pe, threshold = 10)
print(attr(result, "bsr"))
} # }
```

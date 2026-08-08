# EEG Slowing Detection

Detects pathological EEG slowing using spectral analysis. Supports three
methods: theta/delta ratio (TDR), delta-theta/alpha-beta ratio (DTAR),
and peak frequency analysis. Each channel is classified into normal,
mild, moderate, or severe slowing categories.

## Usage

``` r
eegSlowing(
  x,
  method = c("tdr", "dtar", "peak_frequency"),
  bands = NULL,
  assay_name = NULL
)
```

## Arguments

- x:

  A PhysioExperiment object with EEG data.

- method:

  Analysis method: `"tdr"` (theta/delta ratio), `"dtar"` (delta+theta
  over alpha+beta ratio), or `"peak_frequency"` (dominant frequency per
  channel).

- bands:

  Named list of frequency bands for TDR and DTAR methods. Defaults to
  standard EEG bands.

- assay_name:

  Name of the input assay. If `NULL`, the default assay is used.

## Value

A data.frame with columns:

- channel:

  Integer channel index.

- metric:

  Character name of the metric used.

- value:

  Numeric value of the metric.

- classification:

  Character: `"normal"`, `"mild_slowing"`, `"moderate_slowing"`, or
  `"severe_slowing"`.

## References

Nuwer, M. R., et al. (1999). IFCN standards for digital recording of
clinical EEG. Electroencephalography and Clinical Neurophysiology,
106(3), 259-261.

## See also

[`eegSpikeDetect()`](https://x-biosignal.github.io/PhysioEEG/reference/eegSpikeDetect.md),
[`eegQEEG()`](https://x-biosignal.github.io/PhysioEEG/reference/eegQEEG.md),
[`eegAsymmetry()`](https://x-biosignal.github.io/PhysioEEG/reference/eegAsymmetry.md),
[`eegSuppression()`](https://x-biosignal.github.io/PhysioEEG/reference/eegSuppression.md)

## Examples

``` r
if (FALSE) { # \dontrun{
pe <- make_eeg(n_time = 5000, n_channels = 19, sr = 500)
result <- eegSlowing(pe, method = "dtar")
print(result)
} # }
```

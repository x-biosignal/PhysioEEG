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
pe <- make_eeg(n_time = 5000, n_channels = 19, sr = 500)
result <- eegSlowing(pe, method = "dtar")
print(result)
#>    channel metric    value   classification
#> 1        1   dtar 1.281362     mild_slowing
#> 2        2   dtar 1.009824     mild_slowing
#> 3        3   dtar 1.226197     mild_slowing
#> 4        4   dtar 1.216584     mild_slowing
#> 5        5   dtar 1.045893     mild_slowing
#> 6        6   dtar 1.399565     mild_slowing
#> 7        7   dtar 1.359063     mild_slowing
#> 8        8   dtar 1.194493     mild_slowing
#> 9        9   dtar 2.576911 moderate_slowing
#> 10      10   dtar 2.041697 moderate_slowing
#> 11      11   dtar 3.418652 moderate_slowing
#> 12      12   dtar 5.516801   severe_slowing
#> 13      13   dtar 5.362880   severe_slowing
#> 14      14   dtar 5.119384   severe_slowing
#> 15      15   dtar 4.104294   severe_slowing
#> 16      16   dtar 4.415537   severe_slowing
#> 17      17   dtar 3.791294 moderate_slowing
#> 18      18   dtar 2.083023 moderate_slowing
#> 19      19   dtar 1.516108     mild_slowing
```

# Quantitative EEG (QEEG) Analysis

Computes absolute and relative spectral band powers for each channel
using Welch's method (windowed FFT averaging). Results are stored as a
new assay (channels x bands matrix of absolute power) and as relative
power in `metadata(x)$qeeg`.

## Usage

``` r
eegQEEG(
  x,
  bands = NULL,
  window_sec = 2,
  overlap = 0.5,
  assay_name = NULL,
  output_assay = "qeeg"
)
```

## Arguments

- x:

  A PhysioExperiment object with EEG data.

- bands:

  Named list of frequency bands. Each element is a numeric vector of
  length 2 specifying the lower and upper frequency in Hz. Defaults to
  standard EEG bands: delta (1-4), theta (4-8), alpha (8-13), beta
  (13-30), gamma (30-50).

- window_sec:

  Window length in seconds for Welch's method (default: 2).

- overlap:

  Overlap fraction between windows, 0 to 1 (default: 0.5).

- assay_name:

  Name of the input assay. If `NULL`, the default assay is used.

- output_assay:

  Name of the assay to store absolute power results (default: `"qeeg"`).

## Value

Modified PhysioExperiment with:

- Absolute power matrix (n_channels x n_bands) in `output_assay`

- Band definitions and relative power in `metadata(x)$qeeg`, a list
  containing `bands`, `absolute_power`, `relative_power`, `band_names`,
  `window_sec`, and `overlap`.

## References

Nuwer, M. R., et al. (1999). IFCN standards for digital recording of
clinical EEG. Electroencephalography and Clinical Neurophysiology,
106(3), 259-261.

Thatcher, R. W. (2010). Validity and reliability of quantitative
electroencephalography. Journal of Neurotherapy, 14(2), 122-152.

## See also

[`eegSpikeDetect()`](https://x-biosignal.github.io/PhysioEEG/reference/eegSpikeDetect.md),
[`eegAsymmetry()`](https://x-biosignal.github.io/PhysioEEG/reference/eegAsymmetry.md),
[`eegSlowing()`](https://x-biosignal.github.io/PhysioEEG/reference/eegSlowing.md),
[`eegPlotSpectrogram()`](https://x-biosignal.github.io/PhysioEEG/reference/eegPlotSpectrogram.md)

## Examples

``` r
pe <- make_eeg(n_time = 5000, n_channels = 19, sr = 500)
pe_qeeg <- eegQEEG(pe)
qeeg_info <- S4Vectors::metadata(pe_qeeg)$qeeg
print(qeeg_info$relative_power)
#>           delta      theta      alpha       beta       gamma
#>  [1,] 0.5954706 0.04304591 0.29414593 0.06227962 0.005057907
#>  [2,] 0.5634851 0.07171742 0.30719517 0.05157654 0.006025819
#>  [3,] 0.4744163 0.05540047 0.41891176 0.04610726 0.005164219
#>  [4,] 0.3812155 0.06129063 0.50208595 0.04923432 0.006173583
#>  [5,] 0.5089776 0.05274212 0.40125524 0.03211365 0.004911432
#>  [6,] 0.5611147 0.05214036 0.36220007 0.01972321 0.004821617
#>  [7,] 0.5156311 0.04201829 0.41134118 0.02504950 0.005959922
#>  [8,] 0.4725432 0.07901882 0.41077332 0.02905529 0.008609424
#>  [9,] 0.6196409 0.06051027 0.29009771 0.02257753 0.007173592
#> [10,] 0.6511383 0.05000766 0.26404149 0.02722915 0.007583403
#> [11,] 0.6283563 0.09565181 0.23094485 0.03471872 0.010328280
#> [12,] 0.6981922 0.08997026 0.16072007 0.04090029 0.010217171
#> [13,] 0.7981662 0.07503894 0.08639502 0.03330469 0.007095153
#> [14,] 0.6547896 0.10643764 0.15598788 0.07088903 0.011895825
#> [15,] 0.7011096 0.09784320 0.11610948 0.07523716 0.009700563
#> [16,] 0.7317927 0.08424687 0.12074257 0.05580485 0.007413004
#> [17,] 0.6913388 0.06266133 0.16300996 0.07695943 0.006030498
#> [18,] 0.6110704 0.06133510 0.21951613 0.09910003 0.008978304
#> [19,] 0.6688057 0.06392462 0.19725216 0.06445321 0.005564346
```

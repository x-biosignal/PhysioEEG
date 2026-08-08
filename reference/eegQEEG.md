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
if (FALSE) { # \dontrun{
pe <- make_eeg(n_time = 5000, n_channels = 19, sr = 500)
pe_qeeg <- eegQEEG(pe)
qeeg_info <- metadata(pe_qeeg)$qeeg
print(qeeg_info$relative_power)
} # }
```

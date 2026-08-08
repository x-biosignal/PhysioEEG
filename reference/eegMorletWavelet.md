# Morlet Wavelet Transform for EEG

Computes the continuous Morlet wavelet transform for multi-channel EEG
data. For each specified frequency, a complex Morlet wavelet is
constructed and convolved with each channel using FFT-based convolution
for efficiency. Returns time-resolved power (and optionally phase)
across frequencies.

## Usage

``` r
eegMorletWavelet(
  x,
  frequencies = NULL,
  n_cycles = 7,
  assay_name = NULL,
  output_assay = "wavelet_power"
)
```

## Arguments

- x:

  A PhysioExperiment object with EEG data (2D: time x channels).

- frequencies:

  Numeric vector of frequencies in Hz to analyze. If `NULL`, defaults to
  `seq(1, 50, by = 1)`.

- n_cycles:

  Number of cycles in the Morlet wavelet, controlling the trade-off
  between time and frequency resolution (default: 7).

- assay_name:

  Name of the input assay. If `NULL`, the default assay is used.

- output_assay:

  Name of the assay to store wavelet power results (default:
  `"wavelet_power"`).

## Value

Modified PhysioExperiment with:

- 3D power array (time x frequencies x channels) in `output_assay`

- Frequency vector and phase array in `metadata(x)$wavelet`, a list
  containing `frequencies` (numeric vector), `n_cycles` (integer), and
  `phase` (3D array of instantaneous phase values)

## References

Tallon-Baudry, C., et al. (1997). Oscillatory gamma-band activity during
conscious perception. Trends in Cognitive Sciences, 3(4), 151-162.

## See also

[`eegSTFT()`](https://x-biosignal.github.io/PhysioEEG/reference/eegSTFT.md),
[`eegMultitaper()`](https://x-biosignal.github.io/PhysioEEG/reference/eegMultitaper.md),
[`eegPlotSpectrogram()`](https://x-biosignal.github.io/PhysioEEG/reference/eegPlotSpectrogram.md),
[`eegERSP()`](https://x-biosignal.github.io/PhysioEEG/reference/eegERSP.md)

## Examples

``` r
if (FALSE) { # \dontrun{
pe <- make_eeg(n_time = 5000, n_channels = 4, sr = 500)
pe_wt <- eegMorletWavelet(pe, frequencies = seq(5, 40, by = 1))
wp <- SummarizedExperiment::assay(pe_wt, "wavelet_power")
dim(wp)  # time x frequencies x channels
} # }
```

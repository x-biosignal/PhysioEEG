# Short-Time Fourier Transform for EEG

Computes the Short-Time Fourier Transform (STFT) spectrogram for
multi-channel EEG data. Uses a sliding window with configurable overlap
and window function to produce a time-frequency power representation.

## Usage

``` r
eegSTFT(
  x,
  window_sec = 0.5,
  overlap = 0.75,
  window_type = c("hanning", "hamming", "rectangular"),
  assay_name = NULL,
  output_assay = "stft_power"
)
```

## Arguments

- x:

  A PhysioExperiment object with EEG data (2D: time x channels).

- window_sec:

  Window length in seconds (default: 0.5).

- overlap:

  Overlap fraction between adjacent windows, from 0 to 1 exclusive
  (default: 0.75).

- window_type:

  Window function to apply: `"hanning"`, `"hamming"`, or `"rectangular"`
  (default: `"hanning"`).

- assay_name:

  Name of the input assay. If `NULL`, the default assay is used.

- output_assay:

  Name of the assay to store STFT power results (default:
  `"stft_power"`).

## Value

Modified PhysioExperiment with:

- 3D power array (time_bins x frequencies x channels) in `output_assay`

- Time bin centers, frequency vector, and parameters in
  `metadata(x)$stft`, a list containing `time_axis`, `freq_axis`,
  `window_sec`, `overlap`, `window_type`, `window_length`, and
  `hop_size`

## References

Tallon-Baudry, C., et al. (1997). Oscillatory gamma-band activity during
conscious perception. Trends in Cognitive Sciences, 3(4), 151-162.

## See also

[`eegMorletWavelet()`](https://x-biosignal.github.io/PhysioEEG/reference/eegMorletWavelet.md),
[`eegMultitaper()`](https://x-biosignal.github.io/PhysioEEG/reference/eegMultitaper.md),
[`eegPlotSpectrogram()`](https://x-biosignal.github.io/PhysioEEG/reference/eegPlotSpectrogram.md)

## Examples

``` r
if (FALSE) { # \dontrun{
pe <- make_eeg(n_time = 5000, n_channels = 4, sr = 500)
pe_stft <- eegSTFT(pe, window_sec = 0.5, overlap = 0.75)
sp <- SummarizedExperiment::assay(pe_stft, "stft_power")
dim(sp)  # time_bins x frequencies x channels
} # }
```

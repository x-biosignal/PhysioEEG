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

  Name the result is stored under. Despite the argument's name this is a
  key in
  [`metadata`](https://rdrr.io/pkg/S4Vectors/man/Annotated-class.html),
  NOT an assay: after the call `assay(x, "the value given")` does not
  exist while `metadata(x)$the value given` holds the result (default:
  `"the value given"`). (default: `"stft_power"`).

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
pe <- make_eeg(n_time = 2000, n_channels = 4, sr = 250)
pe_stft <- eegSTFT(pe, window_sec = 0.5, overlap = 0.75)
# results are stored in metadata(x)$stft
str(S4Vectors::metadata(pe_stft)$stft, max.level = 1)
#> List of 7
#>  $ time_axis    : num [1:59] 0.252 0.38 0.508 0.636 0.764 ...
#>  $ freq_axis    : num [1:64] 0 1.98 3.97 5.95 7.94 ...
#>  $ window_sec   : num 0.5
#>  $ overlap      : num 0.75
#>  $ window_type  : chr "hanning"
#>  $ window_length: int 126
#>  $ hop_size     : int 32
```

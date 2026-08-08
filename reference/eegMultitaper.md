# Multitaper Power Spectral Density for EEG

Computes the multitaper power spectral density (PSD) estimate for
multi-channel EEG data using Discrete Prolate Spheroidal Sequences
(DPSS, Slepian tapers). The multitaper method provides an optimal
bias-variance trade-off for spectral estimation compared to
single-window methods.

## Usage

``` r
eegMultitaper(
  x,
  bandwidth = 4,
  n_tapers = NULL,
  assay_name = NULL,
  output_assay = "multitaper_psd"
)
```

## Arguments

- x:

  A PhysioExperiment object with EEG data (2D: time x channels).

- bandwidth:

  Time-half-bandwidth product (NW) controlling the spectral
  concentration of the tapers (default: 4). Higher values give smoother
  but lower-resolution estimates.

- n_tapers:

  Number of DPSS tapers to use. If `NULL`, defaults to
  `floor(2 * bandwidth) - 1`.

- assay_name:

  Name of the input assay. If `NULL`, the default assay is used.

- output_assay:

  Name of the assay to store PSD results (default: `"multitaper_psd"`).

## Value

Modified PhysioExperiment with:

- PSD matrix (frequencies x channels) in `output_assay`

- Frequency vector and taper parameters in `metadata(x)$multitaper`, a
  list containing `frequencies` (numeric vector), `bandwidth` (numeric
  NW parameter), and `n_tapers` (integer)

## References

Tallon-Baudry, C., et al. (1997). Oscillatory gamma-band activity during
conscious perception. Trends in Cognitive Sciences, 3(4), 151-162.

Thomson, D. J. (1982). Spectrum estimation and harmonic analysis.
Proceedings of the IEEE, 70(9), 1055-1096.

## See also

[`eegMorletWavelet()`](https://x-biosignal.github.io/PhysioEEG/reference/eegMorletWavelet.md),
[`eegSTFT()`](https://x-biosignal.github.io/PhysioEEG/reference/eegSTFT.md),
[`eegPlotSpectrogram()`](https://x-biosignal.github.io/PhysioEEG/reference/eegPlotSpectrogram.md)

## Examples

``` r
if (FALSE) { # \dontrun{
pe <- make_eeg(n_time = 5000, n_channels = 19, sr = 500)
pe_mt <- eegMultitaper(pe, bandwidth = 4)
psd <- SummarizedExperiment::assay(pe_mt, "multitaper_psd")
dim(psd)  # frequencies x channels
} # }
```

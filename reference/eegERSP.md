# Event-Related Spectral Perturbation (ERSP)

Computes Event-Related Spectral Perturbation for epoched (3D) EEG data.
ERSP quantifies event-related changes in spectral power relative to a
baseline period, expressed in decibels (dB). Uses the Morlet wavelet
transform to compute time-frequency decomposition for each epoch, then
averages power across epochs and normalizes to baseline.

## Usage

``` r
eegERSP(
  x,
  baseline = c(1, 50),
  frequencies = NULL,
  n_cycles = 7,
  assay_name = NULL,
  output_assay = "ersp_data"
)
```

## Arguments

- x:

  A PhysioExperiment object with epoched (3D) EEG data (time x channels
  x epochs).

- baseline:

  Numeric vector of length 2 specifying the baseline time window in
  sample indices (e.g., `c(1, 50)` for the first 50 samples). Default is
  `c(1, 50)`.

- frequencies:

  Numeric vector of frequencies in Hz to analyze. If `NULL`, defaults to
  `seq(1, 50, by = 1)`.

- n_cycles:

  Number of cycles for the Morlet wavelet (default: 7).

- assay_name:

  Name of the input assay. If `NULL`, the default assay is used.

- output_assay:

  Name of the assay to store ERSP results (default: `"ersp"`).

## Value

Modified PhysioExperiment with:

- 3D ERSP array (time x frequencies x channels) in dB in `output_assay`

- Baseline info and frequency vector in `metadata(x)$ersp`, a list
  containing `frequencies` (numeric vector), `baseline` (numeric
  vector), `n_cycles` (integer), and `n_epochs` (integer)

## References

Tallon-Baudry, C., et al. (1997). Oscillatory gamma-band activity during
conscious perception. Trends in Cognitive Sciences, 3(4), 151-162.

Makeig, S. (1993). Auditory event-related dynamics of the EEG spectrum
and effects of exposure to tones. Electroencephalography and Clinical
Neurophysiology, 86(4), 283-293.

## See also

[`eegITC()`](https://x-biosignal.github.io/PhysioEEG/reference/eegITC.md),
[`eegMorletWavelet()`](https://x-biosignal.github.io/PhysioEEG/reference/eegMorletWavelet.md),
[`eegPlotSpectrogram()`](https://x-biosignal.github.io/PhysioEEG/reference/eegPlotSpectrogram.md)

## Examples

``` r
if (FALSE) { # \dontrun{
pe <- make_eeg_erp(n_epochs = 20, n_channels = 2, sr = 250)
pe_ersp <- eegERSP(pe, baseline = c(1, 50), frequencies = seq(5, 30, by = 5))
ersp_data <- SummarizedExperiment::assay(pe_ersp, "ersp")
dim(ersp_data)  # time x frequencies x channels
} # }
```

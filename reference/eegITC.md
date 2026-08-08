# Inter-Trial Coherence (ITC)

Computes Inter-Trial Coherence (also known as phase-locking factor or
phase-locking value) for epoched (3D) EEG data. ITC measures the
consistency of oscillatory phase across trials at each time-frequency
point. Values range from 0 (completely random phase) to 1 (perfectly
phase-locked across all trials).

## Usage

``` r
eegITC(
  x,
  frequencies = NULL,
  n_cycles = 7,
  assay_name = NULL,
  output_assay = "itc_data"
)
```

## Arguments

- x:

  A PhysioExperiment object with epoched (3D) EEG data (time x channels
  x epochs).

- frequencies:

  Numeric vector of frequencies in Hz to analyze. If `NULL`, defaults to
  `seq(1, 50, by = 1)`.

- n_cycles:

  Number of cycles for the Morlet wavelet (default: 7).

- assay_name:

  Name of the input assay. If `NULL`, the default assay is used.

- output_assay:

  Name of the assay to store ITC results (default: `"itc"`).

## Value

Modified PhysioExperiment with:

- 3D ITC array (time x frequencies x channels) in `output_assay`, values
  in `[0, 1]`

- Frequency vector and parameters in `metadata(x)$itc`, a list
  containing `frequencies` (numeric vector), `n_cycles` (integer), and
  `n_epochs` (integer)

## References

Tallon-Baudry, C., et al. (1997). Oscillatory gamma-band activity during
conscious perception. Trends in Cognitive Sciences, 3(4), 151-162.

Makeig, S. (1993). Auditory event-related dynamics of the EEG spectrum
and effects of exposure to tones. Electroencephalography and Clinical
Neurophysiology, 86(4), 283-293.

## See also

[`eegERSP()`](https://x-biosignal.github.io/PhysioEEG/reference/eegERSP.md),
[`eegMorletWavelet()`](https://x-biosignal.github.io/PhysioEEG/reference/eegMorletWavelet.md),
[`eegPlotSpectrogram()`](https://x-biosignal.github.io/PhysioEEG/reference/eegPlotSpectrogram.md)

## Examples

``` r
if (FALSE) { # \dontrun{
pe <- make_eeg_erp(n_epochs = 40, n_channels = 2, sr = 250)
pe_itc <- eegITC(pe, frequencies = seq(5, 30, by = 5))
itc_data <- SummarizedExperiment::assay(pe_itc, "itc")
dim(itc_data)  # time x frequencies x channels
} # }
```

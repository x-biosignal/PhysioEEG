# Filter EEG signals

Applies frequency-domain filtering to EEG data stored in a
PhysioExperiment object. Supports bandpass, highpass, lowpass, and notch
filtering using either FIR (windowed-sinc) or IIR (Butterworth) methods.

## Usage

``` r
eegFilter(
  x,
  lowcut = NULL,
  highcut = NULL,
  notch = NULL,
  method = c("fir", "iir"),
  order = NULL,
  assay_name = NULL,
  output_assay = "filtered"
)
```

## Arguments

- x:

  A PhysioExperiment object.

- lowcut:

  Low cutoff frequency in Hz. If NULL, no highpass filtering is applied
  (set both lowcut and highcut for bandpass).

- highcut:

  High cutoff frequency in Hz. If NULL, no lowpass filtering is applied
  (set both lowcut and highcut for bandpass).

- notch:

  Notch filter center frequency in Hz (e.g., 50 or 60 for powerline
  noise). If NULL, no notch filter is applied. Bandwidth is +/- 2 Hz
  around center.

- method:

  Filtering method: `"fir"` (default) for windowed-sinc FIR filter, or
  `"iir"` for Butterworth IIR filter.

- order:

  Filter order. For FIR, this is the number of taps (auto-selected if
  NULL). For IIR, this is the Butterworth order (default: 4).

- assay_name:

  Name of the assay to filter. If NULL, uses `defaultAssay(x)`.

- output_assay:

  Name of the output assay (default: `"filtered"`).

## Value

A PhysioExperiment object with filtered data in the specified output
assay.

## Details

For FIR mode, the function uses windowed-sinc filters with a Hamming
window. For IIR mode, zero-phase Butterworth filtering is applied via
[`signal::filtfilt()`](https://rdrr.io/pkg/signal/man/filtfilt.html),
with automatic fallback to FIR if the signal package is not available.

## Examples

``` r
if (FALSE) { # \dontrun{
pe <- make_eeg(n_time = 5000, n_channels = 19, sr = 500)
# Bandpass filter
pe_filt <- eegFilter(pe, lowcut = 1, highcut = 40)
# Highpass only
pe_hp <- eegFilter(pe, lowcut = 0.1)
# With notch at 50 Hz
pe_notch <- eegFilter(pe, lowcut = 1, highcut = 40, notch = 50)
# IIR Butterworth
pe_iir <- eegFilter(pe, lowcut = 1, highcut = 40, method = "iir", order = 4)
} # }
```

# Windowed-sinc FIR notch (band-stop) filter

Implements a band-stop filter as allpass minus bandpass. The bandpass
kernel at the specified center frequency +/- bandwidth/2 is subtracted
from a unit impulse to produce a notch response.

## Usage

``` r
.fir_notch(signal, sr, center_freq, bandwidth = 4, order = NULL)
```

## Arguments

- signal:

  Numeric vector of the input signal.

- sr:

  Sampling rate in Hz.

- center_freq:

  Center frequency to reject in Hz.

- bandwidth:

  Full bandwidth of the notch in Hz (default: 4).

- order:

  Filter order. If `NULL`, auto-selected.

## Value

Numeric vector of the filtered signal (same length as input).

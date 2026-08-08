# Windowed-sinc FIR lowpass filter

Applies a windowed-sinc FIR lowpass filter using a Hamming window. The
filter is applied causally via
[`stats::filter`](https://rdrr.io/r/stats/filter.html) with group-delay
compensation. This replaces the previous brick-wall FFT filter which
produced Gibbs ringing artifacts.

## Usage

``` r
.fir_lowpass(signal, sr, cutoff, order = NULL)
```

## Arguments

- signal:

  Numeric vector of the input signal.

- sr:

  Sampling rate in Hz.

- cutoff:

  Cutoff frequency in Hz.

- order:

  Filter order. If `NULL`, auto-selected as approximately 3 cycles of
  the cutoff frequency, clamped to signal length and forced odd.

## Value

Numeric vector of the filtered signal (same length as input).

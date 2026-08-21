# Comodulogram (phase-frequency x amplitude-frequency PAC) for EEG

Computes phase-amplitude coupling over a grid of phase and amplitude
frequencies for one EEG channel, yielding a comodulogram whose peak
locates the dominant coupling. Delegates to
[`PhysioCrossModal::comodulogram()`](https://x-biosignal.github.io/PhysioCrossModal/reference/comodulogram.html).

## Usage

``` r
eegComodulogram(
  pe,
  channel = 1L,
  phase_freqs = seq(2, 14, by = 2),
  amp_freqs = seq(20, 100, by = 10),
  method = c("tort", "canolty", "ozkurt", "plv"),
  phase_bw = 2,
  amp_bw = 10,
  n_bins = 18L,
  assay_name = NULL
)
```

## Arguments

- pe:

  A `PhysioExperiment`.

- channel:

  Channel (index or label) to analyse.

- phase_freqs:

  Phase centre frequencies in Hz (default `seq(2, 14, 2)`).

- amp_freqs:

  Amplitude centre frequencies in Hz (default `seq(20, 100, 10)`).

- method:

  PAC estimator (see
  [`eegPAC()`](https://x-biosignal.github.io/PhysioEEG/reference/eegPAC.md)).

- phase_bw, amp_bw:

  Half-bandwidths in Hz (defaults 2 and 10).

- n_bins:

  Phase bins for the Tort modulation index (default 18).

- assay_name:

  Assay to use (default: the object's default assay).

## Value

The comodulogram list from
[`PhysioCrossModal::comodulogram()`](https://x-biosignal.github.io/PhysioCrossModal/reference/comodulogram.html)
(a `matrix` of phase x amplitude frequencies, plus its `peak`).

## See also

[`eegPAC()`](https://x-biosignal.github.io/PhysioEEG/reference/eegPAC.md)

## Examples

``` r
if (FALSE) { # \dontrun{
pe <- make_eeg(n_time = 7500, n_channels = 2, sr = 250)
cm <- eegComodulogram(pe, channel = 1)
cm$peak
} # }
```

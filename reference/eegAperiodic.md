# Aperiodic (1/f) spectral parameterization of EEG

Separates each channel's power spectrum into an **aperiodic** (1/f)
component and **periodic** (oscillatory) peaks, following the specparam
/ FOOOF model (Donoghue et al. 2020). The aperiodic exponent is a widely
used index of the excitation/inhibition balance and cortical state.
Delegates the fit to
[`PhysioAnalysis::specparam()`](https://x-biosignal.r-universe.dev/PhysioAnalysis/reference/specparam.html).

## Usage

``` r
eegAperiodic(
  pe,
  freq_range = c(1, 45),
  aperiodic_mode = c("fixed", "knee"),
  max_n_peaks = 6L,
  peak_width_limits = c(1, 12),
  min_peak_height = 0.05,
  peak_threshold = 2,
  assay_name = NULL
)
```

## Arguments

- pe:

  A `PhysioExperiment`.

- freq_range:

  Frequency range to fit, in Hz (default `c(1, 45)`).

- aperiodic_mode:

  `"fixed"` (offset + exponent) or `"knee"` (offset + knee + exponent),
  the latter for spectra with a bend in log-log space.

- max_n_peaks:

  Maximum number of oscillatory peaks per channel (default 6).

- peak_width_limits:

  Min/max peak width in Hz (default `c(1, 12)`).

- min_peak_height:

  Minimum peak height above the aperiodic fit (default 0.05).

- peak_threshold:

  Peak detection threshold in SD of the flattened spectrum (default 2).

- assay_name:

  Assay to use (default: the object's default assay).

## Value

An `eeg_aperiodic` object: a list with `aperiodic` (per-channel data
frame: `channel`, `exponent`, `offset`, optionally `knee`, `r_squared`,
`error`), `peaks` (per-channel `CF`/`PW`/`BW`), `exponent` (a named
per-channel vector, e.g. for
[`eegPlotTopomap()`](https://x-biosignal.github.io/PhysioEEG/reference/eegPlotTopomap.md)),
and the underlying `specparam_result`.

## References

Donoghue et al. 2020, Nat Neurosci (specparam / FOOOF).

## See also

[`eegQEEG()`](https://x-biosignal.github.io/PhysioEEG/reference/eegQEEG.md),
[`eegComplexity()`](https://x-biosignal.github.io/PhysioEEG/reference/eegComplexity.md)

## Examples

``` r
if (FALSE) { # \dontrun{
pe <- make_eeg(n_time = 2500, n_channels = 8, sr = 250)
ap <- eegAperiodic(pe, freq_range = c(2, 40))
ap$aperiodic          # per-channel exponent / offset
} # }
```

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
if (requireNamespace("PhysioAnalysis", quietly = TRUE)) {
  pe <- make_eeg(n_time = 2500, n_channels = 8, sr = 250)
  ap <- eegAperiodic(pe, freq_range = c(2, 40))
  ap$aperiodic          # per-channel exponent / offset
}
#>   channel    offset exponent r_squared     error
#> 1     Fp1 1.0569143 1.899793 0.9517044 0.1794116
#> 2     Fp2 0.9976620 1.893552 0.9729035 0.1329855
#> 3      F7 0.7323077 1.689029 0.9699991 0.1272416
#> 4      F3 0.9299269 1.865744 0.9517095 0.1647729
#> 5      Fz 0.7507181 1.747007 0.9560722 0.1427750
#> 6      F4 0.9321072 1.853310 0.9602022 0.1416420
#> 7      F8 0.9810790 1.895442 0.9554800 0.1608447
#> 8      T3 1.0408772 1.930276 0.9532819 0.1750132
```

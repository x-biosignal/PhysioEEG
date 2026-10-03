# Phase-amplitude coupling (cross-frequency coupling) for EEG

Measures how the phase of a slow oscillation modulates the amplitude of
a fast oscillation — within a single channel (`amp_channel = NULL`) or
across two channels. Delegates the estimator to
[`PhysioCrossModal::phaseAmplitudeCoupling()`](https://x-biosignal.github.io/PhysioCrossModal/reference/phaseAmplitudeCoupling.html).

## Usage

``` r
eegPAC(
  pe,
  phase_channel = 1L,
  amp_channel = NULL,
  phase_band = c(4, 8),
  amp_band = c(30, 80),
  method = c("tort", "canolty", "ozkurt", "plv"),
  n_bins = 18L,
  assay_name = NULL
)
```

## Arguments

- pe:

  A `PhysioExperiment` (continuous data; epoched data is concatenated).

- phase_channel:

  Channel (index or label) supplying the modulating phase.

- amp_channel:

  Channel supplying the modulated amplitude; `NULL` (default) uses
  `phase_channel` (within-channel PAC).

- phase_band:

  Phase frequency band in Hz (default theta `c(4, 8)`).

- amp_band:

  Amplitude frequency band in Hz (default gamma `c(30, 80)`).

- method:

  `"tort"` (modulation index), `"canolty"` (mean vector length),
  `"ozkurt"`, or `"plv"`.

- n_bins:

  Number of phase bins for the Tort modulation index (default 18).

- assay_name:

  Assay to use (default: the object's default assay).

## Value

The PAC result from
[`PhysioCrossModal::phaseAmplitudeCoupling()`](https://x-biosignal.github.io/PhysioCrossModal/reference/phaseAmplitudeCoupling.html)
(a list with `pac` and, for `"tort"`, the phase-amplitude
`distribution`).

## References

Tort et al. 2010; Canolty et al. 2006.

## See also

[`eegComodulogram()`](https://x-biosignal.github.io/PhysioEEG/reference/eegComodulogram.md),
[`eegConnectivityMatrix()`](https://x-biosignal.github.io/PhysioEEG/reference/eegConnectivityMatrix.md)

## Examples

``` r
if (requireNamespace("PhysioCrossModal", quietly = TRUE)) {
  pe <- make_eeg(n_time = 5000, n_channels = 4, sr = 250)
  eegPAC(pe, phase_channel = 1, phase_band = c(4, 8), amp_band = c(30, 80))
}
#> $pac
#> [1] 0.0004274329
#> 
#> $method
#> [1] "tort"
#> 
#> $phase_band
#> [1] 4 8
#> 
#> $amp_band
#> [1] 30 80
#> 
#> $n_bins
#> [1] 18
#> 
#> $distribution
#>  [1] 0.06013918 0.05539763 0.05877512 0.05768208 0.05505245 0.05362357
#>  [7] 0.05160298 0.05318179 0.05622009 0.05821859 0.05694298 0.05811456
#> [13] 0.05636442 0.05273079 0.05119354 0.05119147 0.05455322 0.05901553
#> 
#> $bin_centers
#>  [1] -2.9670597 -2.6179939 -2.2689280 -1.9198622 -1.5707963 -1.2217305
#>  [7] -0.8726646 -0.5235988 -0.1745329  0.1745329  0.5235988  0.8726646
#> [13]  1.2217305  1.5707963  1.9198622  2.2689280  2.6179939  2.9670597
#> 
```

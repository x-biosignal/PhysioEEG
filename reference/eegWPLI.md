# EEG Weighted Phase Lag Index (wPLI)

Computes the Weighted Phase Lag Index and its debiased variant between
all EEG channel pairs. The wPLI reduces the influence of volume
conduction by weighting the phase differences by the magnitude of the
imaginary part of the cross-spectrum, computed via Welch's method.

## Usage

``` r
eegWPLI(
  x,
  band = c(8, 13),
  window_sec = 2,
  overlap = 0.5,
  assay_name = NULL,
  debiased = TRUE
)
```

## Arguments

- x:

  A PhysioExperiment object with EEG data (2D: time x channels).

- band:

  Numeric vector of length 2 specifying the frequency band in Hz over
  which to compute wPLI (default: `c(8, 13)` for alpha band).

- window_sec:

  Window length in seconds for spectral estimation (default: 2).

- overlap:

  Overlap fraction between adjacent windows, from 0 to 1 exclusive
  (default: 0.5).

- assay_name:

  Name of the input assay. If `NULL`, the default assay is used.

- debiased:

  Logical; also compute the unbiased debiased wPLI (default TRUE).

## Value

A data.frame with columns:

- channel1:

  Character or integer identifier of the first channel.

- channel2:

  Character or integer identifier of the second channel.

- wpli:

  Numeric wPLI value in `[0, 1]`.

- wpli_debiased:

  Numeric debiased wPLI (Vinck et al., 2011, Eq. 6), computed by the
  shared
  [`wpliEstimate`](https://x-biosignal.github.io/PhysioCore//reference/wpliEstimate.html)
  estimator so it matches
  [`PhysioCrossModal::weightedPLI`](https://x-biosignal.github.io/PhysioCrossModal/reference/weightedPLI.html).
  It is unbiased and therefore distributes around 0 for independent
  signals (it is *not* clamped to be non-negative).

## Details

\$\$wPLI = \|mean(Im(S\_{xy}))\| / mean(\|Im(S\_{xy})\|)\$\$

The debiased wPLI corrects for sample-size bias: \$\$wPLI^2\_{debiased}
= (N \cdot wPLI^2 - 1) / (N - 1)\$\$

## References

Lachaux, J. P., et al. (1999). Measuring phase synchrony in brain
signals. Human Brain Mapping, 8(4), 194-208.

Vinck, M., et al. (2011). An improved index of phase-synchronization for
electrophysiological data in the presence of volume-conduction, noise
and sample-size bias. NeuroImage, 55(4), 1548-1565.

## See also

[`eegCoherence()`](https://x-biosignal.github.io/PhysioEEG/reference/eegCoherence.md),
[`eegPLV()`](https://x-biosignal.github.io/PhysioEEG/reference/eegPLV.md),
[`eegConnectivityMatrix()`](https://x-biosignal.github.io/PhysioEEG/reference/eegConnectivityMatrix.md),
[`eegPlotConnectivity()`](https://x-biosignal.github.io/PhysioEEG/reference/eegPlotConnectivity.md),
[`wpliEstimate`](https://x-biosignal.github.io/PhysioCore//reference/wpliEstimate.html)

## Examples

``` r
if (FALSE) { # \dontrun{
pe <- make_eeg(n_time = 5000, n_channels = 4, sr = 500)
wpli_df <- eegWPLI(pe, band = c(8, 13))
head(wpli_df)
} # }
```

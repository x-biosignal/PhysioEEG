# EEG Spectral Granger Causality

Computes spectral Granger causality between all EEG channel pairs by
fitting bivariate autoregressive (AR) models using Yule-Walker equations
and computing the transfer function in the frequency domain. Unlike
coherence-based measures, Granger causality is directional: GC from
channel A to channel B is generally different from GC from B to A.

## Usage

``` r
eegGrangerCausality(x, order = 5, band = NULL, assay_name = NULL)
```

## Arguments

- x:

  A PhysioExperiment object with EEG data (2D: time x channels).

- order:

  Integer AR model order for the bivariate model (default: 5).

- band:

  Numeric vector of length 2 specifying the frequency band in Hz over
  which to average GC. If `NULL`, returns the average over all positive
  frequencies.

- assay_name:

  Name of the input assay. If `NULL`, the default assay is used.

## Value

A data.frame with columns:

- from_channel:

  Character or integer identifier of the source channel.

- to_channel:

  Character or integer identifier of the target channel.

- gc_value:

  Numeric Granger causality value (\>= 0). Higher values indicate
  stronger directed influence.

## Details

The spectral Granger causality from channel x to channel y at frequency
f is defined as: \$\$GC\_{x \to y}(f) = \log(S\_{yy}(f) / (S\_{yy}(f) -
\|H\_{xy}(f)\|^2 \cdot \Sigma\_{xx}))\$\$

## References

Lachaux, J. P., et al. (1999). Measuring phase synchrony in brain
signals. Human Brain Mapping, 8(4), 194-208.

Granger, C. W. J. (1969). Investigating causal relations by econometric
models and cross-spectral methods. Econometrica, 37(3), 424-438.

## See also

[`eegCoherence()`](https://x-biosignal.github.io/PhysioEEG/reference/eegCoherence.md),
[`eegPLV()`](https://x-biosignal.github.io/PhysioEEG/reference/eegPLV.md),
[`eegWPLI()`](https://x-biosignal.github.io/PhysioEEG/reference/eegWPLI.md),
[`eegConnectivityMatrix()`](https://x-biosignal.github.io/PhysioEEG/reference/eegConnectivityMatrix.md)

## Examples

``` r
if (FALSE) { # \dontrun{
pe <- make_eeg(n_time = 5000, n_channels = 4, sr = 500)
gc_df <- eegGrangerCausality(pe, order = 5, band = c(8, 13))
head(gc_df)
} # }
```

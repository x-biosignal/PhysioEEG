# EEG Connectivity Matrix

Convenience wrapper that computes a symmetric n_channels x n_channels
connectivity matrix using the specified method. Dispatches to the
appropriate connectivity function and returns a named matrix suitable
for visualization or graph-theoretic analysis.

## Usage

``` r
eegConnectivityMatrix(
  x,
  method = c("coherence", "plv", "wpli", "dtf", "pdc"),
  band = c(8, 13),
  ...
)
```

## Arguments

- x:

  A PhysioExperiment object with EEG data (2D: time x channels).

- method:

  Connectivity method: `"coherence"` for magnitude-squared coherence,
  `"plv"` for Phase Locking Value, `"wpli"` for weighted Phase Lag
  Index, or the directed measures `"dtf"`
  ([`eegDTF()`](https://x-biosignal.github.io/PhysioEEG/reference/eegDTF.md))
  and `"pdc"`
  ([`eegPDC()`](https://x-biosignal.github.io/PhysioEEG/reference/eegPDC.md))
  (default: `"coherence"`). For `"dtf"` and `"pdc"` the returned matrix
  is directed (rows = targets, columns = sources).

- band:

  Numeric vector of length 2 specifying the frequency band in Hz
  (default: `c(8, 13)` for alpha band).

- ...:

  Additional arguments passed to the underlying connectivity function
  (e.g., `window_sec`, `overlap`).

## Value

A named numeric matrix of dimension n_channels x n_channels with channel
labels as row and column names. Diagonal elements are 1. Off-diagonal
elements represent the band-averaged connectivity between channel pairs.

## References

Lachaux, J. P., et al. (1999). Measuring phase synchrony in brain
signals. Human Brain Mapping, 8(4), 194-208.

## See also

[`eegCoherence()`](https://x-biosignal.github.io/PhysioEEG/reference/eegCoherence.md),
[`eegPLV()`](https://x-biosignal.github.io/PhysioEEG/reference/eegPLV.md),
[`eegWPLI()`](https://x-biosignal.github.io/PhysioEEG/reference/eegWPLI.md),
[`eegGrangerCausality()`](https://x-biosignal.github.io/PhysioEEG/reference/eegGrangerCausality.md),
[`eegPlotConnectivity()`](https://x-biosignal.github.io/PhysioEEG/reference/eegPlotConnectivity.md)

## Examples

``` r
if (FALSE) { # \dontrun{
pe <- make_eeg(n_time = 5000, n_channels = 4, sr = 500)
conn <- eegConnectivityMatrix(pe, method = "plv", band = c(8, 13))
print(conn)
} # }
```

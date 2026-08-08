# EEG Phase Locking Value (PLV)

Computes the Phase Locking Value between all EEG channel pairs within a
specified frequency band. Each channel is bandpass filtered,
instantaneous phase is extracted via the Hilbert transform, and PLV is
computed as the mean resultant length of the phase difference
distribution.

## Usage

``` r
eegPLV(x, band = c(8, 13), assay_name = NULL)
```

## Arguments

- x:

  A PhysioExperiment object with EEG data (2D: time x channels).

- band:

  Numeric vector of length 2 specifying the frequency band in Hz for
  bandpass filtering before phase extraction (default: `c(8, 13)` for
  alpha band).

- assay_name:

  Name of the input assay. If `NULL`, the default assay is used.

## Value

A data.frame with columns:

- channel1:

  Character or integer identifier of the first channel.

- channel2:

  Character or integer identifier of the second channel.

- plv:

  Numeric PLV value in `[0, 1]`.

## Details

\$\$PLV = \|1/N \sum\_{t=1}^{N} exp(i(\phi_x(t) - \phi_y(t)))\|\$\$

## References

Lachaux, J. P., et al. (1999). Measuring phase synchrony in brain
signals. Human Brain Mapping, 8(4), 194-208.

## See also

[`eegCoherence()`](https://x-biosignal.github.io/PhysioEEG/reference/eegCoherence.md),
[`eegWPLI()`](https://x-biosignal.github.io/PhysioEEG/reference/eegWPLI.md),
[`eegConnectivityMatrix()`](https://x-biosignal.github.io/PhysioEEG/reference/eegConnectivityMatrix.md),
[`eegPlotConnectivity()`](https://x-biosignal.github.io/PhysioEEG/reference/eegPlotConnectivity.md)

## Examples

``` r
if (FALSE) { # \dontrun{
pe <- make_eeg(n_time = 5000, n_channels = 4, sr = 500)
plv_df <- eegPLV(pe, band = c(8, 13))
head(plv_df)
} # }
```

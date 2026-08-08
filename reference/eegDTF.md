# Directed Transfer Function (DTF)

Estimates frequency-resolved directed connectivity between EEG channels
with the Directed Transfer Function (Kaminski & Blinowska 1991),
computed from the transfer function \\H(f)\\ of a multivariate
autoregressive (MVAR) model fitted with
[`PhysioCore::mvarFit()`](https://x-biosignal.r-universe.dev/PhysioCore/reference/mvarFit.html).
Because \\H(f)\\ inverts the whole VAR system, DTF reflects both direct
and indirect pathways. The (default) normalized DTF satisfies \\\sum_j
\mathrm{DTF}\_{ij}(f)^2 = 1\\ for each target \\i\\ (inflow
normalization); `ffDTF` normalizes across all frequencies instead
(Korzeniewska 2003).

## Usage

``` r
eegDTF(
  x,
  order = NULL,
  freqs = NULL,
  normalized = TRUE,
  ffDTF = FALSE,
  band = NULL,
  method = "ols",
  assay_name = NULL
)
```

## Arguments

- x:

  A PhysioExperiment object with 2D EEG data (time x channels).

- order:

  MVAR model order, or `NULL` to select automatically.

- freqs:

  Numeric vector of frequencies in Hz (default: 128 points from 0 to the
  Nyquist frequency).

- normalized:

  Row-normalize the DTF (default: `TRUE`).

- ffDTF:

  Use the full-frequency DTF normalization (default: `FALSE`).

- band:

  Optional numeric length-2 band in Hz over which to average the stored
  connectivity matrix (default: all frequencies).

- method:

  MVAR estimator passed to
  [`PhysioCore::mvarFit()`](https://x-biosignal.r-universe.dev/PhysioCore/reference/mvarFit.html)
  (default: `"ols"`).

- assay_name:

  Input assay name (default: the default assay).

## Value

The PhysioExperiment with `metadata(x)$connectivity` set to a list
containing the band-averaged directed `matrix` (rows = targets, columns
= sources), the frequency-resolved `array` (n_channels x n_channels x
n_freqs), the `frequencies`, and the settings used.

## References

Kaminski, M. J., & Blinowska, K. J. (1991). A new method of the
description of the information flow in the brain structures. Biological
Cybernetics, 65(3), 203-210.

## See also

[`eegPDC()`](https://x-biosignal.github.io/PhysioEEG/reference/eegPDC.md),
[`eegConditionalGC()`](https://x-biosignal.github.io/PhysioEEG/reference/eegConditionalGC.md),
[`eegConnectivityMatrix()`](https://x-biosignal.github.io/PhysioEEG/reference/eegConnectivityMatrix.md),
[`PhysioCore::mvarFit()`](https://x-biosignal.r-universe.dev/PhysioCore/reference/mvarFit.html)

## Examples

``` r
if (FALSE) { # \dontrun{
pe <- make_eeg(n_time = 4000, n_channels = 5, sr = 250)
pe <- eegDTF(pe, order = 5)
dim(metadata(pe)$connectivity$array)
} # }
```

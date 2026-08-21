# Partial Directed Coherence (PDC)

Estimates frequency-resolved directed connectivity with Partial Directed
Coherence (Baccala & Sameshima 2001), computed from the frequency-domain
coefficient matrix \\\bar{A}(f)\\ of an MVAR model
([`PhysioCore::mvarFit()`](https://x-biosignal.github.io/PhysioCore//reference/mvarFit.html)).
Unlike DTF, PDC reflects only *direct* channel-to-channel influences, so
a purely indirect pathway gives PDC near zero. The (default) generalized
PDC weights each row by the inverse residual standard deviation to make
the measure scale-invariant (Baccala 2007). PDC satisfies \\\sum_i
\mathrm{PDC}\_{ij}(f)^2 = 1\\ for each source \\j\\ (outflow
normalization).

## Usage

``` r
eegPDC(
  x,
  order = NULL,
  freqs = NULL,
  generalized = TRUE,
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

- generalized:

  Use generalized PDC (default: `TRUE`).

- band:

  Optional numeric length-2 band in Hz over which to average the stored
  connectivity matrix (default: all frequencies).

- method:

  MVAR estimator passed to
  [`PhysioCore::mvarFit()`](https://x-biosignal.github.io/PhysioCore//reference/mvarFit.html)
  (default: `"ols"`).

- assay_name:

  Input assay name (default: the default assay).

## Value

The PhysioExperiment with `metadata(x)$connectivity` set as in
[`eegDTF()`](https://x-biosignal.github.io/PhysioEEG/reference/eegDTF.md)
(a band-averaged directed `matrix`, the frequency-resolved `array`,
`frequencies`, and settings).

## References

Baccala, L. A., & Sameshima, K. (2001). Partial directed coherence: a
new concept in neural structure determination. Biological Cybernetics,
84(6), 463-474.

## See also

[`eegDTF()`](https://x-biosignal.github.io/PhysioEEG/reference/eegDTF.md),
[`eegConditionalGC()`](https://x-biosignal.github.io/PhysioEEG/reference/eegConditionalGC.md),
[`eegConnectivityMatrix()`](https://x-biosignal.github.io/PhysioEEG/reference/eegConnectivityMatrix.md)

## Examples

``` r
if (FALSE) { # \dontrun{
pe <- make_eeg(n_time = 4000, n_channels = 5, sr = 250)
pe <- eegPDC(pe, order = 5)
metadata(pe)$connectivity$matrix
} # }
```

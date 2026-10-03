# Partial Directed Coherence (PDC)

Estimates frequency-resolved directed connectivity with Partial Directed
Coherence (Baccala & Sameshima 2001), computed from the frequency-domain
coefficient matrix \\\bar{A}(f)\\ of an MVAR model
([`PhysioExperiment::mvarFit()`](https://x-biosignal.r-universe.dev/PhysioExperiment/reference/mvarFit.html)).
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
  [`PhysioExperiment::mvarFit()`](https://x-biosignal.r-universe.dev/PhysioExperiment/reference/mvarFit.html)
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
pe <- make_eeg(n_time = 4000, n_channels = 5, sr = 250)
pe <- eegPDC(pe, order = 5)
S4Vectors::metadata(pe)$connectivity$matrix
#>           Fp1        Fp2         F7         F3        Fz
#> Fp1 0.9308886 0.20475110 0.09796020 0.10207939 0.1918640
#> Fp2 0.1188708 0.95332157 0.17014541 0.08760945 0.1391627
#> F7  0.1051502 0.12406314 0.96623214 0.07591603 0.1297131
#> F3  0.1390672 0.04612233 0.04757274 0.96540989 0.1932865
#> Fz  0.2153827 0.09414714 0.06520345 0.16158563 0.9137596
```

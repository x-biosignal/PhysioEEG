# Delta/Alpha Ratio (DAR)

Computes the delta-to-alpha band-power ratio per channel (or per
region), a quantitative-EEG stroke marker whose elevation tracks
ischaemic slowing (Finnigan et al. 2016). Band powers use a Welch
estimator for a stable spectrum. Each ratio is returned as a
reliability-characterised
[`physioBiomarker`](https://x-biosignal.github.io/PhysioCore//reference/physioBiomarker-constructor.html)
carrying provenance (band, method, software version), reliability
placeholders, and - when age (and montage) are supplied - a published
reference range from
[`normativeLookup`](https://x-biosignal.github.io/PhysioCore//reference/normativeLookup.html)
so the value prints its normative percentile.

## Usage

``` r
eegDAR(
  x,
  delta = c(1, 4),
  alpha = c(8, 13),
  regions = NULL,
  method = "welch",
  age = NULL,
  montage = NULL,
  reliability = list(icc = NA_real_, sem = NA_real_),
  assay_name = NULL
)
```

## Arguments

- x:

  A PhysioExperiment object with 2D EEG data.

- delta:

  Numeric length-2 delta band in Hz (default: `c(1, 4)`).

- alpha:

  Numeric length-2 alpha band in Hz (default: `c(8, 13)`).

- regions:

  Optional named list mapping a region name to a character vector of
  channel labels; band powers are averaged within each region. When
  `NULL` (default) one biomarker is returned per channel.

- method:

  Spectral estimator passed to the band-power routine (default:
  `"welch"`).

- age:

  Optional numeric age used to select the normative reference row.

- montage:

  Optional character montage for the normative lookup.

- reliability:

  Named list of reliability indices to attach (default: `icc` and `sem`
  placeholders).

- assay_name:

  Input assay name (default: the default assay).

## Value

A named list of
[`PhysioBiomarker`](https://x-biosignal.github.io/PhysioCore//reference/PhysioBiomarker.html)
objects, one per channel or region.

## References

Finnigan, S., Wong, A., & Read, S. (2016). Defining abnormal slow EEG
activity in acute ischaemic stroke: Delta/alpha ratio as an optimal QEEG
index. Clinical Neurophysiology, 127(2), 1452-1459.

## See also

[`eegBSI()`](https://x-biosignal.github.io/PhysioEEG/reference/eegBSI.md),
[`eegSlowing()`](https://x-biosignal.github.io/PhysioEEG/reference/eegSlowing.md),
[`eegQEEG()`](https://x-biosignal.github.io/PhysioEEG/reference/eegQEEG.md),
[`PhysioCore::physioBiomarker()`](https://x-biosignal.github.io/PhysioCore//reference/physioBiomarker-constructor.html),
[`PhysioCore::normativeLookup()`](https://x-biosignal.github.io/PhysioCore//reference/normativeLookup.html)

## Examples

``` r
if (FALSE) { # \dontrun{
pe <- make_eeg(n_time = 5000, n_channels = 19, sr = 250)
dar <- eegDAR(pe, age = 60)
dar[["C3"]]
} # }
```

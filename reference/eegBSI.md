# Brain Symmetry Index (BSI) and pairwise-derived pdBSI

Quantifies left/right spectral asymmetry across hemispheric electrode
pairs, a quantitative-EEG stroke marker (van Putten 2007; Sheorajpanday
et al. 2011). The global `bsi` is the spectrum-averaged absolute
relative difference between summed right- and left-hemisphere power,
bounded from 0 (symmetric) to 1 (one hemisphere silent). The
pairwise-derived `pdbsi` averages the same relative difference per
hemispheric pair; with `directed = TRUE` it keeps the sign (positive =
right-dominant). Both are returned as reliability-characterised
[`physioBiomarker`](https://x-biosignal.r-universe.dev/PhysioCore/reference/physioBiomarker-constructor.html)
objects with provenance and, for the global BSI, a normative reference
range.

## Usage

``` r
eegBSI(
  x,
  pairs = NULL,
  band = c(1, 25),
  directed = FALSE,
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

- pairs:

  List of length-2 character vectors, each a `c(right, left)`
  hemispheric electrode pair. Defaults to the standard 10-20 hemispheric
  pairs.

- band:

  Numeric length-2 frequency band in Hz (default: `c(1, 25)`).

- directed:

  Logical; when `TRUE` the pairwise index keeps its sign (default:
  `FALSE`).

- method:

  Spectral estimator (default: `"welch"`).

- age:

  Optional numeric age for the normative lookup.

- montage:

  Optional character montage for the normative lookup.

- reliability:

  Named list of reliability indices to attach (default: `icc` and `sem`
  placeholders).

- assay_name:

  Input assay name (default: the default assay).

## Value

A named list with two
[`PhysioBiomarker`](https://x-biosignal.r-universe.dev/PhysioCore/reference/PhysioBiomarker.html)
objects: `bsi` (global, from 0 to 1) and `pdbsi` (pairwise-derived,
signed when `directed = TRUE`).

## References

van Putten, M. J. A. M. (2007). The revised brain symmetry index.
Clinical Neurophysiology, 118(11), 2362-2367.

Sheorajpanday, R. V. A., et al. (2011). Quantitative EEG in ischemic
stroke: correlation with functional status after 6 months. Clinical
Neurophysiology, 122(5), 874-883.

## See also

[`eegDAR()`](https://x-biosignal.github.io/PhysioEEG/reference/eegDAR.md),
[`eegAsymmetry()`](https://x-biosignal.github.io/PhysioEEG/reference/eegAsymmetry.md),
[`PhysioCore::physioBiomarker()`](https://x-biosignal.r-universe.dev/PhysioCore/reference/physioBiomarker-constructor.html)

## Examples

``` r
if (FALSE) { # \dontrun{
pe <- make_eeg(n_time = 5000, n_channels = 19, sr = 250)
res <- eegBSI(pe, age = 60)
res$bsi
} # }
```

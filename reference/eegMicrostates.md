# EEG Microstate Segmentation

Performs microstate analysis on EEG data by identifying dominant scalp
topographies at Global Field Power (GFP) peaks and assigning each time
point to the best-matching microstate map. Supports polarity-invariant
K-means, atomize-and-agglomerate hierarchical clustering (AAHC), and
PCA-based extraction.

## Usage

``` r
eegMicrostates(
  x,
  n_states = 4,
  method = c("kmeans", "aahc", "pca"),
  min_gfp = 1,
  assay_name = NULL
)
```

## Arguments

- x:

  A PhysioExperiment object with EEG data.

- n_states:

  Number of microstate classes to extract (default: 4).

- method:

  Clustering method: `"kmeans"` (polarity-invariant K-means), `"aahc"`
  (atomize and agglomerate hierarchical clustering), or `"pca"`
  (principal component analysis).

- min_gfp:

  Percentile threshold (0-100) for GFP peak selection (default: 1.0).
  Only GFP peaks above this percentile are used for clustering.

- assay_name:

  Name of the input assay. If `NULL`, the default assay is used.

## Value

Modified PhysioExperiment with microstate results stored in
`metadata(x)$microstates`, a list containing:

- maps:

  Numeric matrix of dimensions n_channels x n_states, each column a
  microstate topography.

- labels:

  Integer vector of length n_time, microstate assignment (1 to n_states)
  for each time point.

- gfp:

  Numeric vector of GFP values per time point.

- n_states:

  Integer number of microstate classes.

## References

Michel, C. M., & Koenig, T. (2018). EEG microstates as a tool for
studying the temporal dynamics of whole-brain neuronal networks.
NeuroImage, 180, 577-593.

## See also

[`eegMicrostateStats()`](https://x-biosignal.github.io/PhysioEEG/reference/eegMicrostateStats.md),
[`eegMicrostateBackfit()`](https://x-biosignal.github.io/PhysioEEG/reference/eegMicrostateBackfit.md),
[`eegMicrostateSequence()`](https://x-biosignal.github.io/PhysioEEG/reference/eegMicrostateSequence.md)

## Examples

``` r
if (FALSE) { # \dontrun{
pe <- make_eeg(n_time = 5000, n_channels = 19, sr = 500)
pe <- eegMicrostates(pe, n_states = 4, method = "kmeans")
ms <- metadata(pe)$microstates
} # }
```

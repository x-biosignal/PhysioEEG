# Re-reference EEG data

Applies a re-referencing scheme to EEG data. Re-referencing transforms
the data by subtracting a reference signal from all channels, which can
improve spatial resolution and comparability across studies.

## Usage

``` r
eegRereference(
  x,
  ref_type = c("average", "robust", "median", "mastoids", "cz", "rest", "channel"),
  ref_channels = NULL,
  exclude = NULL,
  assay_name = NULL,
  output_assay = "rereferenced",
  robust_noise_sd = 4,
  robust_max_iter = 5
)
```

## Arguments

- x:

  A PhysioExperiment object.

- ref_type:

  Re-referencing scheme: `"average"` (common average), `"robust"` (PREP
  robust average reference: iteratively detect bad channels and exclude
  them from the average, Bigdely-Shamlo et al., 2015), `"median"`
  (channel-wise median reference, robust to outlier channels),
  `"mastoids"` (linked mastoids), `"cz"` (Cz reference), `"rest"`
  (Reference Electrode Standardization Technique), or `"channel"`
  (user-specified channels).

- ref_channels:

  Character vector of channel labels to use as reference (required for
  `ref_type = "channel"`).

- exclude:

  Character vector of channel labels to exclude from the average
  reference calculation (only used for `ref_type = "average"`).

- assay_name:

  Name of the assay to re-reference. If NULL, uses `defaultAssay(x)`.

- output_assay:

  Name of the output assay (default: `"rereferenced"`).

- robust_noise_sd:

  Noise threshold (SDs above the median channel variance) for
  bad-channel detection in `ref_type = "robust"` (default: 4).

- robust_max_iter:

  Maximum PREP iterations for `ref_type = "robust"` (default: 5).

## Value

A PhysioExperiment object with re-referenced data in the specified
output assay. Stores reference info (including the channels excluded by
the robust reference) in `metadata(x)$reference` and logs the step in
the object's provenance.

## References

Bigdely-Shamlo, N., et al. (2015). "The PREP pipeline: standardized
preprocessing for large-scale EEG analysis." *Frontiers in
Neuroinformatics*, 9, 16.
[doi:10.3389/fninf.2015.00016](https://doi.org/10.3389/fninf.2015.00016)

## Examples

``` r
if (FALSE) { # \dontrun{
pe <- make_eeg(n_time = 5000, n_channels = 19, sr = 500)
pe_avg    <- eegRereference(pe, ref_type = "average")
pe_robust <- eegRereference(pe, ref_type = "robust")
pe_median <- eegRereference(pe, ref_type = "median")
} # }
```

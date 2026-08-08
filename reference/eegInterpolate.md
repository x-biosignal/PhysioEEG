# Interpolate bad EEG channels

Replaces data in bad channels by interpolating from remaining good
channels using either spherical spline interpolation (Perrin et al.,
1989) or nearest-neighbor weighted averaging.

## Usage

``` r
eegInterpolate(
  x,
  bad_channels,
  method = c("spline", "nearest"),
  assay_name = NULL,
  output_assay = "interpolated"
)
```

## Arguments

- x:

  A PhysioExperiment object.

- bad_channels:

  Character vector of channel labels to interpolate.

- method:

  Interpolation method: `"spline"` for spherical spline (default) or
  `"nearest"` for inverse-distance weighted nearest neighbors.

- assay_name:

  Name of the assay to interpolate. If NULL, uses `defaultAssay(x)`.

- output_assay:

  Name of the output assay (default: `"interpolated"`).

## Value

A PhysioExperiment object with interpolated channels in the specified
output assay.

## Details

Requires electrode positions (pos_x, pos_y, pos_z) in colData. Apply
[`eegMontage`](https://x-biosignal.github.io/PhysioEEG/reference/eegMontage.md)
first if positions are not set.

## Examples

``` r
if (FALSE) { # \dontrun{
pe <- make_eeg(n_time = 5000, n_channels = 19, sr = 500)
pe <- eegMontage(pe, system = "10-20")
bad_df <- eegBadChannels(pe)
bad_labels <- bad_df$channel[bad_df$is_bad]
if (length(bad_labels) > 0) {
  pe_clean <- eegInterpolate(pe, bad_labels)
}
} # }
```

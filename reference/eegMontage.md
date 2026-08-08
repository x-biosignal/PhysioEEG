# Assign electrode montage positions

Maps channel labels to standard electrode positions from a known montage
system (10-20, 10-10, BioSemi64) or from a custom positions data.frame.
Sets `pos_x`, `pos_y`, `pos_z` columns in colData for use by
interpolation, topographic mapping, and source localization functions.

## Usage

``` r
eegMontage(
  x,
  system = c("10-20", "10-10", "biosemi64", "custom"),
  positions = NULL
)
```

## Arguments

- x:

  A PhysioExperiment object.

- system:

  Montage system: `"10-20"` (19 channels), `"10-10"` (~64 channels),
  `"biosemi64"` (64 BioSemi channels), or `"custom"` (user-provided
  positions).

- positions:

  For `system = "custom"`, a data.frame with columns `label`, `pos_x`,
  `pos_y`, `pos_z`.

## Value

A PhysioExperiment object with updated colData containing position
columns.

## Examples

``` r
if (FALSE) { # \dontrun{
pe <- make_eeg(n_time = 5000, n_channels = 19, sr = 500)
pe <- eegMontage(pe, system = "10-20")
head(SummarizedExperiment::colData(pe))
} # }
```

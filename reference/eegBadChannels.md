# Detect bad EEG channels

Identifies bad (noisy, flat, or poorly correlated) EEG channels using
multiple automated criteria. Channels flagged as bad can subsequently be
interpolated using
[`eegInterpolate`](https://x-biosignal.github.io/PhysioEEG/reference/eegInterpolate.md).

## Usage

``` r
eegBadChannels(
  x,
  method = c("all", "flat", "noise", "correlation"),
  flat_threshold = 1e-06,
  noise_threshold = 4,
  corr_threshold = 0.4,
  assay_name = NULL
)
```

## Arguments

- x:

  A PhysioExperiment object.

- method:

  Detection method(s) to apply: `"all"` runs all checks, or specify one
  or more of `"flat"`, `"noise"`, `"correlation"`.

- flat_threshold:

  Variance threshold below which a channel is considered flat (default:
  1e-6).

- noise_threshold:

  Number of standard deviations above median variance to flag a channel
  as noisy (default: 4).

- corr_threshold:

  Minimum mean correlation with other channels. Channels below this are
  flagged (default: 0.4).

- assay_name:

  Name of the assay to analyze. If NULL, uses `defaultAssay(x)`.

## Value

A data.frame with columns: `channel` (label), `is_bad` (logical),
`reason` (character description), `score` (numeric metric value).

## Examples

``` r
if (FALSE) { # \dontrun{
pe <- make_eeg(n_time = 5000, n_channels = 19, sr = 500)
bad_df <- eegBadChannels(pe)
bad_labels <- bad_df$channel[bad_df$is_bad]
} # }
```

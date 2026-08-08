# Plot Sleep Hypnogram

Displays a sleep hypnogram showing sleep stage transitions over time.
Stages are ordered with Wake at the top and N3 at the bottom, with a
characteristic staircase pattern.

## Usage

``` r
eegPlotHypnogram(x, stages = NULL, epoch_sec = 30, colors = NULL)
```

## Arguments

- x:

  A PhysioExperiment object (used for metadata access).

- stages:

  A data.frame with columns `epoch` and `stage`, or `NULL` to read from
  `metadata(x)$sleep_stages`.

- epoch_sec:

  Epoch duration in seconds (default: 30).

- colors:

  Named character vector of colors for each stage. If `NULL`, default
  colors are used.

## Value

A ggplot2 object.

## Examples

``` r
if (FALSE) { # \dontrun{
pe <- make_eeg_sleep(n_time = 150000, n_channels = 2, sr = 500)
stages <- eegSleepStage(pe)
metadata(pe)$sleep_stages <- stages
eegPlotHypnogram(pe)
} # }
```

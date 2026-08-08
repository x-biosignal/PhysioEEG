# Epoch continuous EEG data

Segments continuous 2D EEG data into fixed-length epochs around events,
producing a 3D array (time x channels x epochs). Optionally performs
baseline correction by subtracting the mean of a pre-stimulus window.

## Usage

``` r
eegEpoch(
  x,
  events,
  limits = c(-0.2, 0.8),
  baseline = c(-0.2, 0),
  assay_name = NULL,
  output_assay = "epoched"
)
```

## Arguments

- x:

  A PhysioExperiment object with continuous (2D) data.

- events:

  A data.frame with an `onset_sec` column specifying event times in
  seconds, OR an integer vector of sample indices.

- limits:

  Numeric vector of length 2 specifying the epoch window relative to
  each event in seconds (default: `c(-0.2, 0.8)`).

- baseline:

  Numeric vector of length 2 specifying the baseline window relative to
  each event in seconds (default: `c(-0.2, 0)`). Set to NULL to skip
  baseline correction.

- assay_name:

  Name of the assay to epoch. If NULL, uses `defaultAssay(x)`.

- output_assay:

  Name of the output assay (default: `"epoched"`).

## Value

A PhysioExperiment object with a 3D array (time x channels x epochs) in
the specified output assay. Event information is stored in
`metadata(x)$epoch_events`.

## Examples

``` r
if (FALSE) { # \dontrun{
pe <- make_eeg(n_time = 10000, n_channels = 19, sr = 500)
events <- data.frame(onset_sec = c(1.0, 3.0, 5.0, 7.0))
pe_ep <- eegEpoch(pe, events, limits = c(-0.2, 0.8))
dim(SummarizedExperiment::assay(pe_ep, "epoched"))
# time x channels x epochs
} # }
```

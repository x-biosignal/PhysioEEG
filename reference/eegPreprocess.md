# Full EEG preprocessing pipeline

Convenience wrapper that runs a complete EEG preprocessing pipeline in
sequence: filtering, re-referencing, bad channel detection and
interpolation, optional ICA, epoching, and artifact rejection. Each step
can be toggled on or off.

## Usage

``` r
eegPreprocess(
  x,
  filter = TRUE,
  lowcut = 0.1,
  highcut = 40,
  notch = NULL,
  rereference = TRUE,
  ref_type = "average",
  bad_channels = TRUE,
  interpolate = TRUE,
  ica = FALSE,
  epoch = FALSE,
  events = NULL,
  epoch_limits = c(-0.2, 0.8),
  baseline = c(-0.2, 0),
  artifact_reject = FALSE,
  threshold_uv = 100,
  assay_name = NULL,
  verbose = TRUE
)
```

## Arguments

- x:

  A PhysioExperiment object with continuous (2D) EEG data.

- filter:

  Logical; apply frequency filtering (default: TRUE).

- lowcut:

  Low cutoff frequency in Hz for filtering (default: 0.1).

- highcut:

  High cutoff frequency in Hz for filtering (default: 40).

- notch:

  Notch filter center frequency in Hz (NULL for none).

- rereference:

  Logical; apply re-referencing (default: TRUE).

- ref_type:

  Re-referencing type (default: `"average"`).

- bad_channels:

  Logical; detect and report bad channels (default: TRUE).

- interpolate:

  Logical; interpolate detected bad channels (default: TRUE). Requires
  electrode positions in colData (apply `eegMontage` first or this will
  be skipped with a warning).

- ica:

  Logical; apply ICA-based artifact removal (default: FALSE). Requires
  the `eegICA` function to be available.

- epoch:

  Logical; epoch the data around events (default: FALSE).

- events:

  Event data for epoching (data.frame or integer vector). Required if
  `epoch = TRUE`.

- epoch_limits:

  Epoch window in seconds (default: `c(-0.2, 0.8)`).

- baseline:

  Baseline window in seconds (default: `c(-0.2, 0)`).

- artifact_reject:

  Logical; reject artifact epochs (default: FALSE). Only applies if
  `epoch = TRUE`.

- threshold_uv:

  Amplitude threshold for artifact rejection (default: 100).

- assay_name:

  Starting assay name. If NULL, uses `defaultAssay(x)`.

- verbose:

  Logical; print progress messages (default: TRUE).

## Value

A PhysioExperiment object with processed data. The final assay name
depends on which steps are enabled. Processing log is stored in
`metadata(x)$preprocess_log`.

## Examples

``` r
if (FALSE) { # \dontrun{
pe <- make_eeg(n_time = 10000, n_channels = 19, sr = 500)
pe <- eegMontage(pe, system = "10-20")
events <- data.frame(onset_sec = c(1, 3, 5, 7, 9))
pe_proc <- eegPreprocess(pe, lowcut = 1, highcut = 40, notch = 50,
                          epoch = TRUE, events = events,
                          artifact_reject = TRUE)
} # }
```

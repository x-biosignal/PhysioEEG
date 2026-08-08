# Reject artifacts in epoched EEG data

Identifies and removes epochs contaminated by artifacts from 3D epoched
EEG data. Supports multiple detection criteria: amplitude threshold,
gradient (point-to-point voltage change), and joint probability.

## Usage

``` r
eegArtifactReject(
  x,
  method = c("threshold", "gradient", "joint_probability"),
  threshold_uv = 100,
  gradient_uv_ms = 50,
  jp_threshold = 3,
  assay_name = NULL,
  output_assay = "clean"
)
```

## Arguments

- x:

  A PhysioExperiment object with 3D epoched data.

- method:

  Artifact detection method: `"threshold"` (amplitude), `"gradient"`
  (point-to-point change), or `"joint_probability"` (log-power z-score).

- threshold_uv:

  Maximum absolute amplitude in microvolts for threshold rejection
  (default: 100).

- gradient_uv_ms:

  Maximum point-to-point change in uV/ms for gradient rejection
  (default: 50).

- jp_threshold:

  Number of standard deviations for joint probability rejection
  (default: 3).

- assay_name:

  Name of the assay to check. If NULL, uses `defaultAssay(x)`.

- output_assay:

  Name of the output assay (default: `"clean"`).

## Value

A PhysioExperiment object with clean epochs in the specified output
assay. Artifact log stored in `metadata(x)$artifact_log`.

## Examples

``` r
if (FALSE) { # \dontrun{
pe <- make_eeg(n_time = 10000, n_channels = 19, sr = 500)
events <- data.frame(onset_sec = c(1, 3, 5, 7, 9))
pe_ep <- eegEpoch(pe, events, limits = c(-0.2, 0.8))
pe_clean <- eegArtifactReject(pe_ep, method = "threshold",
                               threshold_uv = 100, assay_name = "epoched")
} # }
```

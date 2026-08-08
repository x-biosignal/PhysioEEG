# Create EEG with Motor Imagery / SSVEP Patterns for BCI

Generates epoched EEG with two classes: left motor imagery (mu ERD over
C4) and right motor imagery (mu ERD over C3). Also includes SSVEP at
specified frequencies at occipital channels. Data is stored as a 3D
array (time x channels x trials).

## Usage

``` r
make_eeg_bci(n_trials = 30, n_channels = 8, sr = 256, trial_sec = 4)
```

## Arguments

- n_trials:

  Number of trials per class (default: 30). Total trials will be
  `2 * n_trials`.

- n_channels:

  Number of channels (default: 8).

- sr:

  Sampling rate in Hz (default: 256).

- trial_sec:

  Trial duration in seconds (default: 4.0).

## Value

A PhysioExperiment with a 3D `"raw"` assay (time x channels x
total_trials). Sets `metadata(x)$labels` with a character vector of
class labels (`"left"` or `"right"`) for each trial.

## See also

[`make_eeg()`](https://x-biosignal.github.io/PhysioEEG/reference/make_eeg.md),
[`eegCSP()`](https://x-biosignal.github.io/PhysioEEG/reference/eegCSP.md),
[`eegBCIfeatures()`](https://x-biosignal.github.io/PhysioEEG/reference/eegBCIfeatures.md),
[`eegBCIclassify()`](https://x-biosignal.github.io/PhysioEEG/reference/eegBCIclassify.md),
[`eegMotorImagery()`](https://x-biosignal.github.io/PhysioEEG/reference/eegMotorImagery.md)

## Examples

``` r
pe <- make_eeg_bci(n_trials = 5, n_channels = 4, sr = 128)
dim(SummarizedExperiment::assay(pe, "raw"))  # 512 x 4 x 10
#> [1] 512   4  10
```

# Create EEG with Embedded ERP Components

Generates epoched EEG (3D array: time x channels x epochs) with known
ERP components: N100 (negative peak at ~100ms) and P300 (positive peak
at ~300ms). Target epochs have larger P300 amplitude than standard
epochs.

## Usage

``` r
make_eeg_erp(n_epochs = 40, n_channels = 19, sr = 250, epoch_sec = 1)
```

## Arguments

- n_epochs:

  Number of epochs (default: 40).

- n_channels:

  Number of channels (default: 19).

- sr:

  Sampling rate in Hz (default: 250).

- epoch_sec:

  Epoch duration in seconds (default: 1.0).

## Value

A PhysioExperiment with a 3D `"raw"` assay (time x channels x epochs)
containing simulated ERP data. Sets `metadata(x)$conditions` with a
character vector of epoch conditions (`"target"` or `"standard"`).

## See also

[`make_eeg()`](https://x-biosignal.github.io/PhysioEEG/reference/make_eeg.md),
[`eegERPdetect()`](https://x-biosignal.github.io/PhysioEEG/reference/eegERPdetect.md),
[`eegERPmeasure()`](https://x-biosignal.github.io/PhysioEEG/reference/eegERPmeasure.md),
[`eegERPbaseline()`](https://x-biosignal.github.io/PhysioEEG/reference/eegERPbaseline.md),
[`eegERPtest()`](https://x-biosignal.github.io/PhysioEEG/reference/eegERPtest.md)

## Examples

``` r
pe <- make_eeg_erp(n_epochs = 10, n_channels = 4, sr = 250)
dim(SummarizedExperiment::assay(pe, "raw"))  # 250 x 4 x 10
#> [1] 250   4  10
```

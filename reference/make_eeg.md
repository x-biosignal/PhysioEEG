# Create Simulated EEG Data

Generates multi-channel EEG with alpha (10 Hz) and beta (20 Hz)
oscillations plus pink noise. Channels are labeled using the
International 10-20 system.

## Usage

``` r
make_eeg(n_time = 5000, n_channels = 19, sr = 500)
```

## Arguments

- n_time:

  Number of time points (default: 5000 = 10s at 500 Hz).

- n_channels:

  Number of EEG channels (default: 19, standard 10-20).

- sr:

  Sampling rate in Hz (default: 500).

## Value

A PhysioExperiment with simulated EEG in the `"raw"` assay. Channel
labels follow the International 10-20 system (Fp1, Fp2, F7, ..., O1,
O2). Column data contains `label` and `type` fields.

## See also

[`make_eeg_erp()`](https://x-biosignal.github.io/PhysioEEG/reference/make_eeg_erp.md),
[`make_eeg_sleep()`](https://x-biosignal.github.io/PhysioEEG/reference/make_eeg_sleep.md),
[`make_eeg_bci()`](https://x-biosignal.github.io/PhysioEEG/reference/make_eeg_bci.md),
[`make_eeg_spikes()`](https://x-biosignal.github.io/PhysioEEG/reference/make_eeg_spikes.md),
[`eegFilter()`](https://x-biosignal.github.io/PhysioEEG/reference/eegFilter.md),
[`eegCoherence()`](https://x-biosignal.github.io/PhysioEEG/reference/eegCoherence.md)

## Examples

``` r
pe <- make_eeg(n_time = 2500, n_channels = 4, sr = 250)
dim(SummarizedExperiment::assay(pe, "raw"))  # 2500 x 4
#> [1] 2500    4
```

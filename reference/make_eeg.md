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

## What this generator does and does not model

Each channel is generated independently, so the channels are only weakly
correlated: on a 19-channel, 10 s record the median absolute
between-channel correlation is about 0.26. Real EEG is spatially smooth
and neighbouring electrodes correlate strongly, which several methods
rely on. In particular
[`eegBadChannels`](https://x-biosignal.github.io/PhysioEEG/reference/eegBadChannels.md)
judges a channel partly by its correlation with the others against
`corr_threshold` (default 0.4), so on this generator's output it flags
every channel as bad – and anything downstream that needs good channels,
such as
[`eegInterpolate`](https://x-biosignal.github.io/PhysioEEG/reference/eegInterpolate.md),
then has none to work with. That is the generator being independent, not
the detector being wrong. Use it for per-channel operations – filtering,
spectra, epoching, per-channel measures – and use real data, or a
generator with a shared spatial component, when exercising anything that
assumes channels covary.

pe \<- make_eeg(n_time = 2500, n_channels = 4, sr = 250)
dim(SummarizedExperiment::assay(pe, "raw")) \# 2500 x 4

## See also

[`make_eeg_erp()`](https://x-biosignal.github.io/PhysioEEG/reference/make_eeg_erp.md),
[`make_eeg_sleep()`](https://x-biosignal.github.io/PhysioEEG/reference/make_eeg_sleep.md),
[`make_eeg_bci()`](https://x-biosignal.github.io/PhysioEEG/reference/make_eeg_bci.md),
[`make_eeg_spikes()`](https://x-biosignal.github.io/PhysioEEG/reference/make_eeg_spikes.md),
[`eegFilter()`](https://x-biosignal.github.io/PhysioEEG/reference/eegFilter.md),
[`eegCoherence()`](https://x-biosignal.github.io/PhysioEEG/reference/eegCoherence.md)

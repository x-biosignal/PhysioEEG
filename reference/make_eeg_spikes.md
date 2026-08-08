# Create EEG with Embedded Epileptic Spikes

Generates multi-channel EEG with sharp transient spikes at known
locations. Spikes are ~50-100ms duration with high amplitude and sharp
morphology, inserted at random time points across a subset of channels.

## Usage

``` r
make_eeg_spikes(n_time = 30000, n_channels = 19, sr = 500, n_spikes = 15)
```

## Arguments

- n_time:

  Number of time points (default: 30000 = 60s at 500Hz).

- n_channels:

  Number of channels (default: 19).

- sr:

  Sampling rate in Hz (default: 500).

- n_spikes:

  Number of epileptic spikes to insert (default: 15).

## Value

A PhysioExperiment with simulated clinical EEG in the `"raw"` assay.
Sets `metadata(x)$spike_locations` with a data.frame containing columns:
`spike_id`, `sample` (sample index), `time_sec` (time in seconds), and
`channels` (list column of affected channel indices).

## See also

[`make_eeg()`](https://x-biosignal.github.io/PhysioEEG/reference/make_eeg.md),
[`eegSpikeDetect()`](https://x-biosignal.github.io/PhysioEEG/reference/eegSpikeDetect.md),
[`eegQEEG()`](https://x-biosignal.github.io/PhysioEEG/reference/eegQEEG.md),
[`eegSlowing()`](https://x-biosignal.github.io/PhysioEEG/reference/eegSlowing.md)

## Examples

``` r
pe <- make_eeg_spikes(n_time = 5000, n_channels = 4, sr = 250, n_spikes = 3)
dim(SummarizedExperiment::assay(pe, "raw"))  # 5000 x 4
#> [1] 5000    4
```

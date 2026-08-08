# Extract Microstate Sequence as Character Labels

Converts integer microstate labels to character labels ("A", "B", "C",
...).

## Usage

``` r
eegMicrostateSequence(x)
```

## Arguments

- x:

  A PhysioExperiment object with microstate labels in
  `metadata(x)$microstates$labels`.

## Value

Character vector of length n_time with labels "A", "B", "C", etc. Each
element corresponds to the microstate class assigned to that time point.

## References

Michel, C. M., & Koenig, T. (2018). EEG microstates as a tool for
studying the temporal dynamics of whole-brain neuronal networks.
NeuroImage, 180, 577-593.

## See also

[`eegMicrostates()`](https://x-biosignal.github.io/PhysioEEG/reference/eegMicrostates.md),
[`eegMicrostateStats()`](https://x-biosignal.github.io/PhysioEEG/reference/eegMicrostateStats.md),
[`eegMicrostateBackfit()`](https://x-biosignal.github.io/PhysioEEG/reference/eegMicrostateBackfit.md)

## Examples

``` r
if (FALSE) { # \dontrun{
pe <- make_eeg(n_time = 2000, n_channels = 19, sr = 500)
pe <- eegMicrostates(pe, n_states = 4)
seq_labels <- eegMicrostateSequence(pe)
table(seq_labels)
} # }
```

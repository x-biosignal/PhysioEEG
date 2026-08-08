# Compute Microstate Statistics

Calculates temporal statistics for each microstate class from a
segmented EEG recording: mean duration, occurrence rate, and time
coverage. Also computes the transition probability matrix between
states.

## Usage

``` r
eegMicrostateStats(x)
```

## Arguments

- x:

  A PhysioExperiment object with microstate labels in
  `metadata(x)$microstates` (from
  [`eegMicrostates`](https://x-biosignal.github.io/PhysioEEG/reference/eegMicrostates.md)).

## Value

A data.frame with columns:

- state:

  Integer microstate class (1 to n_states).

- duration_ms:

  Mean duration of consecutive runs in milliseconds.

- occurrence_per_sec:

  Number of state occurrences (runs) per second.

- coverage_pct:

  Percentage of total time spent in this state.

The transition probability matrix (n_states x n_states) is stored as an
attribute `"transition_matrix"`.

## References

Michel, C. M., & Koenig, T. (2018). EEG microstates as a tool for
studying the temporal dynamics of whole-brain neuronal networks.
NeuroImage, 180, 577-593.

## See also

[`eegMicrostates()`](https://x-biosignal.github.io/PhysioEEG/reference/eegMicrostates.md),
[`eegMicrostateBackfit()`](https://x-biosignal.github.io/PhysioEEG/reference/eegMicrostateBackfit.md),
[`eegMicrostateSequence()`](https://x-biosignal.github.io/PhysioEEG/reference/eegMicrostateSequence.md)

## Examples

``` r
if (FALSE) { # \dontrun{
pe <- make_eeg(n_time = 5000, n_channels = 19, sr = 500)
pe <- eegMicrostates(pe, n_states = 4, method = "kmeans")
stats <- eegMicrostateStats(pe)
print(stats)
attr(stats, "transition_matrix")
} # }
```

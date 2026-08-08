# Backfit Microstate Template Maps to EEG Data

Assigns each time point to the microstate template map with the highest
absolute spatial correlation. This allows applying microstate maps
derived from one dataset to another dataset.

## Usage

``` r
eegMicrostateBackfit(x, maps, assay_name = NULL)
```

## Arguments

- x:

  A PhysioExperiment object with EEG data.

- maps:

  Numeric matrix of microstate template maps (n_channels x n_states), as
  returned in `metadata(x)$microstates$maps`.

- assay_name:

  Name of the input assay. If `NULL`, the default assay is used.

## Value

Modified PhysioExperiment with updated microstate labels in
`metadata(x)$microstates$labels`. Also stores the template maps in
`metadata(x)$microstates$maps` and the number of states in
`metadata(x)$microstates$n_states`.

## References

Michel, C. M., & Koenig, T. (2018). EEG microstates as a tool for
studying the temporal dynamics of whole-brain neuronal networks.
NeuroImage, 180, 577-593.

## See also

[`eegMicrostates()`](https://x-biosignal.github.io/PhysioEEG/reference/eegMicrostates.md),
[`eegMicrostateStats()`](https://x-biosignal.github.io/PhysioEEG/reference/eegMicrostateStats.md),
[`eegMicrostateSequence()`](https://x-biosignal.github.io/PhysioEEG/reference/eegMicrostateSequence.md)

## Examples

``` r
if (FALSE) { # \dontrun{
pe1 <- make_eeg(n_time = 5000, n_channels = 19, sr = 500)
pe1 <- eegMicrostates(pe1, n_states = 4)
maps <- metadata(pe1)$microstates$maps

pe2 <- make_eeg(n_time = 3000, n_channels = 19, sr = 500)
pe2 <- eegMicrostateBackfit(pe2, maps)
} # }
```

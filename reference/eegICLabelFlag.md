# Flag Artifact ICA Components from ICLabel Probabilities

Convenience wrapper over
[`eegICLabel()`](https://x-biosignal.github.io/PhysioEEG/reference/eegICLabel.md)
that returns the indices of the components most likely to be artifacts,
ready to pass to
[`eegICAremove()`](https://x-biosignal.github.io/PhysioEEG/reference/eegICAremove.md).
A component is flagged when its probability in any of the requested
artifact `classes` exceeds `prob_threshold`.

## Usage

``` r
eegICLabelFlag(
  x,
  prob_threshold = 0.5,
  classes = c("muscle", "eye", "heart", "line_noise", "channel_noise"),
  ica_assay = "ica_components",
  line_freq = 50
)
```

## Arguments

- x:

  A PhysioExperiment object with ICA results (from
  [`eegICA()`](https://x-biosignal.github.io/PhysioEEG/reference/eegICA.md)).

- prob_threshold:

  Probability above which a component is flagged (default: 0.5).

- classes:

  Character vector of artifact classes to flag (default: all non-brain,
  non-other classes).

- ica_assay:

  Metadata name holding the component activations (default:
  `"ica_components"`).

- line_freq:

  Mains line frequency in Hz (default: 50).

## Value

An integer vector of flagged component indices (possibly empty), with an
attribute `"labels"` giving the corresponding class of each.

## See also

[`eegICLabel()`](https://x-biosignal.github.io/PhysioEEG/reference/eegICLabel.md),
[`eegICAremove()`](https://x-biosignal.github.io/PhysioEEG/reference/eegICAremove.md),
[`eegICAdetect()`](https://x-biosignal.github.io/PhysioEEG/reference/eegICAdetect.md)

## Examples

``` r
if (FALSE) { # \dontrun{
pe <- make_eeg(n_time = 5000, sr = 500)
pe <- eegICA(pe, n_components = 10, method = "fastica")
bad <- eegICLabelFlag(pe, prob_threshold = 0.5)
pe <- eegICAremove(pe, components = bad)
} # }
```

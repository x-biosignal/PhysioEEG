# Motor-Imagery Lateralization Index (hemispheric ERD/ERS)

Quantifies the hemispheric asymmetry of sensorimotor-rhythm activity
during motor imagery, a rehabilitation-relevant index of affected-
versus unaffected-hemisphere engagement. For each trial a laterality
index \\LI = (R - L) / (\|R\| + \|L\|)\\ is computed from a right- and a
left-hemisphere channel (or ROI), bounded from -1 to 1. With
`method = "erd"` the activation is the event-related desynchronization
magnitude from
[`eegMotorImagery()`](https://x-biosignal.github.io/PhysioEEG/reference/eegMotorImagery.md)
(a positive LI means the right hemisphere desynchronizes more, i.e. is
more active); with `method = "power"` it is the raw band power (a
positive LI means the right hemisphere has more power). A symmetric
input gives \\LI \approx 0\\.

## Usage

``` r
eegLateralization(
  x,
  left_ch = "C3",
  right_ch = "C4",
  band = c(8, 13),
  method = c("erd", "power"),
  baseline_fraction = 0.25,
  conf_level = 0.95,
  reliability = list(icc = NA_real_, sem = NA_real_),
  assay_name = NULL
)
```

## Arguments

- x:

  A PhysioExperiment object with epoched (3D) EEG data (time x channels
  x trials).

- left_ch:

  Character label(s) of the left-hemisphere channel or ROI (default:
  `"C3"`).

- right_ch:

  Character label(s) of the right-hemisphere channel or ROI (default:
  `"C4"`).

- band:

  Numeric length-2 frequency band in Hz (default: `c(8, 13)`, the mu
  rhythm).

- method:

  Activation measure: `"erd"` (desynchronization magnitude, via
  [`eegMotorImagery()`](https://x-biosignal.github.io/PhysioEEG/reference/eegMotorImagery.md))
  or `"power"` (raw band power).

- baseline_fraction:

  Fraction of each trial used as the ERD baseline (default: 0.25); only
  used by `method = "erd"`.

- conf_level:

  Confidence level for the summary interval (default: 0.95).

- reliability:

  Named list of reliability indices attached to the summary biomarker
  (default: `icc` and `sem` placeholders).

- assay_name:

  Input assay name (default: the default assay).

## Value

A list with `per_trial` (an `n_trials` x 1 matrix of laterality indices)
and `summary` (a
[`PhysioBiomarker`](https://x-biosignal.github.io/PhysioCore//reference/PhysioBiomarker.html)
holding the mean LI with its confidence interval, reliability, and
provenance).

## References

Pfurtscheller, G., & Lopes da Silva, F. H. (1999). Event-related EEG/MEG
synchronization and desynchronization: basic principles. Clinical
Neurophysiology, 110(11), 1842-1857.

Bai, O., et al. (2005). Asymmetric spatiotemporal patterns of
event-related desynchronization preceding voluntary sequential finger
movements. Clinical Neurophysiology, 116(5), 1213-1221.

## See also

[`eegMotorImagery()`](https://x-biosignal.github.io/PhysioEEG/reference/eegMotorImagery.md),
[`eegCSP()`](https://x-biosignal.github.io/PhysioEEG/reference/eegCSP.md),
[`PhysioCore::physioBiomarker()`](https://x-biosignal.github.io/PhysioCore//reference/physioBiomarker-constructor.html)

## Examples

``` r
if (FALSE) { # \dontrun{
pe <- make_eeg_bci(n_trials = 20, n_channels = 8, sr = 256)
li <- eegLateralization(pe, method = "power")
li$summary
} # }
```

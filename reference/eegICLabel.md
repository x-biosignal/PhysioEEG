# ICLabel-style ICA Component Classification

Classifies independent components (from
[`eegICA()`](https://x-biosignal.github.io/PhysioEEG/reference/eegICA.md))
into seven classes - `brain`, `muscle`, `eye`, `heart`, `line_noise`,
`channel_noise`, and `other` - and returns a calibrated probability for
each class per component together with the argmax label. The approach
follows the ICLabel framework: interpretable spatial, spectral, and
temporal features are extracted from each component and mapped to class
probabilities by a lightweight multinomial-logistic (softmax) head.

## Usage

``` r
eegICLabel(x, ica_assay = "ica_components", assay_name = NULL, line_freq = 50)
```

## Arguments

- x:

  A PhysioExperiment object with ICA results (from
  [`eegICA()`](https://x-biosignal.github.io/PhysioEEG/reference/eegICA.md)).

- ica_assay:

  Metadata name holding the component activations (default:
  `"ica_components"`).

- assay_name:

  Input assay used only to resolve channel labels (default: first
  assay).

- line_freq:

  Mains line frequency in Hz used for the line-noise features (default:
  50).

## Value

A data.frame with one row per component: `component` (integer index),
one numeric column per class (`brain`, `muscle`, `eye`, `heart`,
`line_noise`, `channel_noise`, `other`) holding probabilities that sum
to 1 across the classes, and `label` (character argmax class).

## Details

Features per component:

- Spatial: frontal energy (eye/blink topography), focality (single
  channel dominance), topography kurtosis.

- Spectral: 1/f slope, high- versus low-frequency band ratio,
  low-frequency fraction, a genuine alpha-peak measure (8-12 Hz power
  above its theta/low-beta neighbours), and a line-noise power ratio and
  fraction.

- Temporal: lag-1 autocorrelation, activation kurtosis (spiky blink or
  ECG signatures), and a roughly 1 Hz periodicity measure for cardiac
  components.

The softmax weights are read from `inst/extdata/iclabel_weights.csv`; if
that file is unavailable an identical built-in weight table is used, so
the classifier always works.

## References

Pion-Tonachini, L., Kreutz-Delgado, K., & Makeig, S. (2019). ICLabel: An
automated electroencephalographic independent component classifier,
dataset, and website. NeuroImage, 198, 181-197.

Winkler, I., Haufe, S., & Tangermann, M. (2011). Automatic
classification of artifactual ICA-components for artifact removal in EEG
signals. Behavioral and Brain Functions, 7, 30.

## See also

[`eegICA()`](https://x-biosignal.github.io/PhysioEEG/reference/eegICA.md),
[`eegICLabelFlag()`](https://x-biosignal.github.io/PhysioEEG/reference/eegICLabelFlag.md),
[`eegICAdetect()`](https://x-biosignal.github.io/PhysioEEG/reference/eegICAdetect.md),
[`eegICAremove()`](https://x-biosignal.github.io/PhysioEEG/reference/eegICAremove.md)

## Examples

``` r
if (FALSE) { # \dontrun{
pe <- make_eeg(n_time = 5000, sr = 500)
pe <- eegICA(pe, n_components = 10, method = "fastica")
probs <- eegICLabel(pe)
} # }
```

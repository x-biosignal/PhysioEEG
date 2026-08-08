# Frontal Alpha Asymmetry

Computes frontal alpha asymmetry indices from paired electrode sites.
For each pair (right, left), band power is computed via FFT and the
asymmetry index is calculated as `log(power_right) - log(power_left)`.
Positive values indicate greater right-hemisphere alpha (typically
associated with greater left-hemisphere cortical activity).

## Usage

``` r
eegAsymmetry(x, pairs = NULL, band = c(8, 13), assay_name = NULL)
```

## Arguments

- x:

  A PhysioExperiment object with EEG data.

- pairs:

  A list of character vectors of length 2, each specifying a (right,
  left) electrode pair. Defaults to
  `list(c("F4", "F3"), c("F8", "F7"))`.

- band:

  Numeric vector of length 2 specifying the frequency band in Hz
  (default: `c(8, 13)` for alpha).

- assay_name:

  Name of the input assay. If `NULL`, the default assay is used.

## Value

A data.frame with columns:

- pair:

  Character label for the electrode pair.

- left_channel:

  Character name of the left electrode.

- right_channel:

  Character name of the right electrode.

- left_power:

  Numeric band power for the left electrode.

- right_power:

  Numeric band power for the right electrode.

- asymmetry_index:

  Numeric asymmetry: log(right) - log(left).

## References

Nuwer, M. R., et al. (1999). IFCN standards for digital recording of
clinical EEG. Electroencephalography and Clinical Neurophysiology,
106(3), 259-261.

## See also

[`eegSpikeDetect()`](https://x-biosignal.github.io/PhysioEEG/reference/eegSpikeDetect.md),
[`eegQEEG()`](https://x-biosignal.github.io/PhysioEEG/reference/eegQEEG.md),
[`eegSlowing()`](https://x-biosignal.github.io/PhysioEEG/reference/eegSlowing.md)

## Examples

``` r
if (FALSE) { # \dontrun{
pe <- make_eeg(n_time = 5000, n_channels = 19, sr = 500)
asym <- eegAsymmetry(pe)
print(asym)
} # }
```

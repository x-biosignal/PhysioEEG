# Conditional (multivariate) Granger causality

Computes conditional Granger causality from `source` to `target` given a
set of `conditioning` channels, using the log ratio of the target's
residual variance in a VAR fitted without the source
(`target + conditioning`) versus with it
(`target + source + conditioning`); both VARs share the same order
(Geweke; Barnett & Seth 2014). Unlike the pairwise
[`eegGrangerCausality()`](https://x-biosignal.github.io/PhysioEEG/reference/eegGrangerCausality.md),
conditioning removes indirect influences, so for a chain the conditional
GC along a purely indirect path is near zero.

## Usage

``` r
eegConditionalGC(
  x,
  target,
  source,
  conditioning = NULL,
  order = NULL,
  method = "ols",
  assay_name = NULL
)
```

## Arguments

- x:

  A PhysioExperiment object with 2D EEG data (time x channels).

- target:

  Target channel (label or index) - the effect.

- source:

  Source channel (label or index) - the putative cause.

- conditioning:

  Channels to condition on (labels or indices). Defaults to all
  remaining channels.

- order:

  MVAR model order, or `NULL` to select automatically on the full
  channel set.

- method:

  MVAR estimator passed to
  [`PhysioCore::mvarFit()`](https://x-biosignal.github.io/PhysioCore//reference/mvarFit.html)
  (default: `"ols"`).

- assay_name:

  Input assay name (default: the default assay).

## Value

A list with the conditional GC `value` (in nats, non-negative in
theory), the `target`, `source`, and `conditioning` labels, and the VAR
`order` used.

## References

Barnett, L., & Seth, A. K. (2014). The MVGC multivariate Granger
causality toolbox. Journal of Neuroscience Methods, 223, 50-68.

## See also

[`eegGrangerCausality()`](https://x-biosignal.github.io/PhysioEEG/reference/eegGrangerCausality.md),
[`eegDTF()`](https://x-biosignal.github.io/PhysioEEG/reference/eegDTF.md),
[`eegPDC()`](https://x-biosignal.github.io/PhysioEEG/reference/eegPDC.md)

## Examples

``` r
if (FALSE) { # \dontrun{
pe <- make_eeg(n_time = 4000, n_channels = 3, sr = 250)
eegConditionalGC(pe, target = 3, source = 1, conditioning = 2)
} # }
```

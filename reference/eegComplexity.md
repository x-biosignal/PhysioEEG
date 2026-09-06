# EEG complexity and nonlinear-dynamics measures

Computes per-channel entropy, fractal-dimension and
detrended-fluctuation complexity measures — the nonlinear family that
complements the spectral, connectivity and source tools. For epoched
(3-D) data each measure is computed per epoch and averaged across
epochs.

## Usage

``` r
eegComplexity(
  pe,
  measures = c("sample_entropy", "permutation_entropy", "lempel_ziv", "higuchi_fd",
    "dfa", "hjorth_mobility", "spectral_entropy"),
  assay_name = NULL,
  m = 2L,
  r = 0.2,
  tau = 1L,
  perm_order = 3L,
  mse_scales = 1:8,
  max_samples = 4000L
)
```

## Arguments

- pe:

  A `PhysioExperiment`.

- measures:

  Character vector of measures to compute (see Details).

- assay_name:

  Assay to use (default: the object's default assay).

- m:

  Embedding dimension for the entropy measures (default 2; permutation
  entropy uses `perm_order`).

- r:

  Tolerance as a fraction of each channel's SD for sample/approximate/
  multiscale entropy (default 0.2).

- tau:

  Time delay for permutation entropy (default 1).

- perm_order:

  Embedding order for permutation entropy (default 3).

- mse_scales:

  Scales for multiscale entropy (default `1:8`).

- max_samples:

  Cap on samples per channel/epoch for the O(N^2) entropy measures;
  longer series are truncated with a warning (default 4000; `NULL`
  disables the cap).

## Value

A data frame with one row per channel and one column per requested
measure (plus a `channel` column). The multiscale-entropy per-scale
curves are attached as `attr(, "mse_scales")` / `attr(, "mse_curve")`.

## Details

Available `measures`: `"sample_entropy"` (Richman–Moorman),
`"approximate_entropy"` (Pincus), `"permutation_entropy"` (Bandt–Pompe),
`"multiscale_entropy"` (Costa; returns the scale-averaged SampEn as a
summary plus the per-scale curve in the attribute), `"lempel_ziv"`
(Kaspar–Schuster), `"higuchi_fd"`, `"katz_fd"`, `"dfa"` (scaling
exponent alpha), `"hurst"` (R/S), `"hjorth_mobility"`,
`"hjorth_complexity"`, `"spectral_entropy"`.

## References

Richman & Moorman 2000; Bandt & Pompe 2002; Costa 2002; Higuchi 1988;
Katz 1988; Peng 1994; Hjorth 1970.

## See also

[`eegQEEG()`](https://x-biosignal.github.io/PhysioEEG/reference/eegQEEG.md),
[`eegAperiodic()`](https://x-biosignal.github.io/PhysioEEG/reference/eegAperiodic.md)

## Examples

``` r
pe <- make_eeg(n_time = 512, n_channels = 4, sr = 128)
cx <- eegComplexity(pe, measures = c("permutation_entropy", "hjorth_mobility",
                                     "spectral_entropy"))
cx
#>   channel permutation_entropy hjorth_mobility spectral_entropy
#> 1     Fp1           0.7293826       0.3545342        0.2493087
#> 2     Fp2           0.8205204       0.1419271        0.3126944
#> 3      F7           0.9005079       0.4122412        0.5447787
#> 4      F3           0.8649210       0.4188192        0.3752383
```

# Statistical Testing for ERP Differences

Performs permutation testing or cluster-based permutation testing to
compare ERP waveforms between conditions.

## Usage

``` r
eegERPtest(
  x,
  y,
  method = c("permutation", "cluster"),
  n_perm = 1000,
  alpha = 0.05,
  cluster_alpha = 0.05,
  assay_name = NULL
)
```

## Arguments

- x:

  A PhysioExperiment with epoched data for condition 1.

- y:

  A PhysioExperiment with epoched data for condition 2.

- method:

  Test method: `"permutation"` (pointwise permutation test) or
  `"cluster"` (cluster-based permutation test).

- n_perm:

  Number of permutations (default: 1000).

- alpha:

  Significance level (default: 0.05).

- cluster_alpha:

  Cluster-forming threshold for individual t-tests (default: 0.05). Only
  used for `"cluster"` method.

- assay_name:

  Input assay name (default: first assay).

## Value

A data.frame with columns: `time_sample` (integer), `t_statistic`
(numeric observed t-value), `p_value` (numeric permutation-based
p-value), and `significant` (logical). For the `"cluster"` method, also
includes `cluster_id` (integer cluster assignment) and `cluster_p`
(numeric cluster-level corrected p-value).

## References

Luck, S. J. (2014). An Introduction to the Event-Related Potential
Technique (2nd ed.). MIT Press.

Maris, E., & Oostenveld, R. (2007). Nonparametric statistical testing of
EEG- and MEG-data. Journal of Neuroscience Methods, 164(1), 177-190.

## See also

[`eegERPdetect()`](https://x-biosignal.github.io/PhysioEEG/reference/eegERPdetect.md),
[`eegERPmeasure()`](https://x-biosignal.github.io/PhysioEEG/reference/eegERPmeasure.md),
[`eegERPdifference()`](https://x-biosignal.github.io/PhysioEEG/reference/eegERPdifference.md),
[`eegERPgrandAverage()`](https://x-biosignal.github.io/PhysioEEG/reference/eegERPgrandAverage.md)

## Examples

``` r
if (FALSE) { # \dontrun{
pe1 <- make_eeg_erp(n_epochs = 20, sr = 250)
pe2 <- make_eeg_erp(n_epochs = 20, sr = 250)
result <- eegERPtest(pe1, pe2, method = "permutation", n_perm = 500)
} # }
```

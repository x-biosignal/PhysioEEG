# Remove ICA Components from EEG

Reconstructs EEG data with specified independent components removed. The
removed components are zeroed out in the mixing matrix before
back-projecting to channel space.

## Usage

``` r
eegICAremove(
  x,
  components,
  ica_assay = "ica_components",
  output_assay = "ica_cleaned"
)
```

## Arguments

- x:

  A PhysioExperiment object with ICA results (from `eegICA`).

- components:

  Integer vector of component indices to remove.

- ica_assay:

  Assay name containing ICA component activations (default: `"ica"`).

- output_assay:

  Output assay name (default: `"ica_cleaned"`).

## Value

Modified PhysioExperiment with cleaned data in `output_assay`. The
cleaned assay contains reconstructed channel data with specified
components removed. Dimensions match the original data.

## References

Hyvarinen, A., & Oja, E. (2000). Independent component analysis:
algorithms and applications. Neural Networks, 13(4-5), 411-430.

Bell, A. J., & Sejnowski, T. J. (1995). An information-maximization
approach to blind separation and blind deconvolution. Neural
Computation, 7(6), 1129-1159.

## See also

[`eegICA()`](https://x-biosignal.github.io/PhysioEEG/reference/eegICA.md),
[`eegICAdetect()`](https://x-biosignal.github.io/PhysioEEG/reference/eegICAdetect.md),
[`eegICAmix()`](https://x-biosignal.github.io/PhysioEEG/reference/eegICAmix.md)

## Examples

``` r
if (FALSE) { # \dontrun{
pe <- make_eeg(n_time = 5000, sr = 500)
pe <- eegICA(pe, n_components = 4, method = "fastica")
pe <- eegICAremove(pe, components = c(1, 2))
} # }
```

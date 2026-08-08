# EEG Independent Component Analysis (ICA)

Decomposes multi-channel EEG into independent components using ICA.
Supports FastICA, Infomax, and JADE algorithms. Results are stored in
the output assay (component activations) and in `metadata(x)$ica`
(mixing and unmixing matrices).

## Usage

``` r
eegICA(
  x,
  n_components = NULL,
  method = c("fastica", "infomax", "jade"),
  max_iter = 200L,
  tol = 1e-06,
  assay_name = NULL,
  output_assay = "ica_components"
)
```

## Arguments

- x:

  A PhysioExperiment object with EEG data.

- n_components:

  Number of independent components to extract. Defaults to the number of
  channels.

- method:

  ICA algorithm: `"fastica"`, `"infomax"`, or `"jade"`.

- max_iter:

  Maximum number of iterations (default: 200).

- tol:

  Convergence tolerance (default: 1e-6).

- assay_name:

  Input assay name (default: first assay).

- output_assay:

  Output assay name (default: `"ica"`).

## Value

Modified PhysioExperiment with component activations in `output_assay`
and ICA metadata in `metadata(x)$ica`. The ICA metadata list contains:
`mixing` (mixing matrix A), `unmixing` (unmixing matrix), `mean`
(channel means), `whiten` (whitening matrix), and `method` (algorithm
used). The output assay has dimensions n_time x n_components.

## References

Hyvarinen, A., & Oja, E. (2000). Independent component analysis:
algorithms and applications. Neural Networks, 13(4-5), 411-430.

Bell, A. J., & Sejnowski, T. J. (1995). An information-maximization
approach to blind separation and blind deconvolution. Neural
Computation, 7(6), 1129-1159.

## See also

[`eegICAremove()`](https://x-biosignal.github.io/PhysioEEG/reference/eegICAremove.md),
[`eegICAdetect()`](https://x-biosignal.github.io/PhysioEEG/reference/eegICAdetect.md),
[`eegICAmix()`](https://x-biosignal.github.io/PhysioEEG/reference/eegICAmix.md),
[`eegFilter()`](https://x-biosignal.github.io/PhysioEEG/reference/eegFilter.md),
[`eegPreprocess()`](https://x-biosignal.github.io/PhysioEEG/reference/eegPreprocess.md)

## Examples

``` r
if (FALSE) { # \dontrun{
pe <- make_eeg(n_time = 5000, sr = 500)
result <- eegICA(pe, n_components = 4, method = "fastica")
} # }
```

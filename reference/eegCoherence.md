# EEG Coherence Analysis

Computes magnitude-squared coherence (MSC) or imaginary part of
coherency between all EEG channel pairs using Welch's method. The signal
is segmented into overlapping windows, each windowed with a Hanning
function, and the cross-spectral density is averaged across windows.

## Usage

``` r
eegCoherence(
  x,
  method = c("coherence", "imaginary"),
  window_sec = 2,
  overlap = 0.5,
  band = c(8, 13),
  assay_name = NULL
)
```

## Arguments

- x:

  A PhysioExperiment object with EEG data (2D: time x channels).

- method:

  Coherence method: `"coherence"` for magnitude-squared coherence or
  `"imaginary"` for imaginary part of coherency (default:
  `"coherence"`).

- window_sec:

  Window length in seconds for Welch's method (default: 2).

- overlap:

  Overlap fraction between adjacent windows, from 0 to 1 exclusive
  (default: 0.5).

- band:

  Numeric vector of length 2 specifying the frequency band in Hz over
  which to average coherence (default: `c(8, 13)` for alpha band).

- assay_name:

  Name of the input assay. If `NULL`, the default assay is used.

## Value

The input PhysioExperiment with connectivity results stored in
`metadata(x)$connectivity`, a list containing:

- matrix:

  Numeric n_channels x n_channels matrix of band-averaged coherence
  values.

- method:

  Character string indicating the method used.

- band:

  Numeric vector of the frequency band used.

- freqs:

  Numeric vector of frequency bins.

- spectra:

  3D array (n_freqs x n_channels x n_channels) of frequency-resolved
  coherence values.

## Details

Magnitude-squared coherence is defined as: \$\$MSC(f) = \|S\_{xy}(f)\|^2
/ (S\_{xx}(f) \cdot S\_{yy}(f))\$\$

Imaginary coherence uses the imaginary part of coherency to reduce
volume conduction artifacts: \$\$ICoh(f) = \|Im(S\_{xy}(f) /
\sqrt{S\_{xx}(f) \cdot S\_{yy}(f)})\|\$\$

## References

Lachaux, J. P., et al. (1999). Measuring phase synchrony in brain
signals. Human Brain Mapping, 8(4), 194-208.

## See also

[`eegPLV()`](https://x-biosignal.github.io/PhysioEEG/reference/eegPLV.md),
[`eegWPLI()`](https://x-biosignal.github.io/PhysioEEG/reference/eegWPLI.md),
[`eegConnectivityMatrix()`](https://x-biosignal.github.io/PhysioEEG/reference/eegConnectivityMatrix.md),
[`eegPlotConnectivity()`](https://x-biosignal.github.io/PhysioEEG/reference/eegPlotConnectivity.md)

## Examples

``` r
if (FALSE) { # \dontrun{
pe <- make_eeg(n_time = 5000, n_channels = 4, sr = 500)
pe <- eegCoherence(pe, method = "coherence", band = c(8, 13))
coh_matrix <- metadata(pe)$connectivity$matrix
} # }
```

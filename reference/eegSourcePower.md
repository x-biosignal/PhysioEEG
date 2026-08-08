# Compute Source-Space Band Power

Computes spectral band power for each source from source-estimated EEG
data. Requires prior source estimation via
[`eegSourceEstimate`](https://x-biosignal.github.io/PhysioEEG/reference/eegSourceEstimate.md).

## Usage

``` r
eegSourcePower(x, bands = NULL, source_assay = "source")
```

## Arguments

- x:

  A PhysioExperiment object with source estimates.

- bands:

  Named list of frequency bands, each a numeric vector of length 2
  (lower, upper Hz). If `NULL`, uses standard EEG bands: delta (1-4),
  theta (4-8), alpha (8-13), beta (13-30), gamma (30-50).

- source_assay:

  Name of the assay containing source data (default: `"source"`).

## Value

A data.frame with columns:

- source_id:

  Integer source index.

- band:

  Character name of the frequency band.

- power:

  Numeric spectral power in the band.

## References

Pascual-Marqui, R. D. (2002). Standardized low-resolution brain
electromagnetic tomography (sLORETA). Methods and Findings in
Experimental and Clinical Pharmacology, 24(Suppl D), 5-12.

## See also

[`eegSourceEstimate()`](https://x-biosignal.github.io/PhysioEEG/reference/eegSourceEstimate.md),
[`eegForwardModel()`](https://x-biosignal.github.io/PhysioEEG/reference/eegForwardModel.md),
[`eegBeamformer()`](https://x-biosignal.github.io/PhysioEEG/reference/eegBeamformer.md)

## Examples

``` r
if (FALSE) { # \dontrun{
pe <- make_eeg(n_time = 1000, n_channels = 19, sr = 250)
fm <- eegForwardModel(pe, method = "spherical", n_sources = 20)
pe <- eegSourceEstimate(pe, fm, method = "mne")
bp <- eegSourcePower(pe)
head(bp)
} # }
```

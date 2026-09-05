# Changelog

## PhysioEEG 0.7.6

### New features

- [`petrosianFD()`](https://x-biosignal.github.io/PhysioEEG/reference/petrosianFD.md)
  — Petrosian fractal dimension (Petrosian 1995): a fast O(N)
  waveform-complexity index from the number of sign changes in the
  signal’s derivative (local-extrema density). Smooth oscillation → FD
  near 1; complex/noisy → higher. Complements
  [`svdEntropy()`](https://x-biosignal.github.io/PhysioEEG/reference/svdEntropy.md)
  and
  [`eegComplexity()`](https://x-biosignal.github.io/PhysioEEG/reference/eegComplexity.md);
  reproduces `antropy.petrosian_fd` bit-for-bit (and, unlike the
  Higuchi/Katz FDs, has a single unambiguous closed form).

## PhysioEEG 0.7.5

### New features

- [`svdEntropy()`](https://x-biosignal.github.io/PhysioEEG/reference/svdEntropy.md)
  — Singular Value Decomposition entropy (Roberts et al. 1999):
  time-delay embed a signal, take the normalized singular values of the
  embedding matrix, and return their Shannon entropy — a measure of the
  signal’s dimensionality (low = structured/oscillatory, high =
  complex/noise-like). A linear-algebraic complement to
  [`eegComplexity()`](https://x-biosignal.github.io/PhysioEEG/reference/eegComplexity.md)’s
  entropy/fractal measures; reproduces `antropy.svd_entropy` bit-for-bit
  (both via LAPACK SVD).

## PhysioEEG 0.7.4

- Test-suite performance: `test-eeg-timefreq.R` runs **12× faster**
  (368s → 31s) by right-sizing the synthetic test data (`n_time` 5000 →
  2000; ERSP/ITC `n_epochs` 40/20 → 12/10). Frequency grids and every
  assertion are unchanged, so coverage is preserved (74 tests, all
  passing). The full package test suite now completes in **~227s**
  (previously exceeded a 500s CI budget), all 1282 tests passing.

## PhysioEEG 0.7.3

- `eegICLabel(backend = "iclabel")` runs the **genuine trained ICLabel
  CNN** (Pion-Tonachini et al. 2019), replacing the need to trust the
  pure-R heuristic when the real model is wanted. It delegates to the
  validated `mne-icalabel` implementation through `reticulate` (carrying
  the PhysioEEG ICA mixing/unmixing matrices into an MNE `ICA` object),
  so it uses mne-icalabel’s exact EEGLAB-style feature extraction and
  the trained weights rather than re-porting the network natively (which
  would risk a plausible-but-wrong result). The default
  `backend = "heuristic"` is unchanged and needs no Python; the trained
  backend is optional (needs `reticulate` + a Python env with `mne` and
  `mne-icalabel`, and standard 10-20/10-10 channel labels for the scalp
  map). Verified to return valid 7-class probability distributions from
  the real network.

## PhysioEEG 0.7.2

- `eegSourceEstimate(method = "dspm")` adds **dSPM** (dynamic
  statistical parametric mapping; Dale et al. 2000): the minimum-norm
  estimate noise-normalized by each source’s projected noise sensitivity
  (`sqrt(diag(K K'))` for an identity noise covariance). Verified to be
  MNE scaled per source by a time-invariant factor.
- [`eegDipoleFit()`](https://x-biosignal.github.io/PhysioEEG/reference/eegDipoleFit.md)
  adds **equivalent-current-dipole (ECD) fitting** — the focal
  complement to the distributed inverses. It localizes the single dipole
  that best explains a scalp topography (grid search inside the head
  sphere refined with Nelder-Mead; the dipole moment is the closed-form
  least-squares fit at each position), using the same forward physics as
  [`eegForwardModel()`](https://x-biosignal.github.io/PhysioEEG/reference/eegForwardModel.md).
  Verified to recover a planted dipole’s position (\< 0.02 head radii),
  orientation, and \> 99.9% goodness of fit. Fits the peak
  global-field-power sample by default, or a chosen sample / window.

## PhysioEEG 0.7.1

- [`eegSurfaceLaplacian()`](https://x-biosignal.github.io/PhysioEEG/reference/eegSurfaceLaplacian.md)
  adds the **surface Laplacian / current source density** (CSD) of
  Perrin et al. (1989) — the spherical-spline scalp Laplacian (Kayser &
  Tenke 2006). CSD is reference-free and a spatial high-pass, deblurring
  volume conduction so each channel reflects the radial current beneath
  it. Verified on its defining properties: invariance to a reference
  constant (the spline constant term absorbs it) and spatial-high-pass
  behaviour (focal topographies are emphasized over smooth ones ~8x).
  Reuses the `.legendre_poly` recurrence and standard 10-10/10-20
  position tables; needs electrode positions (`colData`
  `pos_x/pos_y/pos_z`, else matched to a standard montage by label).

## PhysioEEG 0.7.0

Three research-grade method families that were missing from the EEG
stack.

- [`eegComplexity()`](https://x-biosignal.github.io/PhysioEEG/reference/eegComplexity.md)
  adds the **nonlinear / complexity** family (previously absent
  ecosystem-wide for EEG): per-channel sample & approximate entropy
  (Richman–Moorman, Pincus), permutation entropy (Bandt–Pompe),
  multiscale entropy (Costa), Lempel–Ziv complexity, Higuchi & Katz
  fractal dimension, detrended fluctuation analysis (alpha), Hurst
  (R/S), Hjorth parameters, and spectral entropy. The numeric cores are
  validated against known values (white noise Higuchi FD → 2, Brownian →
  1.5; white-noise DFA/Hurst → 0.5). Entropy measures are O(N^2) and cap
  long channels via `max_samples`.
- [`eegAperiodic()`](https://x-biosignal.github.io/PhysioEEG/reference/eegAperiodic.md)
  adds **aperiodic / 1-f spectral parameterization** (specparam / FOOOF;
  Donoghue 2020): per-channel aperiodic exponent / offset / knee plus
  separated oscillatory peaks, ready for a topomap or a qEEG biomarker.
  Delegates the fit to
  [`PhysioAnalysis::specparam()`](https://x-biosignal.r-universe.dev/PhysioAnalysis/reference/specparam.html).
- [`eegPAC()`](https://x-biosignal.github.io/PhysioEEG/reference/eegPAC.md)
  and
  [`eegComodulogram()`](https://x-biosignal.github.io/PhysioEEG/reference/eegComodulogram.md)
  add **phase-amplitude / cross-frequency coupling** (Tort modulation
  index, Canolty mean-vector-length, Ozkurt, PLV), within a channel or
  across two channels, plus the comodulogram grid. Delegate the
  estimators to
  [`PhysioCrossModal::phaseAmplitudeCoupling()`](https://x-biosignal.github.io/PhysioCrossModal/reference/phaseAmplitudeCoupling.html)
  / `comodulogram()`.

These reuse the existing specparam and PAC engines (new `PhysioAnalysis`
Suggests; `PhysioCrossModal` already suggested) rather than
reimplementing them, and add the genuinely new complexity module
natively. No existing behaviour changed.

## PhysioEEG 0.6.5

- [`make_eeg_bci()`](https://x-biosignal.github.io/PhysioEEG/reference/make_eeg_bci.md)
  now pads channel labels with generic names when more channels than the
  eight-entry motor/occipital pool are requested (previously
  `n_channels > 8` produced `NA` labels, which made the per-channel
  artifact tests error with “missing value where TRUE/FALSE needed”).
- Documentation: the `bci-classification`, `connectivity-analysis`,
  `preprocessing-pipeline`, `sleep-analysis`, and
  `time-frequency-analysis` vignettes were modernized to the current
  package API. They had drifted from the code (stale/renamed arguments
  and functions, changed return structures, and 2D-vs-3D data-shape
  assumptions) and no longer ran. Fixes attach `SummarizedExperiment`,
  use the current `make_eeg*()` signatures and the real
  accessor/analysis functions (`eegFilter(method=)`, `eegRereference`,
  `eegMontage`, `eegBadChannels`, `eegInterpolate`, `eegICA`,
  `eegMorletWavelet`, `eegSTFT`, `eegSleepStage`, the
  spindle/K-complex/slow-wave detectors, etc.), and read
  time-frequency/staging output from `metadata()`. All package vignettes
  now run clean under `R CMD check`. No package behaviour changed.

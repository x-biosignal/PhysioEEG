# Changelog

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

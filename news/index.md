# Changelog

## PhysioEEG 0.6.3

- Made the glass-brain outline hash-governance test portable across the
  r-universe build matrix. It previously shelled out to `sha256sum` and
  compared the tool’s whole output line; Windows Rtools’ `sha256sum`
  prints a binary-mode ’ \*’ separator vs coreutils’ ’ ’, so the
  byte-identical hashes mismatched and the check ERRORed on Windows. The
  SHA-256 is now computed in-process with `digest::digest(file=)` and
  compared hash-only. Added `digest` to Suggests and a `.gitattributes`
  pinning the sealed artifacts so no checkout can EOL-rewrite them.
  Governance intent (offline, closed, tamper-evident) is unchanged.

## PhysioEEG 0.6.2

### Bug fixes

- Added the missing `tests/testthat.R` harness. The package’s 26
  testthat files (including the MNE parity, connectivity ground-truth,
  and golden suites) existed but were never executed by `R CMD check` or
  coverage tooling, so checks passed with zero tests run. The full suite
  passes (0 failures); no package code changed.

## PhysioEEG 0.6.1

### Validation

- Added ground-truth validation of phase-based connectivity (VAL-09):
  [`eegPLV()`](https://x-biosignal.github.io/PhysioEEG/reference/eegPLV.md)
  returns ~1 for phase-locked oscillators and decreases when signals
  decouple;
  [`eegWPLI()`](https://x-biosignal.github.io/PhysioEEG/reference/eegWPLI.md)
  rejects zero-lag volume conduction (wPLI ~ 0) while detecting genuine
  lagged coupling (wPLI ~ 1). These analytic checks replace raw numeric
  parity with MNE-Python for phase connectivity, since PhysioEEG’s
  continuous-Hilbert estimator and MNE’s epoch-spectral estimator differ
  by construction (structural agreement r ~ 0.95).

## PhysioEEG 0.6.0

### New Features

- [`eegPlotConnectogram()`](https://x-biosignal.github.io/PhysioEEG/reference/eegPlotConnectogram.md)
  adds deterministic circular connectivity views with hemisphere, lobe,
  or stable hierarchical-cluster ordering; contiguous module-run arcs;
  and optional cubic Bezier edge bundling.
- Symmetric matrices emit one canonical unordered edge. Asymmetric
  matrices retain every ordered entry and follow the package convention
  of matrix rows as targets and columns as sources.

### Validation and Interpretation

- Strict `abs(value) > threshold` filtering is a descriptive display
  choice, not a p-value or corrected significance mask. Exact node,
  edge, path, module arc, method, band, symmetry, and diagnostic data
  remain inspectable on the returned plot.
- EEG anatomy is inferred only from recognized 10-20/10-10 labels.
  Unknown device labels remain explicitly `unknown`, and tied cluster
  layouts fall back to stable UTF-8 byte order with a diagnostic.
  Length-prefixed edge identifiers remain unique even when channel
  labels contain display separators.

## PhysioEEG 0.5.0

### New Features

- [`eegPlotGlassBrain()`](https://x-biosignal.github.io/PhysioEEG/reference/eegPlotGlassBrain.md)
  creates deterministic sagittal, axial, and coronal maximum-intensity
  projections from finite three-dimensional source coordinates. It uses
  one complete-cloud display transform, a shared type-8
  absolute-amplitude threshold, stable collision tie-breaking, and a
  provenance-labelled schematic outline.
- New source-estimate and beamformer results retain the source
  positions, output name, reduction contract, and coordinate provenance
  needed for glass-brain plotting. Free orientations are summarized by
  root-sum-square within time and RMS across time without changing the
  stored inverse output.

### Internal Improvements

- Source plotting metadata names are reserved so a caller-supplied
  output name cannot silently overwrite source values. Legacy source
  objects without real three-dimensional positions now fail with an
  explicit remediation message.

## PhysioEEG 0.4.0

### New Features

- [`eegPlotSpectrogram()`](https://x-biosignal.github.io/PhysioEEG/reference/eegPlotSpectrogram.md)
  accepts strictly aligned logical, p-value, or inference-result masks.
  Non-selected bins are dimmed without modifying power, and optional
  outlines follow non-interpolated time-frequency cell boundaries.
- Stored
  [`eegSTFT()`](https://x-biosignal.github.io/PhysioEEG/reference/eegSTFT.md),
  [`eegMorletWavelet()`](https://x-biosignal.github.io/PhysioEEG/reference/eegMorletWavelet.md),
  and
  [`eegERSP()`](https://x-biosignal.github.io/PhysioEEG/reference/eegERSP.md)
  products now share one validated plotting path with exact channel and
  axis resolution.

### Internal Improvements

- A private time-frequency inference helper provides complete-family
  threshold correction and sign-separated, four-neighbour cluster-mass
  correction over independent replicate arrays. Sign flips are generated
  per replicate, exact enumeration is used when feasible, and caller RNG
  state is preserved.
- Raw 3D epochs are no longer mistaken for precomputed time-frequency
  data.

## PhysioEEG 0.3.0

### New Features

- [`eegPlotButterflyGFP()`](https://x-biosignal.github.io/PhysioEEG/reference/eegPlotButterflyGFP.md)
  composes exact-order continuous EEG butterfly traces with a Global
  Field Power panel, shared event markers, strict 2D/channel/time
  validation, and explicit population (C) or sample (C - 1)
  denominators.
- [`eegPlotERP()`](https://x-biosignal.github.io/PhysioEEG/reference/eegPlotERP.md)
  now uses the shared colorblind-safe PhysioCore palette and
  [`theme_physio()`](https://x-biosignal.r-universe.dev/PhysioCore/reference/theme_physio.html)
  without changing ERP means or confidence intervals.

### Internal Improvements

- Butterfly plots and microstate analysis share one finite,
  time-by-channel GFP implementation. The population definition
  preserves the established
  [`eegMicrostates()`](https://x-biosignal.github.io/PhysioEEG/reference/eegMicrostates.md)
  result.

## PhysioEEG 0.2.0

Initial release of PhysioEEG as a standalone package in the x-biosignal
ecosystem, split out from the PhysioExperiment monolith. Provides a
comprehensive, end-to-end electroencephalography (EEG) analysis toolkit
operating directly on `PhysioExperiment` objects.

### New Features

#### Preprocessing

- [`eegFilter()`](https://x-biosignal.github.io/PhysioEEG/reference/eegFilter.md)
  performs bandpass, highpass, lowpass, and notch filtering with a
  choice of zero-phase windowed-sinc FIR filters or Butterworth IIR
  filters (`method = "fir"` or `"iir"`).
- [`eegRereference()`](https://x-biosignal.github.io/PhysioEEG/reference/eegRereference.md)
  supports average, mastoid, Cz, and custom reference schemes;
  [`eegMontage()`](https://x-biosignal.github.io/PhysioEEG/reference/eegMontage.md)
  applies 10-20, 10-10, and BioSemi-64 electrode layouts.
- [`eegBadChannels()`](https://x-biosignal.github.io/PhysioEEG/reference/eegBadChannels.md)
  flags flat, noisy, and poorly-correlated channels, and
  [`eegInterpolate()`](https://x-biosignal.github.io/PhysioEEG/reference/eegInterpolate.md)
  reconstructs them via spherical-spline or nearest-neighbour
  interpolation.
- [`eegEpoch()`](https://x-biosignal.github.io/PhysioEEG/reference/eegEpoch.md)
  segments continuous data around events with baseline correction, and
  [`eegArtifactReject()`](https://x-biosignal.github.io/PhysioEEG/reference/eegArtifactReject.md)
  rejects epochs by amplitude threshold, gradient, and related criteria.
- [`eegPreprocess()`](https://x-biosignal.github.io/PhysioEEG/reference/eegPreprocess.md)
  chains filtering, re-referencing, bad-channel handling, and artifact
  rejection into a single configurable pipeline.

#### Independent Component Analysis

- [`eegICA()`](https://x-biosignal.github.io/PhysioEEG/reference/eegICA.md)
  decomposes EEG into independent components using FastICA, Infomax, or
  JADE.
- [`eegICAdetect()`](https://x-biosignal.github.io/PhysioEEG/reference/eegICAdetect.md)
  automatically flags artifactual components (correlation, kurtosis, and
  spatial criteria);
  [`eegICAremove()`](https://x-biosignal.github.io/PhysioEEG/reference/eegICAremove.md)
  and
  [`eegICAmix()`](https://x-biosignal.github.io/PhysioEEG/reference/eegICAmix.md)
  reconstruct cleaned signals.

#### Event-Related Potentials

- [`eegERPdetect()`](https://x-biosignal.github.io/PhysioEEG/reference/eegERPdetect.md)
  locates canonical ERP components (N100, P300, N400, P600, MMN, LPP),
  with
  [`eegERPmeasure()`](https://x-biosignal.github.io/PhysioEEG/reference/eegERPmeasure.md)
  and
  [`eegERPlatency()`](https://x-biosignal.github.io/PhysioEEG/reference/eegERPlatency.md)
  quantifying peak/mean amplitude and fractional-area latency.
- [`eegERPbaseline()`](https://x-biosignal.github.io/PhysioEEG/reference/eegERPbaseline.md),
  [`eegERPdifference()`](https://x-biosignal.github.io/PhysioEEG/reference/eegERPdifference.md),
  and
  [`eegERPgrandAverage()`](https://x-biosignal.github.io/PhysioEEG/reference/eegERPgrandAverage.md)
  handle baseline correction, difference waves, and grand averaging;
  [`eegERPtest()`](https://x-biosignal.github.io/PhysioEEG/reference/eegERPtest.md)
  runs permutation and cluster-based statistics.

#### Time-Frequency Analysis

- [`eegMorletWavelet()`](https://x-biosignal.github.io/PhysioEEG/reference/eegMorletWavelet.md),
  [`eegSTFT()`](https://x-biosignal.github.io/PhysioEEG/reference/eegSTFT.md),
  and
  [`eegMultitaper()`](https://x-biosignal.github.io/PhysioEEG/reference/eegMultitaper.md)
  compute time-frequency representations via Morlet wavelets, short-time
  Fourier transform, and multitaper spectral estimation.
- [`eegERSP()`](https://x-biosignal.github.io/PhysioEEG/reference/eegERSP.md)
  and
  [`eegITC()`](https://x-biosignal.github.io/PhysioEEG/reference/eegITC.md)
  derive event-related spectral perturbation and inter-trial coherence.

#### Connectivity

- [`eegCoherence()`](https://x-biosignal.github.io/PhysioEEG/reference/eegCoherence.md),
  [`eegPLV()`](https://x-biosignal.github.io/PhysioEEG/reference/eegPLV.md),
  and
  [`eegWPLI()`](https://x-biosignal.github.io/PhysioEEG/reference/eegWPLI.md)
  estimate magnitude-squared and imaginary coherence, phase locking
  value, and weighted phase lag index.
- [`eegGrangerCausality()`](https://x-biosignal.github.io/PhysioEEG/reference/eegGrangerCausality.md)
  computes spectral Granger causality, and
  [`eegConnectivityMatrix()`](https://x-biosignal.github.io/PhysioEEG/reference/eegConnectivityMatrix.md)
  assembles all-pairs connectivity matrices.

#### Source Localization

- [`eegForwardModel()`](https://x-biosignal.github.io/PhysioEEG/reference/eegForwardModel.md)
  builds spherical and simplified-BEM head models.
- [`eegSourceEstimate()`](https://x-biosignal.github.io/PhysioEEG/reference/eegSourceEstimate.md)
  reconstructs cortical activity via sLORETA, eLORETA, and minimum-norm
  estimation;
  [`eegBeamformer()`](https://x-biosignal.github.io/PhysioEEG/reference/eegBeamformer.md)
  provides an LCMV spatial filter, and
  [`eegSourcePower()`](https://x-biosignal.github.io/PhysioEEG/reference/eegSourcePower.md)
  summarizes band-limited source power.

#### Microstates

- [`eegMicrostates()`](https://x-biosignal.github.io/PhysioEEG/reference/eegMicrostates.md)
  clusters topographies (k-means, AAHC, PCA) into microstate maps, with
  [`eegMicrostateBackfit()`](https://x-biosignal.github.io/PhysioEEG/reference/eegMicrostateBackfit.md),
  [`eegMicrostateStats()`](https://x-biosignal.github.io/PhysioEEG/reference/eegMicrostateStats.md),
  and
  [`eegMicrostateSequence()`](https://x-biosignal.github.io/PhysioEEG/reference/eegMicrostateSequence.md)
  for back-fitting, statistics, and sequence extraction.

#### Sleep Analysis

- [`eegSleepStage()`](https://x-biosignal.github.io/PhysioEEG/reference/eegSleepStage.md)
  classifies epochs into Wake/N1/N2/N3/REM using simplified AASM
  criteria.
- [`eegSpindleDetect()`](https://x-biosignal.github.io/PhysioEEG/reference/eegSpindleDetect.md),
  [`eegKcomplexDetect()`](https://x-biosignal.github.io/PhysioEEG/reference/eegKcomplexDetect.md),
  and
  [`eegSlowWaveDetect()`](https://x-biosignal.github.io/PhysioEEG/reference/eegSlowWaveDetect.md)
  detect sleep spindles, K-complexes, and slow waves;
  [`eegSleepMetrics()`](https://x-biosignal.github.io/PhysioEEG/reference/eegSleepMetrics.md)
  summarizes sleep architecture.

#### Brain-Computer Interface

- [`eegCSP()`](https://x-biosignal.github.io/PhysioEEG/reference/eegCSP.md)
  computes common spatial patterns;
  [`eegSSVEP()`](https://x-biosignal.github.io/PhysioEEG/reference/eegSSVEP.md)
  and
  [`eegMotorImagery()`](https://x-biosignal.github.io/PhysioEEG/reference/eegMotorImagery.md)
  extract steady-state and motor-imagery features.
- [`eegBCIfeatures()`](https://x-biosignal.github.io/PhysioEEG/reference/eegBCIfeatures.md)
  and
  [`eegBCIclassify()`](https://x-biosignal.github.io/PhysioEEG/reference/eegBCIclassify.md)
  provide band-power, CSP, and Riemannian feature extraction with
  classification.

#### Clinical EEG

- [`eegSpikeDetect()`](https://x-biosignal.github.io/PhysioEEG/reference/eegSpikeDetect.md)
  detects epileptiform spikes (morphology and template methods);
  [`eegQEEG()`](https://x-biosignal.github.io/PhysioEEG/reference/eegQEEG.md)
  produces quantitative band-power profiles.
- [`eegAsymmetry()`](https://x-biosignal.github.io/PhysioEEG/reference/eegAsymmetry.md),
  [`eegSlowing()`](https://x-biosignal.github.io/PhysioEEG/reference/eegSlowing.md),
  and
  [`eegSuppression()`](https://x-biosignal.github.io/PhysioEEG/reference/eegSuppression.md)
  compute hemispheric asymmetry indices, focal/diffuse slowing, and
  burst-suppression detection.

#### Visualization

- Plotting helpers cover signal traces
  ([`eegPlotSignal()`](https://x-biosignal.github.io/PhysioEEG/reference/eegPlotSignal.md)),
  ERP waveforms
  ([`eegPlotERP()`](https://x-biosignal.github.io/PhysioEEG/reference/eegPlotERP.md)),
  topographic maps
  ([`eegPlotTopomap()`](https://x-biosignal.github.io/PhysioEEG/reference/eegPlotTopomap.md),
  [`eegPlotTopomapSeries()`](https://x-biosignal.github.io/PhysioEEG/reference/eegPlotTopomapSeries.md)),
  spectrograms
  ([`eegPlotSpectrogram()`](https://x-biosignal.github.io/PhysioEEG/reference/eegPlotSpectrogram.md)),
  connectivity
  ([`eegPlotConnectivity()`](https://x-biosignal.github.io/PhysioEEG/reference/eegPlotConnectivity.md)),
  hypnograms
  ([`eegPlotHypnogram()`](https://x-biosignal.github.io/PhysioEEG/reference/eegPlotHypnogram.md)),
  ICA components
  ([`eegPlotICA()`](https://x-biosignal.github.io/PhysioEEG/reference/eegPlotICA.md)),
  and source maps
  ([`eegPlotSource()`](https://x-biosignal.github.io/PhysioEEG/reference/eegPlotSource.md)).

#### Simulated Data and Onboarding

- [`make_eeg()`](https://x-biosignal.github.io/PhysioEEG/reference/make_eeg.md),
  [`make_eeg_erp()`](https://x-biosignal.github.io/PhysioEEG/reference/make_eeg_erp.md),
  [`make_eeg_sleep()`](https://x-biosignal.github.io/PhysioEEG/reference/make_eeg_sleep.md),
  [`make_eeg_bci()`](https://x-biosignal.github.io/PhysioEEG/reference/make_eeg_bci.md),
  and
  [`make_eeg_spikes()`](https://x-biosignal.github.io/PhysioEEG/reference/make_eeg_spikes.md)
  generate realistic synthetic EEG for testing and demonstration.
- [`eegQuickStart()`](https://x-biosignal.github.io/PhysioEEG/reference/eegQuickStart.md)
  provides ready-to-run example workflows.

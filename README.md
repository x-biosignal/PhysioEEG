# PhysioEEG <img src="man/figures/logo.png" align="right" height="139" alt="PhysioEEG logo" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/x-biosignal/PhysioEEG/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/x-biosignal/PhysioEEG/actions/workflows/R-CMD-check.yaml)
[![CRAN status](https://www.r-pkg.org/badges/version/PhysioEEG)](https://CRAN.R-project.org/package=PhysioEEG)
[![r-universe](https://x-biosignal.r-universe.dev/badges/PhysioEEG)](https://x-biosignal.r-universe.dev/PhysioEEG)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

**EEG Analysis Functions for PhysioExperiment Objects**

PhysioEEG provides 79 exported functions covering the full
electroencephalography analysis pipeline. Built on top of PhysioCore, it
supports preprocessing, independent component analysis (ICA), event-related
potential (ERP) extraction and measurement, source localization, microstate
analysis, sleep staging, brain-computer interface (BCI) feature extraction,
clinical EEG analysis, time-frequency decomposition, connectivity analysis,
and publication-quality visualization -- all operating directly on
`PhysioExperiment` objects.

## Installation

You can install PhysioEEG from [r-universe](https://x-biosignal.r-universe.dev):

```r
install.packages("PhysioEEG",
  repos = c("https://x-biosignal.r-universe.dev", "https://cloud.r-project.org"))
```

Or install the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("x-biosignal/PhysioEEG")
```

## Quick Start

```r
library(PhysioEEG)

# Generate simulated EEG data with ERP components
pe <- make_eeg_erp(n_epochs = 40, n_channels = 19, sr = 250)

# Run full preprocessing pipeline (filter, re-reference, artifact rejection)
pe <- eegPreprocess(pe, lowcut = 1, highcut = 40, ref = "average")

# Detect and measure P300 component
erp <- eegERPdetect(pe, component = "P300")
measures <- eegERPmeasure(pe, component = "P300")

# Plot ERP waveform and topographic map
eegPlotERP(pe, channels = c("Fz", "Cz", "Pz"))
eegPlotTopomap(pe, time = 0.35)

# Inspect continuous channel traces and reference-dependent GFP together
continuous <- make_eeg(n_time = 2500, n_channels = 8, sr = 500)
eegPlotButterflyGFP(continuous, channels = c("Fz", "F3", "F4"))
```

## Features

### Preprocessing Pipeline

A complete, configurable preprocessing pipeline from raw recordings to analysis-ready data:

- `eegPreprocess()` -- full pipeline with configurable steps
- `eegFilter()` -- bandpass, highpass, lowpass, and notch filtering
- `eegMontage()` -- apply standard electrode montages (10-20, 10-10)
- `eegRereference()` -- re-reference to average, linked mastoids, or custom reference
- `eegBadChannels()` -- automatic bad channel detection by correlation, variance, and spectral criteria
- `eegInterpolate()` -- spherical spline interpolation of bad channels
- `eegEpoch()` -- segment continuous data into event-locked epochs
- `eegArtifactReject()` -- reject epochs by amplitude threshold, gradient, or statistical criteria

### Independent Component Analysis (ICA)

Blind source separation for artifact removal and source identification:

- `eegICA()` -- decompose signals using FastICA, Infomax, or JADE algorithms
- `eegICAdetect()` -- automatic classification of artifact components (eye blinks, saccades, muscle, cardiac)
- `eegICAmix()` -- inspect mixing and unmixing matrices
- `eegICAremove()` -- remove selected components and reconstruct clean signals
- `eegPlotICA()` -- visualize component topographies, time courses, and spectra

### Event-Related Potential (ERP) Analysis

End-to-end ERP component extraction, measurement, and statistical testing:

- `eegERPdetect()` -- detect standard components (N100, P300, N400, P600, MMN) with adaptive time windows
- `eegERPmeasure()` -- measure peak amplitude, peak latency, mean amplitude, and area
- `eegERPlatency()` -- fractional area latency and onset latency estimation
- `eegERPbaseline()` -- apply baseline correction with configurable windows
- `eegERPdifference()` -- compute difference waveforms between conditions
- `eegERPtest()` -- point-by-point and cluster-based permutation tests
- `eegERPgrandAverage()` -- compute grand average across subjects or sessions
- `eegPlotERP()` -- plot ERP waveforms with confidence intervals and condition overlays

### Source Localization

Estimate the cortical sources underlying scalp EEG:

- `eegForwardModel()` -- compute forward model with boundary element method (BEM)
- `eegSourceEstimate()` -- distributed source estimation with eLORETA or sLORETA
- `eegBeamformer()` -- LCMV beamformer for focal source localization
- `eegSourcePower()` -- compute source power maps for frequency bands
- `eegPlotSource()` -- visualize legacy or explicit two-dimensional source maps
- `eegPlotGlassBrain()` -- project real 3D source coordinates into deterministic
  sagittal, axial, and coronal maximum-intensity views

`eegPlotGlassBrain()` requires finite three-dimensional source coordinates; it
never turns a bare amplitude vector into anatomical positions. Its shared
type-8 threshold is computed from absolute amplitude before projection.
Source-estimate orientation components are reduced by root-sum-square within
time and RMS across time. The bundled outline is a schematic display frame,
not MRI registration, patient anatomy, or validation of the inverse solution.

### Microstate Analysis

Characterize the temporal dynamics of global brain states:

- `eegMicrostates()` -- segment EEG into microstates using K-means or AAHC (atomize and agglomerate hierarchical clustering)
- `eegMicrostateBackfit()` -- backfit microstate labels to continuous EEG
- `eegMicrostateSequence()` -- extract microstate transition sequences
- `eegMicrostateStats()` -- compute duration, occurrence, coverage, and transition probabilities

### Sleep Staging and Event Detection

Automated sleep analysis following AASM guidelines:

- `eegSleepStage()` -- automatic sleep staging according to AASM criteria (W, N1, N2, N3, REM)
- `eegSpindleDetect()` -- detect sleep spindles (sigma band bursts)
- `eegKcomplexDetect()` -- detect K-complexes
- `eegSlowWaveDetect()` -- detect slow-wave activity (delta oscillations)
- `eegSleepMetrics()` -- compute sleep efficiency, latency, WASO, and stage percentages
- `eegPlotHypnogram()` -- plot hypnogram with detected events overlay

### Brain-Computer Interface (BCI) Features

Feature extraction and classification for BCI paradigms:

- `eegBCIfeatures()` -- extract multi-domain feature vectors (time, frequency, spatial)
- `eegCSP()` -- common spatial patterns for motor imagery discrimination
- `eegSSVEP()` -- steady-state visually evoked potential detection and classification
- `eegMotorImagery()` -- motor imagery feature extraction (mu/beta ERD/ERS)
- `eegBCIclassify()` -- classify BCI features with LDA, SVM, or random forest

### Clinical EEG

Functions for clinical neurophysiology and QEEG:

- `eegSpikeDetect()` -- detect epileptiform spikes and sharp waves
- `eegQEEG()` -- quantitative EEG analysis (absolute/relative power, peak frequency)
- `eegAsymmetry()` -- compute inter-hemispheric asymmetry indices
- `eegSlowing()` -- detect and quantify EEG slowing (theta/alpha ratio)
- `eegSuppression()` -- detect burst-suppression patterns

### Time-Frequency Analysis

Spectral decomposition with multiple methods:

- `eegMorletWavelet()` -- continuous wavelet transform with Morlet wavelets
- `eegSTFT()` -- short-time Fourier transform
- `eegMultitaper()` -- multitaper spectral estimation (DPSS tapers)
- `eegERSP()` -- event-related spectral perturbation
- `eegITC()` -- inter-trial coherence (phase-locking across trials)
- `eegPlotSpectrogram()` -- plot validated STFT, Morlet, ERSP, or legacy FFT
  grids with strictly aligned display masks

`eegPlotSpectrogram()` keeps numeric power unchanged when a mask is supplied:
non-selected bins are dimmed with display alpha and logical boundaries are
drawn without interpolating the mask. Logical masks and p-value matrices must
match the complete pre-filter time and frequency axes exactly. Cluster masks
require repeated, independent and exchangeable observations; the correction
is for the maximum sign-separated four-neighbour cluster mass, so a corrected
cluster is not bin-wise or mechanistic evidence.

### Connectivity Analysis

Functional and effective connectivity between EEG channels:

- `eegCoherence()` -- magnitude-squared coherence
- `eegPLV()` -- phase-locking value
- `eegWPLI()` -- weighted phase lag index (robust to volume conduction)
- `eegGrangerCausality()` -- Granger causality for directed connectivity
- `eegConnectivityMatrix()` -- compute full connectivity matrices for any metric
- `eegPlotConnectivity()` -- legacy connectivity heatmaps and circle plots
- `eegPlotConnectogram()` -- deterministic directed/undirected circular views
  with anatomical or cluster ordering, module arcs, and optional edge bundling.
  Its strict absolute threshold is descriptive and does not imply corrected
  statistical significance.

### Visualization

Publication-quality plots for every analysis stage:

- `eegPlotSignal()` -- multi-channel signal traces with event markers
- `eegPlotButterflyGFP()` -- exact-order butterfly traces with a shared-time
  Global Field Power panel (population or sample denominator)
- `eegPlotERP()` -- ERP waveforms with confidence intervals
- `eegPlotTopomap()` -- topographic scalp maps with interpolation
- `eegPlotTopomapSeries()` -- series of topomaps across time points
- `eegPlotSpectrogram()` -- masked time-frequency grids with strict axis
  alignment and non-destructive opacity
- `eegPlotConnectivity()` -- connectivity matrices and network graphs
- `eegPlotConnectogram()` -- inspectable circular connectivity views preserving
  matrix direction, signed values, channel identity, and estimator metadata
- `eegPlotICA()` -- ICA component topographies and time courses
- `eegPlotSource()` -- legacy and explicit two-dimensional source maps
- `eegPlotGlassBrain()` -- three-view source maximum-intensity projections
- `eegPlotHypnogram()` -- sleep hypnograms with event annotations

### Simulated Data Generators

Ready-to-use data for testing, demonstration, and teaching:

- `make_eeg()` -- multi-channel EEG with realistic spectral properties
- `make_eeg_erp()` -- epoched data with embedded ERP components (N100, P300)
- `make_eeg_sleep()` -- polysomnography data with sleep stages and spindles
- `make_eeg_bci()` -- motor imagery data with lateralized mu/beta patterns
- `make_eeg_spikes()` -- EEG with embedded epileptiform discharges

## Dependencies

- **R** (>= 4.2)
- **[PhysioCore](https://github.com/x-biosignal/PhysioCore)**
- **SummarizedExperiment**
- **S4Vectors**
- **stats**

## PhysioExperiment Ecosystem

PhysioEEG is the EEG analysis layer of the PhysioExperiment ecosystem, a suite of R packages for multi-modal physiological signal analysis:

| Package | Description |
|---------|-------------|
| [PhysioCore](https://github.com/x-biosignal/PhysioCore) | Core data structures and accessors |
| [PhysioIO](https://github.com/x-biosignal/PhysioIO) | File I/O (EDF, HDF5, BIDS, CSV, MAT) |
| [PhysioPreprocess](https://github.com/x-biosignal/PhysioPreprocess) | Preprocessing (filters, ICA, resampling) |
| [PhysioAnalysis](https://github.com/x-biosignal/PhysioAnalysis) | Analysis and visualization |
| **PhysioEEG** | EEG analysis (ICA, ERP, source, BCI, sleep) |
| [PhysioEMG](https://github.com/x-biosignal/PhysioEMG) | EMG analysis (synergy, fatigue, onset) |
| [PhysioECG](https://github.com/x-biosignal/PhysioECG) | ECG and HRV analysis |

Visit the [r-universe page](https://x-biosignal.r-universe.dev) to browse all available packages.

## License

MIT License. See [LICENSE](LICENSE) for details.

## Author

Yusuke Matsui

## Governance & support

Part of the [Physio ecosystem](https://x-biosignal.r-universe.dev). Community and
policy documents live in the umbrella repository:

- [Code of Conduct](https://github.com/x-biosignal/PhysioExperiment/blob/main/CODE_OF_CONDUCT.md)
- [Contributing](https://github.com/x-biosignal/PhysioExperiment/blob/main/CONTRIBUTING.md)
- [Governance](https://github.com/x-biosignal/PhysioExperiment/blob/main/GOVERNANCE.md)
- [Support](https://github.com/x-biosignal/PhysioExperiment/blob/main/SUPPORT.md)
- [Security policy](https://github.com/x-biosignal/PhysioExperiment/blob/main/SECURITY.md)
- [Deprecation & lifecycle policy](https://github.com/x-biosignal/PhysioExperiment/blob/main/DEPRECATION.md)

# PhysioEEG: EEG Analysis Functions for PhysioExperiment Objects

Provides comprehensive electroencephalography (EEG) analysis functions
for PhysioExperiment objects. Includes preprocessing (filtering,
re-referencing, bad channel detection, interpolation, epoching, artifact
rejection), independent component analysis (FastICA, Infomax, JADE) with
automatic artifact detection, event-related potential (ERP) component
detection and measurement (N100, P300, N400, P600, MMN), source
localization (eLORETA, sLORETA, LCMV beamformer), EEG microstate
analysis (K-means, AAHC), sleep staging (AASM criteria,
spindle/K-complex/slow-wave detection), brain-computer interface
features (CSP, SSVEP, motor imagery), clinical EEG analysis (spike
detection, QEEG, asymmetry indices), time-frequency analysis (Morlet
wavelet, STFT, multitaper, ERSP, ITC), connectivity analysis (coherence,
PLV, wPLI, Granger causality), visualization (signal traces, ERP
waveforms, topographic maps, spectrograms, connectivity plots,
hypnograms, source maps), and simulated data generators for testing and
demonstration.

## See also

Useful links:

- <https://github.com/x-biosignal/PhysioEEG>

- <https://x-biosignal.r-universe.dev/PhysioEEG>

- <https://x-biosignal.github.io/PhysioEEG/>

- Report bugs at <https://github.com/x-biosignal/PhysioEEG/issues>

## Author

**Maintainer**: Yusuke Matsui <mail.to.matsui@gmail.com>

## Examples

``` r
# Simulate EEG, band-pass filter, and summarise signal complexity
pe <- make_eeg(n_time = 2500, n_channels = 8, sr = 250)
pe <- eegFilter(pe, lowcut = 1, highcut = 40)
eegComplexity(pe, measures = c("hjorth_mobility", "spectral_entropy"))
#>   channel hjorth_mobility spectral_entropy
#> 1     Fp1      0.10689380        0.2315440
#> 2     Fp2      0.07422128        0.1377136
#> 3      F7      0.08172641        0.3753449
#> 4      F3      0.11333633        0.3157096
#> 5      Fz      0.07352278        0.1915594
#> 6      F4      0.14445514        0.3461614
#> 7      F8      0.07009828        0.2161299
#> 8      T3      0.05397000        0.1489548
```

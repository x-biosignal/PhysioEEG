# ERP Analysis with PhysioEEG

## ERP Analysis with PhysioEEG

This vignette demonstrates a complete ERP analysis workflow using
PhysioEEG.

### Simulating an Oddball Experiment

``` r

library(PhysioEEG)

# Create oddball paradigm data: 40 epochs, target and standard conditions
pe <- make_eeg_erp(n_epochs = 40, n_channels = 19, sr = 250, epoch_sec = 1.0)

# Check conditions
conditions <- metadata(pe)$conditions
table(conditions)
```

### ERP Component Detection

``` r

# Detect N100 (early sensory processing)
n100 <- eegERPdetect(pe, component = "N100")
print(n100)

# Detect P300 (attention/decision making)
p300 <- eegERPdetect(pe, component = "P300")
print(p300)
```

### ERP Measurement Methods

``` r

# Peak amplitude measurement
peak <- eegERPmeasure(pe, window = c(250, 500), method = "peak",
                      polarity = "positive")

# Mean amplitude measurement (more robust to noise)
mean_amp <- eegERPmeasure(pe, window = c(250, 500), method = "mean",
                          polarity = "positive")

# Adaptive mean (window centered on peak)
adaptive <- eegERPmeasure(pe, window = c(250, 500), method = "adaptive_mean",
                          polarity = "positive")

print(data.frame(method = c("peak", "mean", "adaptive_mean"),
                 amplitude = c(peak$amplitude[1], mean_amp$amplitude[1],
                               adaptive$amplitude[1])))
```

### Fractional Area Latency

``` r

# 50% fractional area latency (more reliable than peak latency)
lat <- eegERPlatency(pe, window = c(250, 500), fraction = 0.5,
                     polarity = "positive")
print(lat)
```

### Difference Waves

``` r

# Separate target and standard epochs (requires indexing by condition)
# Compute difference waveform
pe_diff <- eegERPdifference(pe, pe)  # Simplified example
```

### Artifact Removal with ICA Before ERP

``` r

# For real data, apply ICA before epoching:
# 1. Run ICA on continuous data
# 2. Remove artifact components
# 3. Epoch the cleaned data
# 4. Average and measure ERPs
pe_continuous <- make_eeg(n_time = 10000, n_channels = 19, sr = 250)
pe_ica <- eegICA(pe_continuous, method = "fastica")
artifacts <- eegICAdetect(pe_ica, method = "kurtosis")
pe_clean <- eegICAremove(pe_ica,
                         components = artifacts$component[artifacts$type == "artifact"])
```

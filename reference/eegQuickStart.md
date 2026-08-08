# PhysioEEG Quick Start Guide

Prints a guided walkthrough with runnable code examples for selected EEG
analysis workflows. Covers ERP analysis, sleep staging, BCI
classification, and connectivity analysis.

## Usage

``` r
eegQuickStart(workflow = "all")
```

## Arguments

- workflow:

  Character string: one of `"erp"`, `"sleep"`, `"bci"`,
  `"connectivity"`, or `"all"` (default: `"all"`).

## Value

Invisibly returns `NULL`. Prints guide to console.

## See also

[`make_eeg()`](https://x-biosignal.github.io/PhysioEEG/reference/make_eeg.md),
[`make_eeg_erp()`](https://x-biosignal.github.io/PhysioEEG/reference/make_eeg_erp.md),
[`make_eeg_sleep()`](https://x-biosignal.github.io/PhysioEEG/reference/make_eeg_sleep.md),
[`make_eeg_bci()`](https://x-biosignal.github.io/PhysioEEG/reference/make_eeg_bci.md)

## Examples

``` r
eegQuickStart("erp")
#> === PhysioEEG Quick Start Guide ===
#> 
#> --- ERP Analysis Workflow ---
#> # 1. Create or load EEG data
#> pe <- make_eeg_erp(n_epochs = 40, n_channels = 19, sr = 250)
#> 
#> # 2. Preprocess
#> pe <- eegFilter(pe, lowcut = 0.1, highcut = 30)
#> pe <- eegRereference(pe, ref_type = 'average')
#> 
#> # 3. Detect and measure ERP components
#> components <- eegERPdetect(pe)
#> measurements <- eegERPmeasure(pe, components = 'P300')
#> 
#> # 4. Statistical testing
#> results <- eegERPtest(pe, conditions = metadata(pe)$conditions)
#> 
#> # 5. Visualize
#> # eegPlotERP(pe, channels = c('Fz', 'Cz', 'Pz'))
#> 
```

# Extract ICLabel features for one component (internal)

Extract ICLabel features for one component (internal)

## Usage

``` r
.iclabel_features(topo, act, sr, ch_labels, line_freq = 50)
```

## Arguments

- topo:

  Channel weights (a column of the mixing matrix).

- act:

  Component activation (time course).

- sr:

  Sampling rate in Hz.

- ch_labels:

  Channel labels.

- line_freq:

  Mains line frequency in Hz.

## Value

A named numeric vector of interpretable features.

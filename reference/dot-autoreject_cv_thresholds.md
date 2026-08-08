# Cross-validated per-channel peak-to-peak thresholds (internal)

For each channel, picks the threshold minimising the mean over folds of
RMSE( mean of good training epochs, median of test epochs ).

## Usage

``` r
.autoreject_cv_thresholds(data, ptp, grid, folds)
```

## Arguments

- data:

  3D array (time x channels x epochs).

- ptp:

  Channels x epochs peak-to-peak matrix.

- grid:

  Candidate threshold vector.

- folds:

  Integer fold assignment per epoch.

## Value

Numeric vector of per-channel thresholds.

# Cross-validated selection of consensus and n_interpolate (internal)

Grid-searches (consensus, n_interpolate), scoring each by the mean over
folds of RMSE( mean of cleaned training epochs, median of test epochs ).

## Usage

``` r
.autoreject_select_consensus(data, ptp, tau, pos, kappa_grid, rho_grid, folds)
```

## Arguments

- data:

  3D array (time x channels x epochs).

- ptp:

  Channels x epochs peak-to-peak matrix.

- tau:

  Per-channel thresholds.

- pos:

  Electrode positions matrix.

- kappa_grid:

  Candidate consensus fractions.

- rho_grid:

  Candidate interpolation counts.

- folds:

  Integer fold assignment per epoch.

## Value

List with `consensus`, `n_interpolate`, and `err`.

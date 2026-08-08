# Fit CSP spatial filters on a set of training trials

Fit CSP spatial filters on a set of training trials

## Usage

``` r
.csp_fit(data, idx, labels, n_filters = 3)
```

## Arguments

- data:

  3D array (time x channels x trials).

- idx:

  Training trial indices.

- labels:

  Class labels for `idx` (length == `length(idx)`).

- n_filters:

  Number of CSP filter pairs (default: 3).

## Value

A `(2 n_filters) x channels` spatial filter matrix.

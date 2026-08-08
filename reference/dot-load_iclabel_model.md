# Load the shipped ICLabel softmax weight table

Reads the weight table from `inst/extdata/iclabel_weights.csv` (a `term`
column plus one column per class) and caches it. Falls back to an
identical built-in table when the file is unavailable.

## Usage

``` r
.load_iclabel_model()
```

## Value

A data.frame with a `term` column and one column per class.

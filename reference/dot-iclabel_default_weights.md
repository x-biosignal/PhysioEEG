# Built-in ICLabel softmax weights (fallback)

Interpretable multinomial-logistic weights, identical to the shipped
CSV, used when the packaged file cannot be read.

## Usage

``` r
.iclabel_default_weights()
```

## Value

A data.frame with a `term` column and one column per class.

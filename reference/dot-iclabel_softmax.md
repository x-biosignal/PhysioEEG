# Softmax class probabilities from design terms and weight table (internal)

Softmax class probabilities from design terms and weight table
(internal)

## Usage

``` r
.iclabel_softmax(terms, model, classes)
```

## Arguments

- terms:

  Named design-term vector from
  [`.iclabel_terms()`](https://x-biosignal.github.io/PhysioEEG/reference/dot-iclabel_terms.md).

- model:

  Weight table from
  [`.load_iclabel_model()`](https://x-biosignal.github.io/PhysioEEG/reference/dot-load_iclabel_model.md).

- classes:

  Class order.

## Value

Named numeric vector of probabilities summing to 1.

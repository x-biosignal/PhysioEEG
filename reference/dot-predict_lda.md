# Predict using a trained LDA model

Predict using a trained LDA model

## Usage

``` r
.predict_lda(model, features)
```

## Arguments

- model:

  A list as returned by `.train_lda`.

- features:

  Numeric matrix (n_samples x n_features).

## Value

A list with `class` (character vector) and `confidence` (numeric
vector).

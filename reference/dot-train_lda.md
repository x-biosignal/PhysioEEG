# Train Fisher's LDA model

Train Fisher's LDA model

## Usage

``` r
.train_lda(features, labels, method = "lda")
```

## Arguments

- features:

  Numeric matrix (n_trials x n_features).

- labels:

  Character vector of class labels (must have exactly 2 classes).

- method:

  `"lda"` or `"shrinkage_lda"`.

## Value

A list with `weights`, `threshold`, `classes`, and `class_means`.

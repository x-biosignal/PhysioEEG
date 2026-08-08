# FastICA algorithm (internal)

FastICA algorithm (internal)

## Usage

``` r
.fastica(whitened, n_components, max_iter, tol)
```

## Arguments

- whitened:

  Whitened data matrix (n_time x n_components).

- n_components:

  Number of components.

- max_iter:

  Maximum iterations.

- tol:

  Convergence tolerance.

## Value

Unmixing matrix W (n_components x n_components).

# Infomax ICA algorithm (internal)

Natural gradient Infomax ICA. Updates W using the rule: dW = (I + (1 -
2\*sigmoid(u)) \* t(u)) \* W where u = W \* x.

## Usage

``` r
.infomax(whitened, n_components, max_iter, tol)
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

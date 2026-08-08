# Iterative Frechet (Riemannian geometric) mean of SPD matrices

Computes the geometric mean on the manifold of symmetric positive
definite (SPD) matrices using the iterative fixed-point algorithm.
Converges to the true Frechet mean in the Riemannian metric.

## Usage

``` r
.geometric_mean_spd(covs, max_iter = 50, tol = 1e-08)
```

## Arguments

- covs:

  List of SPD matrices (all same dimensions).

- max_iter:

  Maximum iterations (default: 50).

- tol:

  Convergence tolerance on the tangent-space residual (default: 1e-8).

## Value

The geometric mean SPD matrix.

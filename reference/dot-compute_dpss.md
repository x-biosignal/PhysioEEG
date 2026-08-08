# Compute DPSS (Slepian) tapers via tridiagonal eigendecomposition

Approximates Discrete Prolate Spheroidal Sequences by constructing a
symmetric tridiagonal matrix whose eigenvectors corresponding to the
largest eigenvalues are the desired tapers.

## Usage

``` r
.compute_dpss(n, nw, k)
```

## Arguments

- n:

  Integer length of the data.

- nw:

  Numeric time-half-bandwidth product.

- k:

  Integer number of tapers to return.

## Value

Matrix of dimension n x k containing the DPSS tapers as columns.

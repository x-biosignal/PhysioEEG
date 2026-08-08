# Compute Legendre polynomial of degree n at x

Uses the Bonnet recursion formula to compute the Legendre polynomial
P_n(x) for a given degree n and evaluation point(s) x.

## Usage

``` r
.legendre_poly(n, x)
```

## Arguments

- n:

  Non-negative integer degree.

- x:

  Numeric vector of evaluation points in `[-1, 1]`.

## Value

Numeric vector of P_n(x).

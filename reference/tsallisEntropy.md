# Tsallis Entropy

Computes the Tsallis entropy of order `q` of a one-dimensional signal
(Tsallis, 1988) – the *non-additive* one-parameter generalization of the
Shannon entropy of the signal's amplitude distribution, complementary to
the (additive)
[`renyiEntropy`](https://x-biosignal.github.io/PhysioEEG/reference/renyiEntropy.md).
The signal is discretized into `bins` equal-width bins to form a
probability distribution \\p\\, and \\S_q = \frac{1 - \sum_i p_i^q}{q -
1}\\ is returned (equivalently \\\sum_i p_i \ln_q(1/p_i)\\ with the
q-logarithm). \\q \to 1\\ recovers the Shannon entropy (in nats); unlike
Renyi, Tsallis is non-additive (\\S_q(A,B) = S_q(A) + S_q(B) + (1-q)
S_q(A) S_q(B)\\ for independent \\A, B\\), so Renyi and Tsallis agree
only at \\q = 1\\.

## Usage

``` r
tsallisEntropy(x, q = 2, bins = 16L)
```

## Arguments

- x:

  A numeric vector (the time series).

- q:

  Order of the Tsallis entropy (default 2); `q = 1` returns the Shannon
  entropy.

- bins:

  Number of equal-width bins for the amplitude histogram (default 16).

## Value

A single numeric value, the Tsallis entropy.

## References

Tsallis, C. (1988). "Possible generalization of Boltzmann-Gibbs
statistics." *Journal of Statistical Physics*, 52(1-2), 479–487.
[doi:10.1007/BF01016429](https://doi.org/10.1007/BF01016429)

## See also

[`renyiEntropy`](https://x-biosignal.github.io/PhysioEEG/reference/renyiEntropy.md)
for the additive generalization,
[`svdEntropy`](https://x-biosignal.github.io/PhysioEEG/reference/svdEntropy.md),
[`petrosianFD`](https://x-biosignal.github.io/PhysioEEG/reference/petrosianFD.md).

## Examples

``` r
set.seed(1)
tsallisEntropy(rnorm(1000), q = 2)   # non-additive (collision-like)
#> [1] 0.883062
tsallisEntropy(rnorm(1000), q = 1)   # = Shannon entropy
#> [1] 2.300707
```

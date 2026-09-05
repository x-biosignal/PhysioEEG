# Renyi Entropy

Computes the Renyi entropy of order `alpha` of a one-dimensional signal
(Renyi, 1961) – a one-parameter generalization of the Shannon entropy of
the signal's amplitude distribution. The signal is discretized into
`bins` equal-width bins to form a probability distribution \\p\\, and
\\H\_\alpha = \frac{1}{1-\alpha}\ln\left(\sum_i p_i^\alpha\right)\\ is
returned in nats. Special cases: \\\alpha \to 1\\ is the Shannon
entropy, \\\alpha = 2\\ the collision entropy, and \\\alpha \to \infty\\
the min-entropy. \\H\_\alpha\\ is non-increasing in \\\alpha\\.

## Usage

``` r
renyiEntropy(x, alpha = 2, bins = 16L)
```

## Arguments

- x:

  A numeric vector (the time series).

- alpha:

  Order of the Renyi entropy (default 2, the collision entropy);
  `alpha = 1` returns the Shannon entropy.

- bins:

  Number of equal-width bins for the amplitude histogram (default 16).

## Value

A single numeric value, the Renyi entropy in nats.

## References

Renyi, A. (1961). "On measures of entropy and information." In
*Proceedings of the Fourth Berkeley Symposium on Mathematical Statistics
and Probability*, Volume 1, 547–561.

## See also

[`svdEntropy`](https://x-biosignal.github.io/PhysioEEG/reference/svdEntropy.md),
[`petrosianFD`](https://x-biosignal.github.io/PhysioEEG/reference/petrosianFD.md),
[`eegComplexity`](https://x-biosignal.github.io/PhysioEEG/reference/eegComplexity.md)
for other complexity measures.

## Examples

``` r
set.seed(1)
renyiEntropy(rnorm(1000), alpha = 2)   # collision entropy
#> [1] 2.146111
renyiEntropy(rnorm(1000), alpha = 1)   # = Shannon entropy
#> [1] 2.300707
```

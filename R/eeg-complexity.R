# EEG complexity / nonlinear dynamics.
#
# Per-channel entropy, fractal-dimension and detrended-fluctuation measures that
# were missing from the ecosystem's EEG stack. The numeric cores here are
# domain-neutral (they operate on a plain numeric vector) and are candidates for
# future consolidation into PhysioCore alongside PhysioMoCap's sampleEntropy /
# approximateEntropy; they are implemented here because PhysioEEG must not depend
# on a sibling domain package. Entropy measures are O(N^2), so long channels are
# capped (with a warning) via `max_samples`.

# --- numeric cores (operate on a numeric vector x) --------------------------

# Chebyshev-distance template match count for sample entropy (fast inner loop).
.eeg_c_matchcount <- function(x, mm, r_abs, nv) {
  # nv length-mm templates starting at 1..nv
  M <- vapply(seq_len(mm), function(j) x[j:(j + nv - 1L)], numeric(nv))  # nv x mm
  tot <- 0
  for (i in seq_len(nv - 1L)) {
    rows <- (i + 1L):nv
    dmax <- abs(M[rows, 1L] - M[i, 1L])
    if (mm > 1L) for (jj in 2:mm) dmax <- pmax(dmax, abs(M[rows, jj] - M[i, jj]))
    tot <- tot + sum(dmax <= r_abs)
  }
  tot
}

# Sample entropy (Richman & Moorman 2000). r is a fraction of SD unless r_abs set.
.eeg_sampen <- function(x, m = 2L, r = 0.2, r_abs = NULL) {
  N <- length(x); if (N <= m + 1L) return(NA_real_)
  rr <- if (is.null(r_abs)) r * stats::sd(x) else r_abs
  if (!is.finite(rr) || rr <= 0) return(NA_real_)
  nv <- N - m
  B <- .eeg_c_matchcount(x, m, rr, nv)
  A <- .eeg_c_matchcount(x, m + 1L, rr, nv)
  if (A == 0 || B == 0) return(NA_real_)
  -log(A / B)
}

# Approximate entropy (Pincus 1991) -- includes self-matches, uses phi(m)-phi(m+1).
.eeg_apen <- function(x, m = 2L, r = 0.2) {
  N <- length(x); if (N <= m + 1L) return(NA_real_)
  rr <- r * stats::sd(x); if (!is.finite(rr) || rr <= 0) return(NA_real_)
  phi <- function(mm) {
    nv <- N - mm + 1L
    M <- vapply(seq_len(mm), function(j) x[j:(j + nv - 1L)], numeric(nv))
    Ci <- vapply(seq_len(nv), function(i) {
      dmax <- abs(M[, 1L] - M[i, 1L])
      if (mm > 1L) for (jj in 2:mm) dmax <- pmax(dmax, abs(M[, jj] - M[i, jj]))
      sum(dmax <= rr) / nv
    }, numeric(1))
    mean(log(Ci))
  }
  phi(m) - phi(m + 1L)
}

# Permutation entropy (Bandt & Pompe 2002), normalized to [0,1].
.eeg_permen <- function(x, m = 3L, tau = 1L, normalize = TRUE) {
  n <- length(x) - (m - 1L) * tau
  if (n <= 1L) return(NA_real_)
  emb <- vapply(seq_len(m), function(j) x[((j - 1L) * tau + 1L):((j - 1L) * tau + n)], numeric(n))
  patt <- apply(emb, 1L, function(v) paste(order(v), collapse = "-"))
  p <- table(patt) / n
  H <- -sum(p * log(p))
  if (normalize) H / log(factorial(m)) else H
}

# Multiscale entropy (Costa 2002): coarse-grain then sample-entropy, r fixed on
# the original series. Returns per-scale SampEn.
.eeg_mse <- function(x, m = 2L, r = 0.2, scales = 1:8) {
  rr <- r * stats::sd(x)
  vapply(scales, function(s) {
    cg <- if (s == 1L) x else {
      nb <- length(x) %/% s
      colMeans(matrix(x[seq_len(nb * s)], nrow = s))
    }
    .eeg_sampen(cg, m, r_abs = rr)
  }, numeric(1))
}

# Lempel-Ziv complexity (Kaspar & Schuster 1987), median-binarized, normalized.
.eeg_lziv <- function(x, normalize = TRUE) {
  s <- as.integer(x > stats::median(x))
  n <- length(s); if (n < 2L) return(NA_real_)
  i <- 1L; k <- 1L; l <- 2L; c <- 1L; kmax <- 1L
  repeat {
    if (s[i + k - 1L] == s[l + k - 1L]) {
      k <- k + 1L
      if (l + k - 1L > n) { c <- c + 1L; break }
    } else {
      if (k > kmax) kmax <- k
      i <- i + 1L
      if (i == l) {
        c <- c + 1L
        l <- l + kmax
        if (l > n) break
        i <- 1L; k <- 1L; kmax <- 1L
      } else {
        k <- 1L
      }
    }
  }
  if (normalize) c * log2(n) / n else c
}

# Higuchi fractal dimension (Higuchi 1988).
.eeg_higuchi <- function(x, kmax = NULL) {
  N <- length(x); if (N < 10L) return(NA_real_)
  # kmax is capped low (16) to avoid the aliasing pathology when a near-periodic
  # signal's period falls inside the k range. Normalization validated against
  # ground truth: white noise FD -> 2, Brownian motion -> 1.5.
  if (is.null(kmax)) kmax <- max(2L, min(floor(N / 10L), 16L))
  Lk <- vapply(seq_len(kmax), function(k) {
    Lm <- vapply(seq_len(k), function(mstart) {
      idx <- seq(mstart, N, by = k)
      ni <- length(idx) - 1L
      if (ni < 1L) return(NA_real_)
      sum(abs(diff(x[idx]))) * (N - 1) / (ni * k * k)
    }, numeric(1))
    mean(Lm, na.rm = TRUE)
  }, numeric(1))
  k <- seq_len(kmax); ok <- is.finite(Lk) & Lk > 0
  if (sum(ok) < 2L) return(NA_real_)
  fd <- unname(stats::coef(stats::lm(log(Lk[ok]) ~ log(1 / k[ok])))[2L])
  min(max(fd, 1), 2)                        # FD of a planar curve is bounded [1,2]
}

# Katz fractal dimension (Katz 1988).
.eeg_katz <- function(x) {
  n <- length(x) - 1L; if (n < 2L) return(NA_real_)
  L <- sum(sqrt(1 + diff(x)^2))
  d <- max(sqrt((seq_along(x) - 1L)^2 + (x - x[1L])^2))
  if (d <= 0 || L <= 0) return(NA_real_)
  log10(n) / (log10(n) + log10(d / L))
}

# Detrended fluctuation analysis (Peng 1994): scaling exponent alpha.
.eeg_dfa <- function(x, scales = NULL) {
  N <- length(x); if (N < 16L) return(NA_real_)
  y <- cumsum(x - mean(x))
  if (is.null(scales))
    scales <- unique(round(2^seq(log2(4), log2(N / 4), length.out = 16L)))
  scales <- scales[scales >= 4L & scales <= N %/% 2L]
  Fn <- vapply(scales, function(s) {
    nb <- N %/% s
    seg <- matrix(y[seq_len(nb * s)], nrow = s)
    tt <- seq_len(s)
    rms <- apply(seg, 2L, function(col) {
      res <- stats::lm.fit(cbind(1, tt), col)$residuals
      sqrt(mean(res^2))
    })
    sqrt(mean(rms^2))
  }, numeric(1))
  ok <- is.finite(Fn) & Fn > 0
  if (sum(ok) < 2L) return(NA_real_)
  unname(stats::coef(stats::lm(log(Fn[ok]) ~ log(scales[ok])))[2L])
}

# Hurst exponent by rescaled-range (R/S) analysis.
.eeg_hurst <- function(x) {
  N <- length(x); if (N < 16L) return(NA_real_)
  scales <- unique(round(2^seq(log2(8), log2(N / 2), length.out = 12L)))
  RS <- vapply(scales, function(s) {
    nb <- N %/% s; if (nb < 1L) return(NA_real_)
    vals <- vapply(seq_len(nb), function(b) {
      seg <- x[((b - 1L) * s + 1L):(b * s)]
      z <- cumsum(seg - mean(seg))
      S <- stats::sd(seg)
      if (S <= 0) NA_real_ else (max(z) - min(z)) / S
    }, numeric(1))
    mean(vals, na.rm = TRUE)
  }, numeric(1))
  ok <- is.finite(RS) & RS > 0
  if (sum(ok) < 2L) return(NA_real_)
  unname(stats::coef(stats::lm(log(RS[ok]) ~ log(scales[ok])))[2L])
}

# Hjorth parameters (Hjorth 1970): activity, mobility, complexity.
.eeg_hjorth <- function(x) {
  dx <- diff(x); ddx <- diff(dx)
  v0 <- stats::var(x); v1 <- stats::var(dx); v2 <- stats::var(ddx)
  mob <- if (v0 > 0) sqrt(v1 / v0) else NA_real_
  cmplx <- if (v1 > 0 && is.finite(mob) && mob > 0) sqrt(v2 / v1) / mob else NA_real_
  c(activity = v0, mobility = mob, complexity = cmplx)
}

# Spectral entropy (normalized Shannon entropy of the power spectrum).
.eeg_specentropy <- function(x, normalize = TRUE) {
  n <- length(x); if (n < 4L) return(NA_real_)
  p <- Mod(stats::fft(x - mean(x)))^2
  half <- 2:floor(n / 2)
  ps <- p[half]; s <- sum(ps); if (s <= 0) return(NA_real_)
  ps <- ps / s
  H <- -sum(ps * log(ps + 1e-12))
  if (normalize) H / log(length(ps)) else H
}

# per-channel vector list from a 2D (time x ch) or 3D (time x ch x epoch) assay
.eeg_channel_series <- function(data, ch) {
  if (length(dim(data)) == 3L) lapply(seq_len(dim(data)[3]), function(e) data[, ch, e])
  else list(data[, ch])
}

#' Singular Value Decomposition (SVD) Entropy
#'
#' Computes the SVD entropy of a one-dimensional signal (Roberts et al., 1999) --
#' a measure of the signal's dimensionality / the number of eigenvectors needed to
#' explain it. The signal is time-delay embedded into an
#' \code{order}-dimensional matrix; the normalized singular values of that matrix
#' are treated as a probability distribution, and their Shannon entropy is the SVD
#' entropy. Low SVD entropy indicates a low-dimensional (structured, e.g.
#' oscillatory) signal; high SVD entropy indicates a high-dimensional (complex or
#' noise-like) signal.
#'
#' @param x A numeric vector (the time series).
#' @param order Embedding dimension (default 3); the number of singular values.
#' @param delay Embedding delay / lag in samples (default 1).
#' @param normalize If \code{TRUE}, divide by \code{log2(order)} to bound the
#'   result in \eqn{[0, 1]}; otherwise return the entropy in bits (default
#'   \code{FALSE}).
#' @return A single numeric value, the SVD entropy (in bits, or normalized).
#'
#' @references
#' Roberts, S. J., Penny, W., & Rezek, I. (1999). "Temporal and spatial complexity
#'   measures for electroencephalogram based brain-computer interfacing."
#'   \emph{Medical & Biological Engineering & Computing}, 37(1), 93--98.
#'   \doi{10.1007/BF02513272}
#'
#' @seealso \code{\link{eegComplexity}} for the multi-measure complexity wrapper
#'   (permutation entropy, Lempel-Ziv, Hjorth, DFA, ...).
#'
#' @export
#' @examples
#' set.seed(1)
#' svdEntropy(sin(seq(0, 20 * pi, length.out = 500)))  # low (structured)
#' svdEntropy(rnorm(500))                               # high (noise-like)
svdEntropy <- function(x, order = 3L, delay = 1L, normalize = FALSE) {
  x <- as.numeric(x)
  order <- as.integer(order); delay <- as.integer(delay)
  n <- length(x)
  if (order < 2L) stop("`order` must be at least 2.", call. = FALSE)
  if (delay < 1L) stop("`delay` must be at least 1.", call. = FALSE)
  if (order * delay > n) stop("`order * delay` must be <= length(x).", call. = FALSE)
  m <- n - (order - 1L) * delay
  emb <- matrix(0, nrow = m, ncol = order)
  for (i in seq_len(order)) {
    start <- (i - 1L) * delay + 1L
    emb[, i] <- x[start:(start + m - 1L)]
  }
  sv <- svd(emb, nu = 0L, nv = 0L)$d          # singular values (LAPACK, as in numpy)
  sv <- sv / sum(sv)                          # normalized singular values
  h <- -sum(ifelse(sv > 0, sv * log2(sv), 0)) # Shannon entropy of the SV spectrum
  if (normalize) h <- h / log2(order)
  h
}


#' Petrosian Fractal Dimension
#'
#' Computes the Petrosian fractal dimension of a one-dimensional signal
#' (Petrosian, 1995) -- a fast waveform-complexity measure based on the number of
#' sign changes in the signal's derivative (i.e. the density of local extrema).
#' A smooth, oscillatory signal has few derivative sign changes and a Petrosian FD
#' near 1; a complex or noise-like signal has many and a higher FD.
#' \eqn{P = \log_{10} N / (\log_{10} N + \log_{10}(N / (N + 0.4\,N_\delta)))},
#' where \eqn{N} is the signal length and \eqn{N_\delta} the number of sign
#' changes in the first difference.
#'
#' @param x A numeric vector (the time series).
#' @return A single numeric value, the Petrosian fractal dimension (typically in
#'   \eqn{[1, \sim1.1]} for physiological signals).
#'
#' @references
#' Petrosian, A. (1995). "Kolmogorov complexity of finite sequences and
#'   recognition of different preictal EEG patterns." In \emph{Proceedings of the
#'   Eighth IEEE Symposium on Computer-Based Medical Systems}, 212--217.
#'   \doi{10.1109/CBMS.1995.465426}
#'
#' @seealso \code{\link{svdEntropy}} for the singular-value complexity view,
#'   \code{\link{eegComplexity}} for the multi-measure wrapper.
#'
#' @export
#' @examples
#' petrosianFD(sin(seq(0, 20 * pi, length.out = 500)))  # low (smooth oscillation)
#' set.seed(1); petrosianFD(rnorm(500))                 # higher (noise-like)
petrosianFD <- function(x) {
  x <- as.numeric(x)
  n <- length(x)
  if (n < 3L) stop("`x` must have length >= 3.", call. = FALSE)
  d <- diff(x)
  sb <- d < 0                              # sign of the derivative (as in np.signbit)
  nzc <- sum(sb[-1] != sb[-length(sb)])    # sign changes in the derivative
  log10(n) / (log10(n) + log10(n / (n + 0.4 * nzc)))
}


#' Renyi Entropy
#'
#' Computes the Renyi entropy of order \code{alpha} of a one-dimensional signal
#' (Renyi, 1961) -- a one-parameter generalization of the Shannon entropy of the
#' signal's amplitude distribution. The signal is discretized into \code{bins}
#' equal-width bins to form a probability distribution \eqn{p}, and
#' \eqn{H_\alpha = \frac{1}{1-\alpha}\ln\left(\sum_i p_i^\alpha\right)} is returned
#' in nats. Special cases: \eqn{\alpha \to 1} is the Shannon entropy,
#' \eqn{\alpha = 2} the collision entropy, and \eqn{\alpha \to \infty} the
#' min-entropy. \eqn{H_\alpha} is non-increasing in \eqn{\alpha}.
#'
#' @param x A numeric vector (the time series).
#' @param alpha Order of the Renyi entropy (default 2, the collision entropy);
#'   \code{alpha = 1} returns the Shannon entropy.
#' @param bins Number of equal-width bins for the amplitude histogram (default 16).
#' @return A single numeric value, the Renyi entropy in nats.
#'
#' @references
#' Renyi, A. (1961). "On measures of entropy and information." In \emph{Proceedings
#'   of the Fourth Berkeley Symposium on Mathematical Statistics and Probability},
#'   Volume 1, 547--561.
#'
#' @seealso \code{\link{svdEntropy}}, \code{\link{petrosianFD}},
#'   \code{\link{eegComplexity}} for other complexity measures.
#'
#' @export
#' @examples
#' set.seed(1)
#' renyiEntropy(rnorm(1000), alpha = 2)   # collision entropy
#' renyiEntropy(rnorm(1000), alpha = 1)   # = Shannon entropy
renyiEntropy <- function(x, alpha = 2, bins = 16L) {
  x <- as.numeric(x); bins <- as.integer(bins)
  if (is.na(bins) || bins < 2L) stop("`bins` must be an integer >= 2.", call. = FALSE)
  if (!is.numeric(alpha) || length(alpha) != 1L || alpha < 0) {
    stop("`alpha` must be a single non-negative number.", call. = FALSE)
  }
  b <- cut(x, breaks = seq(min(x), max(x), length.out = bins + 1L),
           include.lowest = TRUE, labels = FALSE)
  p <- as.numeric(table(factor(b, levels = seq_len(bins))))
  p <- p / sum(p); p <- p[p > 0]
  if (abs(alpha - 1) < 1e-12) {
    -sum(p * log(p))                        # Shannon limit (alpha -> 1), in nats
  } else {
    (1 / (1 - alpha)) * log(sum(p^alpha))   # Renyi entropy, in nats
  }
}


#' Tsallis Entropy
#'
#' Computes the Tsallis entropy of order \code{q} of a one-dimensional signal
#' (Tsallis, 1988) -- the \emph{non-additive} one-parameter generalization of the
#' Shannon entropy of the signal's amplitude distribution, complementary to the
#' (additive) \code{\link{renyiEntropy}}. The signal is discretized into
#' \code{bins} equal-width bins to form a probability distribution \eqn{p}, and
#' \eqn{S_q = \frac{1 - \sum_i p_i^q}{q - 1}} is returned (equivalently
#' \eqn{\sum_i p_i \ln_q(1/p_i)} with the q-logarithm). \eqn{q \to 1} recovers the
#' Shannon entropy (in nats); unlike Renyi, Tsallis is non-additive
#' (\eqn{S_q(A,B) = S_q(A) + S_q(B) + (1-q) S_q(A) S_q(B)} for independent
#' \eqn{A, B}), so Renyi and Tsallis agree only at \eqn{q = 1}.
#'
#' @param x A numeric vector (the time series).
#' @param q Order of the Tsallis entropy (default 2); \code{q = 1} returns the
#'   Shannon entropy.
#' @param bins Number of equal-width bins for the amplitude histogram (default 16).
#' @return A single numeric value, the Tsallis entropy.
#'
#' @references
#' Tsallis, C. (1988). "Possible generalization of Boltzmann-Gibbs statistics."
#'   \emph{Journal of Statistical Physics}, 52(1-2), 479--487.
#'   \doi{10.1007/BF01016429}
#'
#' @seealso \code{\link{renyiEntropy}} for the additive generalization,
#'   \code{\link{svdEntropy}}, \code{\link{petrosianFD}}.
#'
#' @export
#' @examples
#' set.seed(1)
#' tsallisEntropy(rnorm(1000), q = 2)   # non-additive (collision-like)
#' tsallisEntropy(rnorm(1000), q = 1)   # = Shannon entropy
tsallisEntropy <- function(x, q = 2, bins = 16L) {
  x <- as.numeric(x); bins <- as.integer(bins)
  if (is.na(bins) || bins < 2L) stop("`bins` must be an integer >= 2.", call. = FALSE)
  if (!is.numeric(q) || length(q) != 1L || q < 0) {
    stop("`q` must be a single non-negative number.", call. = FALSE)
  }
  b <- cut(x, breaks = seq(min(x), max(x), length.out = bins + 1L),
           include.lowest = TRUE, labels = FALSE)
  p <- as.numeric(table(factor(b, levels = seq_len(bins))))
  p <- p / sum(p); p <- p[p > 0]
  lnq <- if (abs(q - 1) < 1e-12) log(1 / p) else ((1 / p)^(1 - q) - 1) / (1 - q)  # q-logarithm
  sum(p * lnq)                                    # = (1 - sum p^q)/(q - 1); Shannon at q -> 1
}


#' Dispersion Entropy
#'
#' Computes the dispersion entropy (DispEn) of a one-dimensional signal
#' (Rostaghi & Azami, 2016) -- a fast, robust symbolic-dynamics complexity
#' measure. Unlike the amplitude-histogram entropies (\code{\link{renyiEntropy}},
#' \code{\link{tsallisEntropy}}) and the ordinal permutation entropy, DispEn maps
#' the signal to \code{c} amplitude classes through the normal cumulative
#' distribution function (NCDF), forms embedding vectors of length
#' \code{dimension} (the \emph{dispersion patterns}), and takes the Shannon
#' entropy of the pattern distribution, normalized by \eqn{\ln(c^{m})}. The NCDF
#' mapping makes it insensitive to outliers and much cheaper than sample entropy,
#' with no threshold-crossing boundary effects. Set \code{reverse = TRUE} to
#' return the reverse dispersion entropy (RDEn), the mean-square deviation of the
#' pattern distribution from uniformity (larger = more regular / less complex).
#'
#' This implementation reproduces \code{NeuroKit2}'s \code{entropy_dispersion}
#' bit-for-bit, including its normalization convention (a base-2 Shannon term
#' divided by the natural-log \eqn{\ln(c^{m})}, so DispEn ranges in
#' \eqn{[0, 1/\ln 2]}).
#'
#' @param x A numeric vector (the time series).
#' @param c Number of amplitude classes / symbols (default 6; Rostaghi & Azami
#'   recommend 4--8).
#' @param dimension Embedding dimension \eqn{m} of the dispersion patterns
#'   (default 3).
#' @param delay Time delay (lag) in samples (default 1).
#' @param reverse If \code{TRUE}, return the reverse dispersion entropy (RDEn)
#'   instead of DispEn (default \code{FALSE}).
#' @return A single numeric value: the dispersion entropy (or, if
#'   \code{reverse = TRUE}, the reverse dispersion entropy).
#'
#' @references
#' Rostaghi, M., & Azami, H. (2016). "Dispersion entropy: A measure for
#'   time-series analysis." \emph{IEEE Signal Processing Letters}, 23(5),
#'   610--614. \doi{10.1109/LSP.2016.2542881}
#'
#' @seealso \code{\link{renyiEntropy}}, \code{\link{tsallisEntropy}},
#'   \code{\link{svdEntropy}}, \code{\link{petrosianFD}}.
#'
#' @export
#' @examples
#' set.seed(1)
#' dispersionEntropy(rnorm(1000))                 # DispEn (c = 6, m = 3)
#' dispersionEntropy(rnorm(1000), reverse = TRUE) # reverse dispersion entropy
dispersionEntropy <- function(x, c = 6L, dimension = 3L, delay = 1L,
                              reverse = FALSE) {
  x <- as.numeric(x)
  c <- as.integer(c); dimension <- as.integer(dimension); delay <- as.integer(delay)
  if (is.na(c) || c < 2L) stop("`c` must be an integer >= 2.", call. = FALSE)
  if (is.na(dimension) || dimension < 1L) stop("`dimension` must be an integer >= 1.", call. = FALSE)
  if (is.na(delay) || delay < 1L) stop("`delay` must be an integer >= 1.", call. = FALSE)
  s <- stats::sd(x)
  if (is.na(s) || s == 0) stop("`x` has zero variance; dispersion entropy is undefined.", call. = FALSE)
  # NCDF symbolization: standardize (sample sd) -> normal CDF -> c equal-probability classes
  y <- stats::pnorm((x - mean(x)) / s)
  sym <- findInterval(y, (0:(c - 1L)) / c)         # dispatch to classes 1..c
  n_rows <- length(sym) - (dimension - 1L) * delay
  if (n_rows < 1L) {
    stop("Signal too short for the requested `dimension` and `delay`.", call. = FALSE)
  }
  # time-delay embedding -> dispersion patterns (rows)
  emb <- vapply(seq_len(dimension),
                function(k) sym[((k - 1L) * delay + 1L):((k - 1L) * delay + n_rows)],
                integer(n_rows))
  key <- do.call(paste, c(as.data.frame(emb), sep = ","))
  p <- as.numeric(table(key)); p <- p / sum(p)
  n_possible <- c^dimension
  if (reverse) {
    unname(sum((p - 1 / n_possible)^2) / (1 - 1 / n_possible))   # RDEn (deviation from uniform)
  } else {
    unname((-sum(p * log2(p))) / log(n_possible))                # DispEn (NeuroKit2 normalization)
  }
}


#' Fuzzy Entropy
#'
#' Computes the fuzzy entropy (FuzzyEn) of a one-dimensional signal (Chen et al.,
#' 2007) -- the \emph{fuzzy} generalization of the sample entropy. Sample entropy
#' counts template matches with a crisp threshold (a Heaviside step: two vectors
#' either match within tolerance \eqn{r} or do not); fuzzy entropy replaces that
#' hard step with a smooth exponential membership
#' \eqn{\mu = \exp(-(d^{n})/r)}, so near-matches contribute partially. This makes
#' FuzzyEn markedly more robust and continuous than sample entropy, especially on
#' short or noisy segments where the crisp count is unstable. Each embedding
#' vector has its own baseline (mean) removed before the Chebyshev distance is
#' taken, following Chen et al.
#'
#' This implementation reproduces \code{NeuroKit2}'s \code{entropy_fuzzy}
#' bit-for-bit (the default \code{tolerance = "sd"} there corresponds to
#' \code{r = 0.2}).
#'
#' @param x A numeric vector (the time series).
#' @param dimension Embedding dimension \eqn{m} (default 2).
#' @param r Tolerance as a fraction of the signal's standard deviation; the
#'   distance tolerance is \code{r * sd(x)} (default 0.2, matching NeuroKit2's
#'   \code{tolerance = "sd"}).
#' @param delay Time delay (lag) in samples (default 1).
#' @param n Fuzzy power (exponent) of the membership function (default 1).
#' @return A single numeric value, the fuzzy entropy.
#'
#' @references
#' Chen, W., Zhuang, J., Yu, W., & Wang, Z. (2009). "Measuring complexity using
#'   FuzzyEn, ApEn, and SampEn." \emph{Medical Engineering & Physics}, 31(1),
#'   61--68. \doi{10.1016/j.medengphy.2008.04.005}
#'
#' @seealso \code{\link{dispersionEntropy}}, \code{\link{svdEntropy}},
#'   \code{\link{renyiEntropy}}, \code{\link{tsallisEntropy}}.
#'
#' @export
#' @examples
#' set.seed(1)
#' fuzzyEntropy(rnorm(300), dimension = 2)   # fuzzy generalization of sample entropy
fuzzyEntropy <- function(x, dimension = 2L, r = 0.2, delay = 1L, n = 1) {
  x <- as.numeric(x)
  dimension <- as.integer(dimension); delay <- as.integer(delay)
  if (is.na(dimension) || dimension < 1L) stop("`dimension` must be an integer >= 1.", call. = FALSE)
  if (is.na(delay) || delay < 1L) stop("`delay` must be an integer >= 1.", call. = FALSE)
  if (!is.numeric(r) || length(r) != 1L || r <= 0) stop("`r` must be a single positive number.", call. = FALSE)
  if (!is.numeric(n) || length(n) != 1L || n <= 0) stop("`n` must be a single positive number.", call. = FALSE)
  s <- stats::sd(x)
  if (is.na(s) || s == 0) stop("`x` has zero variance; fuzzy entropy is undefined.", call. = FALSE)
  r_tol <- r * s
  N <- length(x)
  embed <- function(d) {
    n_rows <- N - (d - 1L) * delay
    if (n_rows < 2L) stop("Signal too short for the requested `dimension` and `delay`.", call. = FALSE)
    vapply(seq_len(d), function(k) x[((k - 1L) * delay + 1L):((k - 1L) * delay + n_rows)], numeric(n_rows))
  }
  # phi at embedding dimension d; `drop_last` mirrors NeuroKit2's approximate=FALSE (m) vs TRUE (m+1)
  phi <- function(d, drop_last) {
    e <- embed(d)
    if (drop_last) e <- e[-nrow(e), , drop = FALSE]
    e <- e - rowMeans(e)                                   # remove each vector's local baseline (fuzzy)
    dmat <- as.matrix(stats::dist(e, method = "maximum"))  # Chebyshev distance
    sim <- exp(-(dmat^n) / r_tol)                          # fuzzy membership (crisp step -> smooth)
    count <- colSums(sim)                                  # includes the self-match (mu = 1)
    mean((count - 1) / (nrow(e) - 1))                      # remove self-match, as in sample entropy
  }
  -log(phi(dimension + 1L, FALSE) / phi(dimension, TRUE))
}


#' Increment Entropy
#'
#' Computes the increment entropy (IncrEn) of a one-dimensional signal (Liu et
#' al., 2016) -- a symbolic complexity measure defined on the signal's
#' \emph{increments} (successive differences) rather than its amplitudes. Each
#' increment in an embedding vector is mapped to a word letter combining its
#' \emph{sign} (\eqn{-, 0, +}) with a \emph{magnitude class} in
#' \eqn{\{0, 1, \dots, q\}} obtained by quantizing \eqn{|{\rm increment}|}
#' relative to the vector's own standard deviation; IncrEn is the Shannon entropy
#' of the resulting word distribution, normalized by \eqn{m - 1}. Because it
#' encodes both the direction and the graded size of changes, it is sensitive to
#' dynamics that amplitude-based entropies miss, and is a distinct mechanism from
#' the amplitude-histogram (\code{\link{renyiEntropy}}, \code{\link{tsallisEntropy}}),
#' NCDF-pattern (\code{\link{dispersionEntropy}}) and template-matching
#' (\code{\link{fuzzyEntropy}}) entropies.
#'
#' This implementation reproduces \code{NeuroKit2}'s \code{entropy_increment}
#' bit-for-bit.
#'
#' @param x A numeric vector (the time series).
#' @param dimension Embedding dimension \eqn{m} of the increment words
#'   (default 2).
#' @param q Number of magnitude quantization levels (default 4); each increment
#'   is binned into one of \code{q + 1} magnitude classes by its size relative to
#'   the embedding vector's standard deviation.
#' @return A single numeric value, the increment entropy.
#'
#' @references
#' Liu, X., Jiang, A., Xu, N., & Xue, J. (2016). "Increment Entropy as a Measure
#'   of Complexity for Time Series." \emph{Entropy}, 18(1), 22.
#'   \doi{10.3390/e18010022}
#'
#' @seealso \code{\link{dispersionEntropy}}, \code{\link{fuzzyEntropy}},
#'   \code{\link{renyiEntropy}}, \code{\link{tsallisEntropy}}.
#'
#' @export
#' @examples
#' set.seed(1)
#' incrementEntropy(rnorm(1000), dimension = 2, q = 4)
incrementEntropy <- function(x, dimension = 2L, q = 4L) {
  x <- as.numeric(x)
  dimension <- as.integer(dimension); q <- as.integer(q)
  if (is.na(dimension) || dimension < 2L) stop("`dimension` must be an integer >= 2.", call. = FALSE)
  if (is.na(q) || q < 1L) stop("`q` must be an integer >= 1.", call. = FALSE)
  dx <- diff(x)
  n_rows <- length(dx) - (dimension - 1L)
  if (n_rows < 2L) stop("Signal too short for the requested `dimension`.", call. = FALSE)
  # embed the increments (delay 1)
  emb <- vapply(seq_len(dimension),
                function(k) dx[k:(k + n_rows - 1L)], numeric(n_rows))
  sgn <- sign(emb)
  row_sd <- apply(emb, 1L, stats::sd)                       # per-vector sd (ddof = 1)
  tempm <- matrix(row_sd, n_rows, dimension)
  ratio <- abs(emb) * q / tempm
  ratio[!is.finite(ratio)] <- 0                            # sd == 0 rows -> reset below
  size <- floor(ratio); size[size > q] <- q                # magnitude class in 0..q
  size[row_sd == 0, ] <- 0
  words <- sgn * size                                      # letters: sign x magnitude class
  key <- do.call(paste, c(as.data.frame(words), sep = ","))
  p <- as.numeric(table(key)); p <- p / sum(p)
  (-sum(p * log2(p))) / (dimension - 1L)                   # normalized Shannon (bits)
}


#' Slope Entropy
#'
#' Computes the slope entropy (SlopEn) of a one-dimensional signal (Cuesta-Frau,
#' 2019) -- a symbolic complexity measure defined on the \emph{slopes} between
#' successive samples. Each first difference is turned into the \emph{angle}
#' \eqn{\theta = \arctan(\Delta x)} (in degrees) and mapped to one of five
#' symbols by two angle thresholds: steep-up (\eqn{> \gamma}), gentle-up
#' (\eqn{(\delta, \gamma]}), flat (\eqn{[-\delta, \delta]}), gentle-down and
#' steep-down; SlopEn is the Shannon entropy of the resulting length-(m-1) slope
#' patterns. By combining a coarse (steep/gentle) and fine (flat) angular
#' resolution it captures the shape of local trends, a mechanism distinct from
#' the amplitude, NCDF-pattern, template-matching and increment-magnitude
#' entropies.
#'
#' This implementation reproduces \code{NeuroKit2}'s \code{entropy_slope}
#' bit-for-bit.
#'
#' @param x A numeric vector (the time series).
#' @param dimension Embedding dimension \eqn{m} (default 3); slope patterns have
#'   length \code{dimension - 1}.
#' @param thresholds Increasing positive angle thresholds (degrees) that split
#'   the slope into symbol classes (default \code{c(0.1, 45)}: a flat band within
#'   +/-0.1 deg, a gentle band up to +/-45 deg, and a steep class beyond).
#' @param delay Time delay (lag) in samples for the difference (default 1).
#' @return A single numeric value, the slope entropy (bits).
#'
#' @references
#' Cuesta-Frau, D. (2019). "Slope Entropy: A New Time Series Complexity Estimator
#'   Based on Both Symbolic Patterns and Amplitude Information." \emph{Entropy},
#'   21(12), 1167. \doi{10.3390/e21121167}
#'
#' @seealso \code{\link{incrementEntropy}}, \code{\link{dispersionEntropy}},
#'   \code{\link{fuzzyEntropy}}.
#'
#' @export
#' @examples
#' set.seed(1)
#' slopeEntropy(rnorm(1000), dimension = 3)
slopeEntropy <- function(x, dimension = 3L, thresholds = c(0.1, 45), delay = 1L) {
  x <- as.numeric(x)
  dimension <- as.integer(dimension); delay <- as.integer(delay)
  th <- as.numeric(thresholds)
  if (is.na(dimension) || dimension < 2L) stop("`dimension` must be an integer >= 2.", call. = FALSE)
  if (is.na(delay) || delay < 1L) stop("`delay` must be an integer >= 1.", call. = FALSE)
  if (length(th) < 2L || any(th <= 0) || is.unsorted(th)) {
    stop("`thresholds` must be >= 2 increasing positive values.", call. = FALSE)
  }
  n <- length(x)
  if (n <= delay + dimension) stop("Signal too short for the requested `dimension` and `delay`.", call. = FALSE)
  tx <- atan(x[(delay + 1L):n] - x[1:(n - delay)]) * 180 / pi   # slope angle in degrees
  ntx <- length(tx)
  symbols <- numeric(ntx)
  for (q in 2:length(th)) {                                     # thresholds -> {-k..k} symbols
    symbols[tx <= th[q] & tx > th[q - 1L]] <- q - 1L
    symbols[tx >= -th[q] & tx < -th[q - 1L]] <- -(q - 1L)
    if (q == length(th)) {
      symbols[tx > th[q]]  <- q
      symbols[tx < -th[q]] <- -q
    }
  }
  n_pat <- ntx - dimension + 1L                                 # count of length-(m-1) slope patterns
  emb <- vapply(seq_len(dimension - 1L),
                function(k) symbols[k:(k + n_pat - 1L)], numeric(n_pat))
  if (is.null(dim(emb))) emb <- matrix(emb, ncol = 1L)
  key <- do.call(paste, c(as.data.frame(emb), sep = ","))
  p <- as.numeric(table(key)); p <- p / sum(p)
  -sum(p * log2(p))                                             # Shannon entropy (bits)
}


#' Phase Entropy
#'
#' Computes the phase entropy (PhasEn) of a one-dimensional signal (Rohila &
#' Sharma, 2019) from its Second-Order Difference Plot (SODP). Each point plots
#' the first difference \eqn{x(n+\tau)-x(n)} against the second difference
#' \eqn{x(n+2\tau)-x(n+\tau)}; the polar angle of each point is binned into
#' \eqn{k} equal sectors of \eqn{[0, 2\pi)}, and PhasEn is the (normalized)
#' Shannon entropy of the angle-weighted sector distribution. It summarizes how
#' the signal's local rate-of-change and its change-of-rate co-vary -- a
#' phase-plane geometry distinct from the amplitude, template-matching and
#' symbolic (dispersion/increment/slope) entropies.
#'
#' This implementation reproduces \code{NeuroKit2}'s \code{entropy_phase}
#' bit-for-bit, including its normalization (a base-2 Shannon term divided by the
#' natural-log \eqn{\ln k}, so PhasEn ranges in \eqn{[0, 1/\ln 2]}).
#'
#' @param x A numeric vector (the time series).
#' @param delay Time delay (lag) in samples for the differences (default 1).
#' @param k Number of angular sectors of the SODP (default 4).
#' @return A single numeric value, the phase entropy.
#'
#' @references
#' Rohila, A., & Sharma, A. (2019). "Phase entropy: a new complexity measure for
#'   heart rate variability." \emph{Physiological Measurement}, 40(10), 105006.
#'   \doi{10.1088/1361-6579/ab499e}
#'
#' @seealso \code{\link{slopeEntropy}}, \code{\link{incrementEntropy}},
#'   \code{\link{dispersionEntropy}}.
#'
#' @export
#' @examples
#' set.seed(1)
#' phaseEntropy(rnorm(1000), k = 4)
phaseEntropy <- function(x, delay = 1L, k = 4L) {
  x <- as.numeric(x)
  delay <- as.integer(delay); k <- as.integer(k)
  if (is.na(delay) || delay < 1L) stop("`delay` must be an integer >= 1.", call. = FALSE)
  if (is.na(k) || k < 2L) stop("`k` must be an integer >= 2.", call. = FALSE)
  n <- length(x)
  if (n <= 2L * delay) stop("Signal too short for the requested `delay`.", call. = FALSE)
  yy <- x[(2L * delay + 1L):n] - x[(delay + 1L):(n - delay)]   # second difference
  xx <- x[(delay + 1L):(n - delay)] - x[1:(n - 2L * delay)]    # first difference
  theta <- atan(yy / xx)                                       # SODP polar angle (quadrant-corrected)
  theta[yy < 0 & xx < 0] <- theta[yy < 0 & xx < 0] + pi
  theta[yy < 0 & xx > 0] <- theta[yy < 0 & xx > 0] + 2 * pi
  theta[yy > 0 & xx < 0] <- theta[yy > 0 & xx < 0] + pi
  angles <- seq(0, 2 * pi, length.out = k + 1L)
  freq <- vapply(seq_len(k),
                 function(i) sum(theta[theta > angles[i] & theta < angles[i + 1L]]),
                 numeric(1))                                    # angle-weighted sector mass
  freq <- freq / sum(freq)
  p <- freq[freq > 0]
  (-sum(p * log2(p))) / log(k)                                 # NeuroKit2 normalization (base-2 / ln k)
}


#' EEG complexity and nonlinear-dynamics measures
#'
#' Computes per-channel entropy, fractal-dimension and detrended-fluctuation
#' complexity measures — the nonlinear family that complements the spectral,
#' connectivity and source tools. For epoched (3-D) data each measure is computed
#' per epoch and averaged across epochs.
#'
#' Available `measures`: `"sample_entropy"` (Richman–Moorman), `"approximate_entropy"`
#' (Pincus), `"permutation_entropy"` (Bandt–Pompe), `"multiscale_entropy"` (Costa;
#' returns the scale-averaged SampEn as a summary plus the per-scale curve in the
#' attribute), `"lempel_ziv"` (Kaspar–Schuster), `"higuchi_fd"`, `"katz_fd"`,
#' `"dfa"` (scaling exponent alpha), `"hurst"` (R/S), `"hjorth_mobility"`,
#' `"hjorth_complexity"`, `"spectral_entropy"`.
#'
#' @param pe A `PhysioExperiment`.
#' @param measures Character vector of measures to compute (see Details).
#' @param assay_name Assay to use (default: the object's default assay).
#' @param m Embedding dimension for the entropy measures (default 2; permutation
#'   entropy uses `perm_order`).
#' @param r Tolerance as a fraction of each channel's SD for sample/approximate/
#'   multiscale entropy (default 0.2).
#' @param tau Time delay for permutation entropy (default 1).
#' @param perm_order Embedding order for permutation entropy (default 3).
#' @param mse_scales Scales for multiscale entropy (default `1:8`).
#' @param max_samples Cap on samples per channel/epoch for the O(N^2) entropy
#'   measures; longer series are truncated with a warning (default 4000; `NULL`
#'   disables the cap).
#' @return A data frame with one row per channel and one column per requested
#'   measure (plus a `channel` column). The multiscale-entropy per-scale curves
#'   are attached as `attr(, "mse_scales")` / `attr(, "mse_curve")`.
#' @references Richman & Moorman 2000; Bandt & Pompe 2002; Costa 2002; Higuchi
#'   1988; Katz 1988; Peng 1994; Hjorth 1970.
#' @seealso [eegQEEG()], [eegAperiodic()]
#' @export
#' @examples
#' pe <- make_eeg(n_time = 512, n_channels = 4, sr = 128)
#' cx <- eegComplexity(pe, measures = c("permutation_entropy", "hjorth_mobility",
#'                                      "spectral_entropy"))
#' cx
eegComplexity <- function(pe,
                          measures = c("sample_entropy", "permutation_entropy",
                                       "lempel_ziv", "higuchi_fd", "dfa",
                                       "hjorth_mobility", "spectral_entropy"),
                          assay_name = NULL, m = 2L, r = 0.2, tau = 1L,
                          perm_order = 3L, mse_scales = 1:8, max_samples = 4000L) {
  stopifnot(inherits(pe, "PhysioExperiment"))
  known <- c("sample_entropy", "approximate_entropy", "permutation_entropy",
             "multiscale_entropy", "lempel_ziv", "higuchi_fd", "katz_fd", "dfa",
             "hurst", "hjorth_mobility", "hjorth_complexity", "spectral_entropy")
  bad <- setdiff(measures, known)
  if (length(bad)) stop("Unknown measure(s): ", paste(bad, collapse = ", "),
                        ". Available: ", paste(known, collapse = ", "), call. = FALSE)
  if (is.null(assay_name)) assay_name <- defaultAssay(pe)
  data <- SummarizedExperiment::assay(pe, assay_name)
  cd <- SummarizedExperiment::colData(pe)
  n_ch <- if (length(dim(data)) >= 2L) dim(data)[2] else 1L
  ch_labels <- if ("label" %in% colnames(cd)) as.character(cd$label) else
    paste0("ch", seq_len(n_ch))

  entropy_measures <- c("sample_entropy", "approximate_entropy", "multiscale_entropy")
  truncated <- FALSE
  cap <- function(v) {
    if (!is.null(max_samples) && length(v) > max_samples) { truncated <<- TRUE; v[seq_len(max_samples)] }
    else v
  }

  one_series <- function(v) {
    out <- list()
    for (meas in measures) {
      vv <- if (meas %in% entropy_measures) cap(v) else v
      out[[meas]] <- switch(meas,
        sample_entropy      = .eeg_sampen(vv, m, r),
        approximate_entropy = .eeg_apen(vv, m, r),
        permutation_entropy = .eeg_permen(v, perm_order, tau),
        multiscale_entropy  = mean(.eeg_mse(vv, m, r, mse_scales), na.rm = TRUE),
        lempel_ziv          = .eeg_lziv(v),
        higuchi_fd          = .eeg_higuchi(v),
        katz_fd             = .eeg_katz(v),
        dfa                 = .eeg_dfa(v),
        hurst               = .eeg_hurst(v),
        hjorth_mobility     = .eeg_hjorth(v)[["mobility"]],
        hjorth_complexity   = .eeg_hjorth(v)[["complexity"]],
        spectral_entropy    = .eeg_specentropy(v))
    }
    out
  }

  mse_curves <- vector("list", n_ch)
  rows <- lapply(seq_len(n_ch), function(ch) {
    series <- .eeg_channel_series(data, ch)
    per <- lapply(series, one_series)
    agg <- lapply(measures, function(meas)
      mean(vapply(per, function(p) p[[meas]], numeric(1)), na.rm = TRUE))
    names(agg) <- measures
    if ("multiscale_entropy" %in% measures) {
      curves <- vapply(series, function(v) .eeg_mse(cap(v), m, r, mse_scales),
                       numeric(length(mse_scales)))
      mse_curves[[ch]] <<- rowMeans(matrix(curves, nrow = length(mse_scales)), na.rm = TRUE)
    }
    as.data.frame(agg, stringsAsFactors = FALSE)
  })
  df <- cbind(channel = ch_labels, do.call(rbind, rows), stringsAsFactors = FALSE)
  if (truncated)
    warning("Some channels exceeded max_samples (", max_samples,
            ") and were truncated for the O(N^2) entropy measures.", call. = FALSE)
  if ("multiscale_entropy" %in% measures) {
    attr(df, "mse_scales") <- mse_scales
    attr(df, "mse_curve") <- stats::setNames(mse_curves, ch_labels)
  }
  rownames(df) <- NULL
  df
}

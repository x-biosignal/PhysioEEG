library(testthat)
library(PhysioEEG)

mk_pe <- function(mat, sr = 200) {
  labs <- colnames(mat)
  dimnames(mat) <- NULL
  pe <- make_eeg(n_time = nrow(mat), n_channels = ncol(mat), sr = sr)
  SummarizedExperiment::assay(pe, defaultAssay(pe)) <- mat
  if (!is.null(labs)) SummarizedExperiment::colData(pe)$label <- labs
  pe
}

test_that("fractal / scaling cores recover known values on synthetic signals", {
  set.seed(1); N <- 4000
  wn <- rnorm(N); brown <- cumsum(rnorm(N))
  # Higuchi: white noise -> ~2, Brownian motion -> ~1.5
  expect_equal(PhysioEEG:::.eeg_higuchi(wn), 2, tolerance = 0.12)
  expect_equal(PhysioEEG:::.eeg_higuchi(brown), 1.5, tolerance = 0.18)
  expect_true(PhysioEEG:::.eeg_higuchi(wn) <= 2 && PhysioEEG:::.eeg_higuchi(wn) >= 1)
  # DFA / Hurst of white noise -> ~0.5
  expect_equal(PhysioEEG:::.eeg_dfa(wn), 0.5, tolerance = 0.12)
  expect_equal(PhysioEEG:::.eeg_hurst(wn), 0.5, tolerance = 0.15)
  # correlated AR(0.9) -> DFA alpha clearly above 0.5
  ar <- as.numeric(stats::filter(rnorm(N), 0.9, method = "recursive"))
  expect_gt(PhysioEEG:::.eeg_dfa(ar), 0.7)
})

test_that("entropy cores behave (white high, sine low; permutation normalized)", {
  set.seed(2); N <- 2000; t <- (0:(N - 1)) / 200
  wn <- rnorm(N); sine <- sin(2 * pi * 10 * t)
  expect_gt(PhysioEEG:::.eeg_sampen(wn), PhysioEEG:::.eeg_sampen(sine))
  expect_gt(PhysioEEG:::.eeg_permen(wn), PhysioEEG:::.eeg_permen(sine))
  expect_lte(PhysioEEG:::.eeg_permen(wn), 1)          # normalized to [0,1]
  expect_gt(PhysioEEG:::.eeg_lziv(wn), PhysioEEG:::.eeg_lziv(sine))
  expect_lt(PhysioEEG:::.eeg_specentropy(sine), 0.2)  # single tone -> low
  expect_gt(PhysioEEG:::.eeg_specentropy(wn), 0.8)    # broadband -> high
})

test_that("eegComplexity returns per-channel table and discriminates signals", {
  set.seed(3); N <- 2000; sr <- 200; t <- (0:(N - 1)) / sr
  pe <- mk_pe(cbind(white = rnorm(N), sine = sin(2 * pi * 10 * t)), sr)
  cx <- eegComplexity(pe, measures = c("spectral_entropy", "lempel_ziv",
        "higuchi_fd", "permutation_entropy", "hjorth_mobility", "sample_entropy"))
  expect_s3_class(cx, "data.frame")
  expect_equal(nrow(cx), 2L)
  expect_setequal(cx$channel, c("white", "sine"))
  w <- cx[cx$channel == "white", ]; s <- cx[cx$channel == "sine", ]
  expect_gt(w$spectral_entropy, s$spectral_entropy)
  expect_gt(w$lempel_ziv, s$lempel_ziv)
  expect_gt(w$higuchi_fd, s$higuchi_fd)
  expect_gt(w$permutation_entropy, s$permutation_entropy)
  expect_true(all(cx$higuchi_fd >= 1 & cx$higuchi_fd <= 2))
})

test_that("eegComplexity attaches multiscale-entropy curves", {
  set.seed(4)
  pe <- mk_pe(matrix(rnorm(1500 * 2), ncol = 2), sr = 200)
  cx <- eegComplexity(pe, measures = c("multiscale_entropy", "dfa"), mse_scales = 1:5)
  expect_true(all(c("multiscale_entropy", "dfa") %in% names(cx)))
  expect_equal(length(attr(cx, "mse_scales")), 5L)
  expect_equal(length(attr(cx, "mse_curve")), 2L)
  expect_equal(length(attr(cx, "mse_curve")[[1]]), 5L)
})

test_that("eegComplexity errors on unknown measure and warns on truncation", {
  pe <- mk_pe(matrix(rnorm(500), ncol = 1), sr = 100)
  expect_error(eegComplexity(pe, measures = "nope"), "Unknown measure")
  peL <- mk_pe(matrix(rnorm(2500), ncol = 1), sr = 200)
  expect_warning(eegComplexity(peL, measures = "sample_entropy", max_samples = 1000),
                 "truncated")
})

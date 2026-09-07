library(testthat)
library(PhysioEEG)

# per-source power (sum over the 3 orientation columns) of a source estimate.
.srcpow <- function(S, ns) {
  vapply(seq_len(ns), function(j)
    sum(colMeans(S[, ((j - 1) * 3 + 1):(j * 3)]^2)), numeric(1))
}

test_that("realistic forward model returns a valid leadfield with source normals", {
  skip_if_not_installed("PhysioHeadModels")
  pe <- make_eeg(n_time = 100, n_channels = 19, sr = 250)
  fm <- eegForwardModel(pe, method = "sphere_analytic", n_sources = 150)
  expect_equal(nrow(fm$leadfield), 19L)
  expect_equal(ncol(fm$leadfield), fm$n_sources * 3L)
  expect_false(is.null(fm$source_normals))
  expect_equal(nrow(fm$source_normals), fm$n_sources)
  # average-referenced leadfield has full sensor rank (n_electrodes - 1)
  expect_equal(qr(sweep(fm$leadfield, 2, colMeans(fm$leadfield)))$rank, 18L)
})

test_that("BEM forward runs and yields a finite leadfield of the right shape", {
  skip_if_not_installed("PhysioHeadModels")
  pe <- make_eeg(n_time = 60, n_channels = 19, sr = 250)
  fm <- eegForwardModel(pe, method = "bem", n_sources = 80)
  expect_equal(dim(fm$leadfield), c(19L, 240L))
  expect_true(all(is.finite(fm$leadfield)))
  expect_false(is.null(fm$source_normals))
})

test_that("realistic forward localizes a seeded cortical source at the correct depth", {
  skip_if_not_installed("PhysioHeadModels")
  pe0 <- make_eeg(n_time = 120, n_channels = 19, sr = 250)
  fm_real <- eegForwardModel(pe0, method = "sphere_analytic", n_sources = 200)
  fm_sph <- eegForwardModel(pe0, method = "spherical", n_sources = 200)
  loc <- function(fm, topo, seed) {
    set.seed(seed)
    pe <- pe0
    SummarizedExperiment::assay(pe, "raw") <-
      matrix(topo, 120, 19, byrow = TRUE) +
      matrix(stats::rnorm(120 * 19, sd = 0.02 * stats::sd(topo)), 120, 19)
    S <- S4Vectors::metadata(eegSourceEstimate(pe, fm, method = "sloreta"))$source
    as.numeric(fm$source_positions[which.max(.srcpow(S, fm$n_sources)), ])
  }
  radius <- function(p) sqrt(sum(p^2))
  set.seed(100); ks <- sample(seq_len(200), 12)
  E <- t(vapply(ks, function(k) {
    pt <- as.numeric(fm_real$source_positions[k, ])
    ori <- fm_real$source_normals[k, ]
    topo <- as.numeric(fm_real$leadfield[, ((k - 1) * 3 + 1):(k * 3)] %*% ori)
    er <- loc(fm_real, topo, k); es <- loc(fm_sph, topo, k)
    c(tot = sqrt(sum((pt - er)^2)),
      rad_real = abs(radius(er) - radius(pt)),
      rad_sph = abs(radius(es) - radius(pt)))
  }, numeric(3)))
  # localizes within a few cm (unit sphere ~ head radius)
  expect_lt(mean(E[, "tot"]), 0.5)
  # the anatomically-constrained cortical source space removes the random
  # cloud's depth ambiguity -> smaller radial (depth) error
  expect_lt(mean(E[, "rad_real"]), mean(E[, "rad_sph"]))
})

test_that("nyhead / fsaverage methods are gated on the downloaded data", {
  skip_if_not_installed("PhysioHeadModels")
  pe <- make_eeg(n_time = 60, n_channels = 19, sr = 250)
  expect_error(eegForwardModel(pe, method = "nyhead"),
               "not found|hdf5r|New York")
  expect_error(eegForwardModel(pe, method = "fsaverage"), "fsaverage")
})

test_that("legacy spherical / bem_simplified forward still work unchanged", {
  pe <- make_eeg(n_time = 60, n_channels = 19, sr = 250)
  fm <- eegForwardModel(pe, method = "spherical", n_sources = 50)
  expect_equal(dim(fm$leadfield), c(19L, 150L))
  expect_null(fm$source_normals)
  fm2 <- eegForwardModel(pe, method = "bem_simplified", n_sources = 50)
  expect_equal(dim(fm2$leadfield), c(19L, 150L))
})

library(testthat)
library(PhysioEEG)

test_that("eegDipoleFit recovers a planted dipole (position + orientation + GOF)", {
  pe <- make_eeg(n_time = 1, n_channels = 19, sr = 100)
  elec <- PhysioEEG:::.eeg_electrode_xyz(pe)
  elec <- elec / sqrt(rowSums(elec^2))
  q <- c(0.3, 0.2, 0.45); mo <- c(1, -0.5, 0.8)                 # planted dipole
  V <- as.vector(PhysioEEG:::.eeg_leadfield_at(elec, q, "spherical") %*% mo)
  SummarizedExperiment::assay(pe, defaultAssay(pe), withDimnames = FALSE) <-
    matrix(V, nrow = 1)

  fit <- eegDipoleFit(pe, n_grid = 8L)
  expect_s3_class(fit, "eeg_dipole_fit")
  expect_lt(sqrt(sum((fit$position - q)^2)), 0.1)               # position recovered
  expect_gt(fit$gof, 0.98)                                      # near-perfect fit
  # orientation recovered up to sign
  cosang <- abs(sum(fit$orientation * (mo / sqrt(sum(mo^2)))))
  expect_gt(cosang, 0.95)
  expect_output(print(fit), "eeg_dipole_fit")
})

test_that("eegDipoleFit picks the peak-GFP sample by default and runs bem_simplified", {
  pe <- make_eeg(n_time = 80, n_channels = 19, sr = 100)
  fit <- eegDipoleFit(pe, method = "bem_simplified")
  expect_s3_class(fit, "eeg_dipole_fit")
  expect_length(fit$time, 1L)
  expect_true(is.finite(fit$gof) && sqrt(sum(fit$position^2)) <= 1)
})

test_that("eegSourceEstimate dspm runs, is finite, and differs from MNE", {
  pe <- make_eeg(n_time = 40, n_channels = 19, sr = 100)
  fm <- eegForwardModel(pe, method = "spherical", n_sources = 50)
  sd <- eegSourceEstimate(pe, fm, method = "dspm")
  Sd <- S4Vectors::metadata(sd)$source
  Sm <- S4Vectors::metadata(eegSourceEstimate(pe, fm, method = "mne"))$source
  expect_true(is.matrix(Sd) && all(is.finite(Sd)))
  expect_equal(dim(Sd), c(40L, 50L * 3L))
  expect_false(isTRUE(all.equal(Sd, Sm)))
})

test_that("dSPM is MNE scaled per-source by a time-invariant factor", {
  pe <- make_eeg(n_time = 30, n_channels = 19, sr = 100)
  fm <- eegForwardModel(pe, method = "spherical", n_sources = 40)
  Sd <- t(S4Vectors::metadata(eegSourceEstimate(pe, fm, method = "dspm"))$source)
  Sm <- t(S4Vectors::metadata(eegSourceEstimate(pe, fm, method = "mne"))$source)
  rows <- which(rowSums(abs(Sm)) > 1e-6)[1:6]
  for (i in rows) {
    ratio <- Sd[i, ] / Sm[i, ]
    expect_lt(stats::sd(ratio) / (abs(mean(ratio)) + 1e-12), 1e-6)   # constant across time
  }
})

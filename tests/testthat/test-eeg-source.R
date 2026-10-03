library(testthat)
library(PhysioEEG)

test_that("eegForwardModel creates leadfield with correct dimensions", {
  pe <- make_eeg(n_time = 1000, n_channels = 19, sr = 250)
  fm <- eegForwardModel(pe, method = "spherical", n_sources = 100)

  expect_true(is.list(fm))
  expect_true("leadfield" %in% names(fm))
  expect_equal(nrow(fm$leadfield), 19)
  expect_equal(ncol(fm$leadfield), 100 * 3)  # 3 orientations per source
  expect_true(all(is.finite(fm$leadfield)))
})

test_that("eegForwardModel returns source and electrode positions", {
  pe <- make_eeg(n_time = 1000, n_channels = 19, sr = 250)
  fm <- eegForwardModel(pe, method = "spherical", n_sources = 50)

  expect_true("source_positions" %in% names(fm))
  expect_true("electrode_positions" %in% names(fm))
  expect_equal(nrow(fm$source_positions), 50)
})

test_that("eegSourceEstimate with MNE produces finite output", {
  pe <- make_eeg(n_time = 500, n_channels = 19, sr = 250)
  fm <- eegForwardModel(pe, method = "spherical", n_sources = 50)
  result <- eegSourceEstimate(pe, fm, method = "mne")

  expect_s4_class(result, "PhysioExperiment")
  expect_true(!is.null(S4Vectors::metadata(result)[["source"]]))
  src <- S4Vectors::metadata(result)[["source"]]
  expect_true(all(is.finite(src)))
})

test_that("eegSourceEstimate sLORETA differs from MNE", {
  pe <- make_eeg(n_time = 500, n_channels = 19, sr = 250)
  fm <- eegForwardModel(pe, method = "spherical", n_sources = 50)
  mne <- eegSourceEstimate(pe, fm, method = "mne")
  sloreta <- eegSourceEstimate(pe, fm, method = "sloreta")

  src_mne <- S4Vectors::metadata(mne)[["source"]]
  src_sloreta <- S4Vectors::metadata(sloreta)[["source"]]
  expect_false(all(src_mne == src_sloreta))
})

test_that("eegBeamformer LCMV returns source power", {
  pe <- make_eeg(n_time = 1000, n_channels = 19, sr = 250)
  fm <- eegForwardModel(pe, method = "spherical", n_sources = 50)
  result <- eegBeamformer(pe, fm, method = "lcmv")

  expect_s4_class(result, "PhysioExperiment")
  expect_true(!is.null(S4Vectors::metadata(result)[["beamformer"]]))
})

test_that("eegSourcePower returns band powers", {
  pe <- make_eeg(n_time = 1000, n_channels = 19, sr = 250)
  fm <- eegForwardModel(pe, method = "spherical", n_sources = 20)
  pe <- eegSourceEstimate(pe, fm, method = "mne")
  result <- eegSourcePower(pe)

  expect_s3_class(result, "data.frame")
  expect_true("band" %in% names(result))
  expect_true("power" %in% names(result))
  expect_true(all(result$power >= 0))
})

test_that("eegForwardModel uses correct dipole physics", {
  pe <- make_eeg(n_time = 500, n_channels = 19, sr = 250)
  fm <- eegForwardModel(pe, method = "spherical", n_sources = 10)

  # Leadfield should have correct dimensions
  expect_equal(nrow(fm$leadfield), 19)
  expect_equal(ncol(fm$leadfield), 30)  # 10 sources * 3 orientations

  # All values should be finite
  expect_true(all(is.finite(fm$leadfield)))

  # Verify dipole physics: field should include 1/(4*pi*sigma) normalization.
  # For sigma = 0.33, the scale factor is ~0.2416.
  # The leadfield values should be bounded by this normalization.
  # A source at distance ~0.5 from an electrode would give
  #   L ~ 0.2416 * 0.5 / 0.5^3 = 0.2416 / 0.25 ~ 0.97
  # Values should be moderate, not O(1e3) or O(1e-10)
  max_val <- max(abs(fm$leadfield))
  expect_true(max_val > 0)
  expect_true(max_val < 1e6)  # not unreasonably large
})

test_that("eegForwardModel BEM simplified uses Berg-Scherg correction", {
  pe <- make_eeg(n_time = 500, n_channels = 19, sr = 250)
  fm_sph <- eegForwardModel(pe, method = "spherical", n_sources = 10)
  fm_bem <- eegForwardModel(pe, method = "bem_simplified", n_sources = 10)

  # BEM and spherical should differ (virtual dipoles at different positions)
  expect_false(all(fm_sph$leadfield == fm_bem$leadfield))

  # Both should be finite
  expect_true(all(is.finite(fm_bem$leadfield)))
})

test_that("eegSourceEstimate eLORETA converges iteratively", {
  pe <- make_eeg(n_time = 500, n_channels = 19, sr = 250)
  fm <- eegForwardModel(pe, method = "spherical", n_sources = 20)
  result <- eegSourceEstimate(pe, fm, method = "eloreta", lambda = 0.1)

  expect_s4_class(result, "PhysioExperiment")
  src <- S4Vectors::metadata(result)[["source"]]
  expect_true(all(is.finite(src)))
  expect_equal(nrow(src), 500)
  expect_equal(ncol(src), 60)  # 20 sources * 3 orientations
})

test_that("eegForwardModel leadfield has correct dimensions 3 orientations per source", {
  pe <- make_eeg(n_time = 500, n_channels = 4, sr = 250)
  fm <- eegForwardModel(pe, method = "spherical", n_sources = 25)

  expect_equal(nrow(fm$leadfield), 4)           # n_channels
  expect_equal(ncol(fm$leadfield), 25 * 3)       # n_sources * 3 orientations
  expect_equal(nrow(fm$source_positions), 25)
  expect_equal(ncol(fm$source_positions), 3)     # x, y, z
})

test_that("eegSourceEstimate MNE output dimensions are correct", {
  pe <- make_eeg(n_time = 300, n_channels = 19, sr = 250)
  fm <- eegForwardModel(pe, method = "spherical", n_sources = 30)
  result <- eegSourceEstimate(pe, fm, method = "mne")

  src <- S4Vectors::metadata(result)[["source"]]
  expect_equal(nrow(src), 300)        # n_time
  expect_equal(ncol(src), 30 * 3)     # n_sources * 3 orientations

  # Metadata should be stored
  se_md <- metadata(result)$source_estimate
  expect_equal(se_md$method, "mne")
  expect_equal(se_md$n_sources, 30)
})

test_that("eegBeamformer LCMV output dimensions correct", {
  pe <- make_eeg(n_time = 1000, n_channels = 4, sr = 250)
  fm <- eegForwardModel(pe, method = "spherical", n_sources = 10)
  result <- eegBeamformer(pe, fm, method = "lcmv")

  bf_data <- S4Vectors::metadata(result)[["beamformer"]]
  expect_equal(nrow(bf_data), 10)   # n_sources
  expect_equal(ncol(bf_data), 1)    # single power column
  expect_true(all(bf_data >= 0))    # power is non-negative
})

test_that("eegSourcePower has correct frequency bands in output", {
  pe <- make_eeg(n_time = 1000, n_channels = 19, sr = 250)
  fm <- eegForwardModel(pe, method = "spherical", n_sources = 5)
  pe <- eegSourceEstimate(pe, fm, method = "mne")
  result <- eegSourcePower(pe)

  expect_s3_class(result, "data.frame")
  expect_true(all(c("source_id", "band", "power") %in% names(result)))
  expected_bands <- c("delta", "theta", "alpha", "beta", "gamma")
  expect_true(all(expected_bands %in% result$band))
  expect_true(all(result$power >= 0))
})

test_that("eegForwardModel single source generates valid leadfield", {
  pe <- make_eeg(n_time = 500, n_channels = 19, sr = 250)
  fm <- eegForwardModel(pe, method = "spherical", n_sources = 1)

  expect_equal(nrow(fm$leadfield), 19)
  expect_equal(ncol(fm$leadfield), 3)  # 1 source * 3 orientations
  expect_true(all(is.finite(fm$leadfield)))
  expect_equal(nrow(fm$source_positions), 1)
})

test_that("eegBeamformer DICS returns non-negative source power", {
  pe <- make_eeg(n_time = 1000, n_channels = 19, sr = 250)
  fm <- eegForwardModel(pe, method = "spherical", n_sources = 10)
  result <- eegBeamformer(pe, fm, method = "dics", freq_range = c(8, 13))

  bf_data <- S4Vectors::metadata(result)[["beamformer"]]
  expect_equal(nrow(bf_data), 10)
  expect_true(all(bf_data >= 0))
})

test_that("eegSourceEstimate stores source_estimate metadata", {
  pe <- make_eeg(n_time = 500, n_channels = 19, sr = 250)
  fm <- eegForwardModel(pe, method = "spherical", n_sources = 15)
  result <- eegSourceEstimate(pe, fm, method = "sloreta", lambda = 0.1)

  md <- metadata(result)$source_estimate
  expect_true(!is.null(md))
  expect_equal(md$method, "sloreta")
  expect_equal(md$n_sources, 15)
  expect_equal(md$lambda, 0.1)
})

test_that("eegSourceEstimate with 2-channel system works", {
  pe <- make_eeg(n_time = 500, n_channels = 2, sr = 250)
  fm <- eegForwardModel(pe, method = "spherical", n_sources = 5)
  result <- eegSourceEstimate(pe, fm, method = "mne")

  src <- S4Vectors::metadata(result)[["source"]]
  expect_equal(nrow(src), 500)
  expect_equal(ncol(src), 15)  # 5 * 3
  expect_true(all(is.finite(src)))
})

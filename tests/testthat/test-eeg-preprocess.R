library(testthat)
library(PhysioEEG)

# =============================================================================
# eegFilter tests
# =============================================================================

test_that("eegFilter applies bandpass filter correctly", {
  pe <- make_eeg(n_time = 5000, n_channels = 4, sr = 500)
  pe_filt <- eegFilter(pe, lowcut = 1, highcut = 40)

  expect_s4_class(pe_filt, "PhysioExperiment")
  expect_true("filtered" %in% SummarizedExperiment::assayNames(pe_filt))
  expect_equal(dim(SummarizedExperiment::assay(pe_filt, "filtered")),
               dim(SummarizedExperiment::assay(pe, "raw")))
})

test_that("eegFilter applies highpass filter correctly", {
  pe <- make_eeg(n_time = 5000, n_channels = 4, sr = 500)
  pe_hp <- eegFilter(pe, lowcut = 1)

  expect_s4_class(pe_hp, "PhysioExperiment")
  expect_true("filtered" %in% SummarizedExperiment::assayNames(pe_hp))

  raw <- SummarizedExperiment::assay(pe, "raw")
  filt <- SummarizedExperiment::assay(pe_hp, "filtered")
  # Filtered should be different from raw

  expect_false(identical(raw, filt))
})

test_that("eegFilter applies lowpass filter correctly", {
  pe <- make_eeg(n_time = 5000, n_channels = 4, sr = 500)
  pe_lp <- eegFilter(pe, highcut = 40)

  expect_s4_class(pe_lp, "PhysioExperiment")
  expect_true("filtered" %in% SummarizedExperiment::assayNames(pe_lp))
  expect_equal(dim(SummarizedExperiment::assay(pe_lp, "filtered")),
               dim(SummarizedExperiment::assay(pe, "raw")))
})

test_that("eegFilter applies notch filter correctly", {
  pe <- make_eeg(n_time = 5000, n_channels = 4, sr = 500)
  pe_notch <- eegFilter(pe, notch = 50)

  expect_s4_class(pe_notch, "PhysioExperiment")
  expect_true("filtered" %in% SummarizedExperiment::assayNames(pe_notch))
})

test_that("eegFilter supports combined bandpass + notch", {
  pe <- make_eeg(n_time = 5000, n_channels = 4, sr = 500)
  pe_combo <- eegFilter(pe, lowcut = 1, highcut = 40, notch = 50)

  expect_s4_class(pe_combo, "PhysioExperiment")
  filt <- SummarizedExperiment::assay(pe_combo, "filtered")
  expect_equal(nrow(filt), 5000)
  expect_equal(ncol(filt), 4)
})

test_that("eegFilter IIR mode works", {
  pe <- make_eeg(n_time = 5000, n_channels = 4, sr = 500)
  pe_iir <- eegFilter(pe, lowcut = 1, highcut = 40, method = "iir", order = 4)

  expect_s4_class(pe_iir, "PhysioExperiment")
  expect_true("filtered" %in% SummarizedExperiment::assayNames(pe_iir))
})

test_that("eegFilter errors on invalid parameters", {
  pe <- make_eeg(n_time = 5000, n_channels = 4, sr = 500)

  # No filter specified

  expect_error(eegFilter(pe), "At least one of")

  # lowcut >= highcut
  expect_error(eegFilter(pe, lowcut = 40, highcut = 1), "lowcut must be less than highcut")

  # lowcut >= Nyquist
  expect_error(eegFilter(pe, lowcut = 260))
})

test_that("eegFilter writes to custom output assay", {
  pe <- make_eeg(n_time = 5000, n_channels = 4, sr = 500)
  pe_filt <- eegFilter(pe, lowcut = 1, highcut = 40, output_assay = "my_filtered")

  expect_true("my_filtered" %in% SummarizedExperiment::assayNames(pe_filt))
})


# =============================================================================
# eegRereference tests
# =============================================================================

test_that("eegRereference average reference works", {
  pe <- make_eeg(n_time = 2000, n_channels = 19, sr = 500)
  pe_ref <- eegRereference(pe, ref_type = "average")

  expect_s4_class(pe_ref, "PhysioExperiment")
  expect_true("rereferenced" %in% SummarizedExperiment::assayNames(pe_ref))

  reref <- SummarizedExperiment::assay(pe_ref, "rereferenced")
  # After average reference, row means should be approximately zero
  row_means <- rowMeans(reref)
  expect_true(all(abs(row_means) < 1e-10))
})

test_that("eegRereference Cz reference works", {
  pe <- make_eeg(n_time = 2000, n_channels = 19, sr = 500)
  pe_cz <- eegRereference(pe, ref_type = "cz")

  expect_s4_class(pe_cz, "PhysioExperiment")
  reref <- SummarizedExperiment::assay(pe_cz, "rereferenced")

  # Cz channel should be all zeros after Cz reference
  cd <- SummarizedExperiment::colData(pe_cz)
  cz_idx <- which(as.character(cd$label) == "Cz")
  expect_true(all(abs(reref[, cz_idx]) < 1e-10))
})

test_that("eegRereference channel reference works", {
  pe <- make_eeg(n_time = 2000, n_channels = 19, sr = 500)
  pe_ch <- eegRereference(pe, ref_type = "channel", ref_channels = "Fz")

  expect_s4_class(pe_ch, "PhysioExperiment")
  reref <- SummarizedExperiment::assay(pe_ch, "rereferenced")

  # Fz should be zeros
  cd <- SummarizedExperiment::colData(pe_ch)
  fz_idx <- which(as.character(cd$label) == "Fz")
  expect_true(all(abs(reref[, fz_idx]) < 1e-10))
})

test_that("eegRereference stores reference info in metadata", {
  pe <- make_eeg(n_time = 2000, n_channels = 19, sr = 500)
  pe_ref <- eegRereference(pe, ref_type = "average")

  md <- S4Vectors::metadata(pe_ref)
  expect_true(!is.null(md$reference))
  expect_equal(md$reference$type, "average")
})

test_that("eegRereference errors on missing channel", {
  pe <- make_eeg(n_time = 2000, n_channels = 19, sr = 500)

  expect_error(eegRereference(pe, ref_type = "channel"),
               "ref_channels must be specified")
  expect_error(eegRereference(pe, ref_type = "channel",
                               ref_channels = "NONEXISTENT"),
               "No matching channels found")
})


# =============================================================================
# eegBadChannels tests
# =============================================================================

test_that("eegBadChannels detects flat channels", {
  pe <- make_eeg(n_time = 2000, n_channels = 4, sr = 500)
  # Make one channel flat
  data <- SummarizedExperiment::assay(pe, "raw")
  data[, 2] <- 0
  SummarizedExperiment::assay(pe, "raw") <- data

  bad_df <- eegBadChannels(pe, method = "flat")

  expect_s3_class(bad_df, "data.frame")
  expect_true("channel" %in% colnames(bad_df))
  expect_true("is_bad" %in% colnames(bad_df))
  expect_true("reason" %in% colnames(bad_df))
  expect_true("score" %in% colnames(bad_df))

  expect_true(bad_df$is_bad[2])
  expect_true(grepl("flat", bad_df$reason[2]))
})

test_that("eegBadChannels detects noisy channels", {
  pe <- make_eeg(n_time = 2000, n_channels = 4, sr = 500)
  # Make one channel very noisy
  data <- SummarizedExperiment::assay(pe, "raw")
  data[, 3] <- data[, 3] * 100
  SummarizedExperiment::assay(pe, "raw") <- data

  bad_df <- eegBadChannels(pe, method = "noise", noise_threshold = 3)

  expect_true(bad_df$is_bad[3])
  expect_true(grepl("noise", bad_df$reason[3]))
})

test_that("eegBadChannels detects poorly correlated channels", {
  pe <- make_eeg(n_time = 2000, n_channels = 4, sr = 500)
  # Replace one channel with independent noise
  data <- SummarizedExperiment::assay(pe, "raw")
  set.seed(42)
  data[, 1] <- rnorm(2000, sd = 100)
  SummarizedExperiment::assay(pe, "raw") <- data

  bad_df <- eegBadChannels(pe, method = "correlation", corr_threshold = 0.3)

  # Channel 1 should have low correlation with others
  expect_true(bad_df$is_bad[1])
  expect_true(grepl("correlation", bad_df$reason[1]))
})

test_that("eegBadChannels 'all' method runs all checks", {
  pe <- make_eeg(n_time = 2000, n_channels = 4, sr = 500)
  data <- SummarizedExperiment::assay(pe, "raw")
  data[, 2] <- 0  # flat channel
  SummarizedExperiment::assay(pe, "raw") <- data

  bad_df <- eegBadChannels(pe, method = "all")

  expect_equal(nrow(bad_df), 4)
  expect_true(bad_df$is_bad[2])
})


# =============================================================================
# eegInterpolate tests
# =============================================================================

test_that("eegInterpolate spline method works", {
  pe <- make_eeg(n_time = 2000, n_channels = 19, sr = 500)
  pe <- eegMontage(pe, system = "10-20")

  # Corrupt channel and interpolate
  data <- SummarizedExperiment::assay(pe, "raw")
  data[, 3] <- 0  # zero out F7
  SummarizedExperiment::assay(pe, "raw") <- data

  pe_interp <- eegInterpolate(pe, bad_channels = "F7", method = "spline")

  expect_s4_class(pe_interp, "PhysioExperiment")
  expect_true("interpolated" %in% SummarizedExperiment::assayNames(pe_interp))

  interp_data <- SummarizedExperiment::assay(pe_interp, "interpolated")
  # Interpolated channel should no longer be flat
  expect_gt(var(interp_data[, 3]), 1e-6)
})

test_that("eegInterpolate nearest method works", {
  pe <- make_eeg(n_time = 2000, n_channels = 19, sr = 500)
  pe <- eegMontage(pe, system = "10-20")

  data <- SummarizedExperiment::assay(pe, "raw")
  data[, 3] <- 0
  SummarizedExperiment::assay(pe, "raw") <- data

  pe_interp <- eegInterpolate(pe, bad_channels = "F7", method = "nearest")

  interp_data <- SummarizedExperiment::assay(pe_interp, "interpolated")
  expect_gt(var(interp_data[, 3]), 1e-6)
})

test_that("eegInterpolate errors without positions", {
  pe <- make_eeg(n_time = 2000, n_channels = 4, sr = 500)

  expect_error(eegInterpolate(pe, bad_channels = "Fp1"),
               "Electrode positions.*required")
})


# =============================================================================
# eegEpoch tests
# =============================================================================

test_that("eegEpoch creates correct 3D output", {
  pe <- make_eeg(n_time = 10000, n_channels = 4, sr = 500)
  events <- data.frame(onset_sec = c(2, 4, 6, 8))
  pe_ep <- eegEpoch(pe, events, limits = c(-0.2, 0.8))

  expect_s4_class(pe_ep, "PhysioExperiment")
  expect_true(!is.null(S4Vectors::metadata(pe_ep)[["epoched"]]))

  ep_data <- S4Vectors::metadata(pe_ep)[["epoched"]]
  d <- dim(ep_data)
  expect_equal(length(d), 3)  # 3D array

  # Expected epoch length: 0.2 + 0.8 = 1.0 sec at 500 Hz = 501 samples
  expected_len <- as.integer(round(0.8 * 500)) - as.integer(round(-0.2 * 500)) + 1
  expect_equal(d[1], expected_len)
  expect_equal(d[2], 4)  # channels
  expect_equal(d[3], 4)  # epochs
})

test_that("eegEpoch applies baseline correction", {
  pe <- make_eeg(n_time = 10000, n_channels = 4, sr = 500)
  events <- data.frame(onset_sec = c(2, 4, 6, 8))

  pe_bl <- eegEpoch(pe, events, limits = c(-0.2, 0.8), baseline = c(-0.2, 0))
  pe_no_bl <- eegEpoch(pe, events, limits = c(-0.2, 0.8), baseline = NULL,
                         output_assay = "no_bl")

  bl_data <- S4Vectors::metadata(pe_bl)[["epoched"]]
  # Baseline corrected data should differ from non-baseline corrected
  # (unless baseline mean happens to be exactly zero)
  no_bl_data <- S4Vectors::metadata(pe_no_bl)[["no_bl"]]
  expect_false(identical(bl_data, no_bl_data))
})

test_that("eegEpoch handles edge events", {
  pe <- make_eeg(n_time = 5000, n_channels = 4, sr = 500)
  # Events at very beginning and end should be excluded
  events <- data.frame(onset_sec = c(0.01, 2, 4, 9.99))

  expect_message(
    pe_ep <- eegEpoch(pe, events, limits = c(-0.2, 0.8)),
    "Excluding"
  )

  ep_data <- S4Vectors::metadata(pe_ep)[["epoched"]]
  # Only middle events should remain
  expect_lt(dim(ep_data)[3], 4)
})

test_that("eegEpoch accepts integer sample indices", {
  pe <- make_eeg(n_time = 5000, n_channels = 4, sr = 500)
  events <- c(500L, 1500L, 2500L)

  pe_ep <- eegEpoch(pe, events, limits = c(-0.2, 0.8))

  expect_s4_class(pe_ep, "PhysioExperiment")
  ep_data <- S4Vectors::metadata(pe_ep)[["epoched"]]
  expect_equal(dim(ep_data)[3], 3)
})

test_that("eegEpoch stores epoch info in metadata", {
  pe <- make_eeg(n_time = 10000, n_channels = 4, sr = 500)
  events <- data.frame(onset_sec = c(2, 4, 6))
  pe_ep <- eegEpoch(pe, events, limits = c(-0.2, 0.8))

  md <- S4Vectors::metadata(pe_ep)
  expect_true(!is.null(md$epoch_events))
  expect_equal(md$epoch_events$n_epochs, 3)
})


# =============================================================================
# eegArtifactReject tests
# =============================================================================

test_that("eegArtifactReject threshold method works", {
  pe <- make_eeg(n_time = 10000, n_channels = 4, sr = 500)
  events <- data.frame(onset_sec = c(2, 4, 6, 8))
  pe_ep <- eegEpoch(pe, events, limits = c(-0.2, 0.8), baseline = NULL)

  # Inject a large artifact into one epoch
  ep_data <- S4Vectors::metadata(pe_ep)[["epoched"]]
  ep_data[50, 1, 2] <- 500  # epoch 2, huge value
  md <- S4Vectors::metadata(pe_ep)
  md[["epoched"]] <- ep_data
  S4Vectors::metadata(pe_ep) <- md

  pe_clean <- eegArtifactReject(pe_ep, method = "threshold",
                                 threshold_uv = 100,
                                 assay_name = "epoched")

  clean_data <- S4Vectors::metadata(pe_clean)[["clean"]]
  # Should have fewer epochs
  expect_lt(dim(clean_data)[3], dim(ep_data)[3])
})

test_that("eegArtifactReject gradient method works", {
  pe <- make_eeg(n_time = 10000, n_channels = 4, sr = 500)
  events <- data.frame(onset_sec = c(2, 4, 6, 8))
  pe_ep <- eegEpoch(pe, events, limits = c(-0.2, 0.8), baseline = NULL)

  # Inject a sharp transient into one epoch
  ep_data <- S4Vectors::metadata(pe_ep)[["epoched"]]
  ep_data[50, 1, 3] <- 0
  ep_data[51, 1, 3] <- 10000  # huge gradient
  md <- S4Vectors::metadata(pe_ep)
  md[["epoched"]] <- ep_data
  S4Vectors::metadata(pe_ep) <- md

  pe_clean <- eegArtifactReject(pe_ep, method = "gradient",
                                 gradient_uv_ms = 50,
                                 assay_name = "epoched")

  clean_data <- S4Vectors::metadata(pe_clean)[["clean"]]
  expect_lt(dim(clean_data)[3], dim(ep_data)[3])
})

test_that("eegArtifactReject joint probability method works", {
  pe <- make_eeg(n_time = 10000, n_channels = 4, sr = 500)
  events <- data.frame(onset_sec = c(2, 4, 6, 8))
  pe_ep <- eegEpoch(pe, events, limits = c(-0.2, 0.8), baseline = NULL)

  # Make one epoch have very high power
  ep_data <- S4Vectors::metadata(pe_ep)[["epoched"]]
  ep_data[, , 1] <- ep_data[, , 1] * 100
  md <- S4Vectors::metadata(pe_ep)
  md[["epoched"]] <- ep_data
  S4Vectors::metadata(pe_ep) <- md

  pe_clean <- eegArtifactReject(pe_ep, method = "joint_probability",
                                 jp_threshold = 1,
                                 assay_name = "epoched")

  clean_data <- S4Vectors::metadata(pe_clean)[["clean"]]
  expect_lt(dim(clean_data)[3], dim(ep_data)[3])
})

test_that("eegArtifactReject stores artifact log", {
  pe <- make_eeg(n_time = 10000, n_channels = 4, sr = 500)
  events <- data.frame(onset_sec = c(2, 4, 6, 8))
  pe_ep <- eegEpoch(pe, events, limits = c(-0.2, 0.8), baseline = NULL)

  ep_data <- S4Vectors::metadata(pe_ep)[["epoched"]]
  ep_data[50, 1, 2] <- 500
  md <- S4Vectors::metadata(pe_ep)
  md[["epoched"]] <- ep_data
  S4Vectors::metadata(pe_ep) <- md

  pe_clean <- eegArtifactReject(pe_ep, method = "threshold",
                                 threshold_uv = 100,
                                 assay_name = "epoched")

  md <- S4Vectors::metadata(pe_clean)
  expect_true(!is.null(md$artifact_log))
  expect_true(!is.null(md$artifact_summary))
  expect_gt(md$artifact_summary$n_rejected, 0)
})

test_that("eegArtifactReject errors on non-3D data", {
  pe <- make_eeg(n_time = 5000, n_channels = 4, sr = 500)

  expect_error(eegArtifactReject(pe, method = "threshold"),
               "3D epoched data")
})


# =============================================================================
# eegMontage tests
# =============================================================================

test_that("eegMontage assigns 10-20 positions correctly", {
  pe <- make_eeg(n_time = 2000, n_channels = 19, sr = 500)
  pe_m <- eegMontage(pe, system = "10-20")

  cd <- SummarizedExperiment::colData(pe_m)
  expect_true("pos_x" %in% colnames(cd))
  expect_true("pos_y" %in% colnames(cd))
  expect_true("pos_z" %in% colnames(cd))

  # All 19 channels should have positions
  expect_equal(sum(!is.na(cd$pos_x)), 19)
  expect_equal(sum(!is.na(cd$pos_y)), 19)
  expect_equal(sum(!is.na(cd$pos_z)), 19)
})

test_that("eegMontage assigns 10-10 positions", {
  pe <- make_eeg(n_time = 2000, n_channels = 19, sr = 500)
  pe_m <- eegMontage(pe, system = "10-10")

  cd <- SummarizedExperiment::colData(pe_m)
  # All 10-20 channels are in 10-10 as well
  expect_equal(sum(!is.na(cd$pos_x)), 19)
})

test_that("eegMontage stores montage info in metadata", {
  pe <- make_eeg(n_time = 2000, n_channels = 19, sr = 500)
  pe_m <- eegMontage(pe, system = "10-20")

  md <- S4Vectors::metadata(pe_m)
  expect_true(!is.null(md$montage))
  expect_equal(md$montage$system, "10-20")
  expect_equal(md$montage$n_matched, 19)
})

test_that("eegMontage custom positions work", {
  pe <- make_eeg(n_time = 2000, n_channels = 4, sr = 500)
  cd <- SummarizedExperiment::colData(pe)
  labels <- as.character(cd$label)

  custom_pos <- data.frame(
    label = labels,
    pos_x = c(0.1, 0.2, 0.3, 0.4),
    pos_y = c(0.5, 0.6, 0.7, 0.8),
    pos_z = c(0.0, 0.1, 0.2, 0.3),
    stringsAsFactors = FALSE
  )

  pe_m <- eegMontage(pe, system = "custom", positions = custom_pos)
  cd_out <- SummarizedExperiment::colData(pe_m)
  expect_equal(as.numeric(cd_out$pos_x), c(0.1, 0.2, 0.3, 0.4))
})

test_that("eegMontage errors on no matching channels", {
  pe <- make_eeg(n_time = 2000, n_channels = 4, sr = 500)
  # Rename labels to something unrecognizable
  cd <- SummarizedExperiment::colData(pe)
  cd$label <- paste0("XX", seq_len(4))
  SummarizedExperiment::colData(pe) <- cd

  expect_error(eegMontage(pe, system = "10-20"), "No channel labels matched")
})


# =============================================================================
# eegPreprocess (full pipeline) tests
# =============================================================================

test_that("eegPreprocess full pipeline runs without error", {
  pe <- make_eeg(n_time = 10000, n_channels = 19, sr = 500)
  pe <- eegMontage(pe, system = "10-20")

  events <- data.frame(onset_sec = c(2, 4, 6, 8))

  pe_proc <- eegPreprocess(
    pe,
    filter = TRUE, lowcut = 1, highcut = 40,
    rereference = TRUE, ref_type = "average",
    bad_channels = FALSE, interpolate = FALSE,
    epoch = TRUE, events = events,
    epoch_limits = c(-0.2, 0.8),
    baseline = c(-0.2, 0),
    artifact_reject = TRUE, threshold_uv = 200,
    verbose = FALSE
  )

  expect_s4_class(pe_proc, "PhysioExperiment")

  # Check that preprocessing log was created
  md <- S4Vectors::metadata(pe_proc)
  expect_true(!is.null(md$preprocess_log))
  expect_gt(length(md$preprocess_log), 0)
})

test_that("eegPreprocess filter-only pipeline works", {
  pe <- make_eeg(n_time = 5000, n_channels = 4, sr = 500)

  pe_proc <- eegPreprocess(
    pe,
    filter = TRUE, lowcut = 1, highcut = 40,
    rereference = FALSE, bad_channels = FALSE,
    interpolate = FALSE, epoch = FALSE,
    verbose = FALSE
  )

  expect_s4_class(pe_proc, "PhysioExperiment")
  expect_true("filtered" %in% SummarizedExperiment::assayNames(pe_proc))
})

test_that("eegPreprocess errors when epoch=TRUE without events", {
  pe <- make_eeg(n_time = 5000, n_channels = 4, sr = 500)

  expect_error(
    eegPreprocess(pe, epoch = TRUE, events = NULL, verbose = FALSE),
    "events must be provided"
  )
})

test_that("eegPreprocess records final assay in metadata", {
  pe <- make_eeg(n_time = 5000, n_channels = 4, sr = 500)

  pe_proc <- eegPreprocess(
    pe,
    filter = TRUE, lowcut = 1, highcut = 40,
    rereference = TRUE, ref_type = "average",
    bad_channels = FALSE, interpolate = FALSE,
    epoch = FALSE, verbose = FALSE
  )

  md <- S4Vectors::metadata(pe_proc)
  expect_true(!is.null(md$preprocess_final_assay))
  expect_equal(md$preprocess_final_assay, "rereferenced")
})

library(testthat)
library(PhysioEEG)

# EEG matrix with a shared brain rhythm plus one injected high-amplitude channel.
make_bad_channel_pe <- function(seed = 1, nch = 19, nt = 3000, bad = 7,
                                artifact_sd = 200) {
  set.seed(seed)
  base <- matrix(rnorm(nt * nch, sd = 10), nt, nch) +
    20 * sin(2 * pi * 10 * seq_len(nt) / 500)
  base[, bad] <- base[, bad] + rnorm(nt, sd = artifact_sd)
  pe <- PhysioExperiment(
    assays = list(raw = base),
    colData = S4Vectors::DataFrame(label = paste0("E", seq_len(nch)),
                                   type = rep("EEG", nch)),
    samplingRate = 500)
  list(pe = pe, bad = bad, base = base)
}

rms <- function(m) sqrt(mean(m^2))

test_that("robust reference matches the clean average and beats plain average", {
  fx <- make_bad_channel_pe()
  good <- setdiff(seq_len(ncol(fx$base)), fx$bad)
  clean_avg <- rowMeans(fx$base[, good])            # ground-truth reference
  clean_ref <- fx$base[, good] - clean_avg

  rob <- SummarizedExperiment::assay(
    eegRereference(fx$pe, ref_type = "robust"), "rereferenced")
  # good-channel RMS within 1% of the ground-truth clean average reference
  expect_lt(rms(rob[, good] - clean_ref) / rms(clean_ref), 0.01)

  # the plain rowMeans reference is contaminated by the bad channel
  plain <- SummarizedExperiment::assay(
    eegRereference(fx$pe, ref_type = "average"), "rereferenced")
  expect_gt(rms(plain[, good] - clean_ref), rms(rob[, good] - clean_ref))
  # and its deviation exceeds the per-channel spread of the injected artifact
  expect_gt(rms(plain[, good] - clean_ref), 200 / ncol(fx$base) * 0.5)
})

test_that("robust reference excludes the injected channel and records provenance", {
  fx <- make_bad_channel_pe()
  pr <- eegRereference(fx$pe, ref_type = "robust")
  ref <- S4Vectors::metadata(pr)$reference
  expect_equal(ref$type, "robust")
  expect_true(paste0("E", fx$bad) %in% ref$excluded_channels)
  prov <- PhysioCore::provenance(pr)
  expect_true(any(prov$step == "eegRereference"))
})

test_that("robust reference is idempotent (< 1e-6)", {
  fx <- make_bad_channel_pe()
  pr <- eegRereference(fx$pe, ref_type = "robust")
  rob <- SummarizedExperiment::assay(pr, "rereferenced")
  pe2 <- fx$pe
  SummarizedExperiment::assay(pe2, "raw") <- rob
  rob2 <- SummarizedExperiment::assay(
    eegRereference(pe2, ref_type = "robust"), "rereferenced")
  expect_lt(max(abs(rob2 - rob)), 1e-6)
})

test_that("median reference is robust and idempotent", {
  fx <- make_bad_channel_pe()
  med <- SummarizedExperiment::assay(
    eegRereference(fx$pe, ref_type = "median"), "rereferenced")
  # subtracting the row median again changes nothing
  pe2 <- fx$pe
  SummarizedExperiment::assay(pe2, "raw") <- med
  med2 <- SummarizedExperiment::assay(
    eegRereference(pe2, ref_type = "median"), "rereferenced")
  expect_lt(max(abs(med2 - med)), 1e-6)
  # median of a median-referenced signal is ~0 across channels
  expect_lt(max(abs(apply(med, 1, stats::median))), 1e-6)
})

test_that("existing reference types still work", {
  fx <- make_bad_channel_pe()
  av <- eegRereference(fx$pe, ref_type = "average")
  expect_true("rereferenced" %in% SummarizedExperiment::assayNames(av))
})

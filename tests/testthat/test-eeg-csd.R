library(testthat)
library(PhysioEEG)

mkV <- function(v, n_ch = 19, sr = 250) {
  m <- if (is.matrix(v)) v else matrix(v, ncol = n_ch)
  pe <- make_eeg(n_time = nrow(m), n_channels = n_ch, sr = sr)
  SummarizedExperiment::assay(pe, defaultAssay(pe), withDimnames = FALSE) <- m
  pe
}

test_that("eegSurfaceLaplacian is reference-free (invariant to a constant offset)", {
  set.seed(1); V <- matrix(rnorm(19 * 5), ncol = 19)
  c1 <- SummarizedExperiment::assay(eegSurfaceLaplacian(mkV(V)), "csd")
  c2 <- SummarizedExperiment::assay(eegSurfaceLaplacian(mkV(V + 7)), "csd")
  expect_lt(max(abs(c1 - c2)), 1e-8)          # CSD removes the reference constant
})

test_that("eegSurfaceLaplacian acts as a spatial high-pass (focal >> smooth)", {
  pe0 <- make_eeg(n_time = 1, n_channels = 19, sr = 250)
  pos <- PhysioEEG:::.eeg_electrode_xyz(pe0)
  pos <- pos / sqrt(rowSums(pos^2))
  smooth <- matrix(pos[, 3], nrow = 1)                       # linear in z: very smooth
  focal <- matrix(exp(-((seq_len(19) - 10)^2) / 2), nrow = 1) # a bump at one site
  es <- sum(SummarizedExperiment::assay(eegSurfaceLaplacian(mkV(smooth)), "csd")^2)
  ef <- sum(SummarizedExperiment::assay(eegSurfaceLaplacian(mkV(focal)), "csd")^2)
  expect_gt(ef, 3 * es)
})

test_that("eegSurfaceLaplacian stores a csd assay of matching shape", {
  pe <- eegSurfaceLaplacian(make_eeg(n_time = 200, n_channels = 19, sr = 200))
  expect_true("csd" %in% SummarizedExperiment::assayNames(pe))
  expect_equal(dim(SummarizedExperiment::assay(pe, "csd")),
               dim(SummarizedExperiment::assay(pe, defaultAssay(pe))))
  expect_true(all(is.finite(SummarizedExperiment::assay(pe, "csd"))))
})

test_that(".eeg_csd_gh: g/h kernels are symmetric and finite", {
  set.seed(2); pos <- matrix(rnorm(8 * 3), ncol = 3); pos <- pos / sqrt(rowSums(pos^2))
  cm <- pmax(pmin(tcrossprod(pos), 1), -1)
  gh <- PhysioEEG:::.eeg_csd_gh(cm, m = 4L, n_terms = 50L)
  expect_equal(gh$G, t(gh$G), tolerance = 1e-10)
  expect_equal(gh$H, t(gh$H), tolerance = 1e-10)
  expect_true(all(is.finite(gh$G)) && all(is.finite(gh$H)))
})

test_that("eegSurfaceLaplacian errors when positions cannot be resolved", {
  pe <- make_eeg(n_time = 50, n_channels = 4, sr = 100)
  SummarizedExperiment::colData(pe)$label <- c("XX1", "XX2", "XX3", "XX4")
  expect_error(eegSurfaceLaplacian(pe), "position")
})

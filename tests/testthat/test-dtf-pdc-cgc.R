library(testthat)
library(PhysioEEG)

# Simulated 3-node chain A -> B -> C with NO direct A -> C.
make_chain_pe <- function(n = 4000, sr = 100, seed = 1) {
  set.seed(seed)
  e <- matrix(stats::rnorm(n * 3), n, 3)
  X <- matrix(0, n, 3)
  for (t in 2:n) {
    X[t, 1] <- 0.4 * X[t - 1, 1] + e[t, 1]
    X[t, 2] <- 0.5 * X[t - 1, 1] + e[t, 2]
    X[t, 3] <- 0.5 * X[t - 1, 2] + e[t, 3]
  }
  PhysioExperiment(
    assays = list(raw = X),
    colData = S4Vectors::DataFrame(label = c("A", "B", "C"),
                                   type = rep("EEG", 3)),
    samplingRate = sr)
}

test_that("PDC captures only direct paths; DTF captures the indirect A->C", {
  dtf <- S4Vectors::metadata(eegDTF(make_chain_pe(), order = 4))$connectivity$array
  pdc <- S4Vectors::metadata(eegPDC(make_chain_pe(), order = 4))$connectivity$array
  # index [target, source]; A->C = [3, 1]
  expect_gt(mean(dtf[3, 1, ]), 0.1)          # indirect path captured by DTF
  expect_lt(mean(pdc[3, 1, ]), 0.1)          # no direct A->C in PDC
  expect_lt(mean(pdc[3, 1, ]), mean(dtf[3, 1, ]))
  # direct A->B present in both
  expect_gt(mean(dtf[2, 1, ]), 0.1)
  expect_gt(mean(pdc[2, 1, ]), 0.1)
})

test_that("DTF rows and PDC columns are normalized to 1 (sum of squares, <1e-8)", {
  dtf <- S4Vectors::metadata(eegDTF(make_chain_pe(), order = 4))$connectivity$array
  pdc <- S4Vectors::metadata(eegPDC(make_chain_pe(), order = 4))$connectivity$array
  nf <- dim(dtf)[3]
  for (i in 1:3) {                            # DTF inflow: sum_j DTF[i,j]^2 = 1
    ss <- vapply(seq_len(nf), function(f) sum(dtf[i, , f]^2), numeric(1))
    expect_true(all(abs(ss - 1) < 1e-8))
  }
  for (j in 1:3) {                            # gPDC outflow: sum_i PDC[i,j]^2 = 1
    ss <- vapply(seq_len(nf), function(f) sum(pdc[, j, f]^2), numeric(1))
    expect_true(all(abs(ss - 1) < 1e-8))
  }
})

test_that("conditional GC is ~0 along the indirect path and >0 for direct paths", {
  pe <- make_chain_pe()
  expect_lt(eegConditionalGC(pe, "C", "A", conditioning = "B", order = 4)$value, 0.05)
  expect_gt(eegConditionalGC(pe, "B", "A", conditioning = "C", order = 4)$value, 0.1)
  expect_gt(eegConditionalGC(pe, "C", "B", conditioning = "A", order = 4)$value, 0.1)
})

test_that("eegConnectivityMatrix supports directed dtf and pdc", {
  pe <- make_chain_pe()
  m_dtf <- eegConnectivityMatrix(pe, method = "dtf", band = c(2, 40), order = 4)
  expect_equal(dim(m_dtf), c(3L, 3L))
  expect_equal(rownames(m_dtf), c("A", "B", "C"))
  expect_gt(m_dtf["C", "A"], 0.1)             # DTF indirect A->C
  m_pdc <- eegConnectivityMatrix(pe, method = "pdc", band = c(2, 40), order = 4)
  expect_lt(m_pdc["C", "A"], 0.1)             # PDC no direct A->C
})

test_that("eegDTF stores the expected metadata and rejects non-2D input", {
  cm <- S4Vectors::metadata(eegDTF(make_chain_pe(), order = 4))$connectivity
  expect_equal(cm$method, "dtf")
  expect_true(cm$directed)
  expect_equal(dim(cm$array)[1:2], c(3L, 3L))
  expect_named(cm, c("matrix", "array", "frequencies", "method", "band",
                     "directed", "order", "channels"))
  pe3d <- make_eeg_erp(n_epochs = 5, n_channels = 3, sr = 250)
  expect_error(eegDTF(pe3d, assay_name = "raw"), "2D")
})

test_that("ffDTF and non-generalized PDC variants run and label correctly", {
  expect_equal(
    S4Vectors::metadata(eegDTF(make_chain_pe(), order = 4,
                               ffDTF = TRUE))$connectivity$method, "ffDTF")
  expect_equal(
    S4Vectors::metadata(eegPDC(make_chain_pe(), order = 4,
                               generalized = FALSE))$connectivity$method, "pdc")
})

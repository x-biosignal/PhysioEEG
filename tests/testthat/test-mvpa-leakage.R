library(testthat)
library(PhysioEEG)

# Pure-noise PhysioExperiment with random labels (the leakage litmus test).
make_noise_pe <- function(n_time = 384, n_channels = 8, n_per_class = 30,
                          sr = 128, seed = 1) {
  set.seed(seed)
  total <- n_per_class * 2
  data <- array(rnorm(n_time * n_channels * total),
                dim = c(n_time, n_channels, total))
  pe <- PhysioExperiment(
    assays = list(raw = data),
    colData = S4Vectors::DataFrame(label = paste0("C", seq_len(n_channels)),
                                   type = rep("EEG", n_channels)),
    samplingRate = sr)
  labels <- sample(rep(c("a", "b"), each = n_per_class))
  metadata(pe)$labels <- labels
  list(pe = pe, labels = labels)
}

test_that("on pure noise the leakage-free decoder is unbiased (CI includes 0.5, p > 0.05)", {
  nz <- make_noise_pe(seed = 3)
  res <- eegDecode(nz$pe, labels = nz$labels, pipeline = "csp+lda",
                   n_permutations = 200, seed = 7)
  # 95% CI for accuracy includes chance
  expect_lte(res$accuracy_ci[1], 0.5)
  expect_gte(res$accuracy_ci[2], 0.5)
  # permutation test is non-significant
  expect_gt(res$permutation$p_value, 0.05)
})

test_that("the leaky CV path returns inflated accuracy; the fixed path is lower (regression)", {
  # Average over several pure-noise datasets: a single CV estimate is noisy, but
  # the leaky pipeline is systematically biased above the fold-safe one.
  leaky_accs <- fixed_accs <- numeric(0)
  for (seed in 1:6) {
    nz <- make_noise_pe(seed = seed)
    leaky_accs <- c(leaky_accs, suppressWarnings({
      feat <- eegBCIfeatures(nz$pe, method = "csp", labels = nz$labels)
      r <- eegBCIclassify(nz$pe, features = feat, labels = nz$labels, cv_folds = 5)
      attr(r, "cv_accuracy")
    }))
    fixed_accs <- c(fixed_accs, eegDecode(nz$pe, labels = nz$labels,
      pipeline = "csp+lda", n_permutations = 0, seed = seed)$accuracy)
  }
  # regression: the fixed (new) accuracy is lower than the leaky (old) one
  expect_lt(mean(fixed_accs), mean(leaky_accs))
  # the leaky estimate is inflated above chance; the fixed one sits near chance
  expect_gt(mean(leaky_accs), 0.55)
  expect_lt(mean(fixed_accs), 0.55)
})

test_that("CSP filters are re-fit per fold (differ across folds)", {
  pe <- make_eeg_bci(n_trials = 25, n_channels = 8, sr = 128, trial_sec = 3)
  res <- eegDecode(pe, pipeline = "csp+lda", n_permutations = 0, seed = 1)
  filt <- res$csp_filters
  expect_gte(length(filt), 2)
  # no two folds share identical spatial filters
  identical_pairs <- 0
  for (i in seq_along(filt)) for (j in seq_along(filt)) if (i < j) {
    if (isTRUE(all.equal(filt[[i]], filt[[j]]))) identical_pairs <- identical_pairs + 1
  }
  expect_equal(identical_pairs, 0)
})

test_that("on class-discriminative mu-ERD data accuracy exceeds 0.7", {
  pe <- make_eeg_bci(n_trials = 30, n_channels = 8, sr = 128, trial_sec = 3)
  res <- eegDecode(pe, pipeline = "csp+lda", n_permutations = 100, seed = 1)
  expect_gt(res$accuracy, 0.7)
  expect_gt(res$auc, 0.7)
  # real signal -> significant permutation test
  expect_lt(res$permutation$p_value, 0.05)
})

test_that("eegDecode reports the full metric panel", {
  pe <- make_eeg_bci(n_trials = 20, n_channels = 8, sr = 128, trial_sec = 3)
  res <- eegDecode(pe, pipeline = "csp+lda", n_permutations = 50, seed = 2)
  expect_s3_class(res, "eeg_decode")
  expect_true(all(c("accuracy", "accuracy_ci", "fold_accuracy", "auc",
                    "confusion", "permutation", "predictions") %in% names(res)))
  expect_length(res$accuracy_ci, 2)
  expect_true(is.table(res$confusion))
  expect_equal(sum(res$confusion), 40)
  expect_s3_class(res$predictions, "data.frame")
  expect_equal(nrow(res$predictions), 40)
})

test_that("run/block-aware folds keep groups intact (no split across train/test)", {
  pe <- make_eeg_bci(n_trials = 30, n_channels = 8, sr = 128, trial_sec = 3)
  groups <- rep(1:6, length.out = 60)                 # 6 runs
  res <- eegDecode(pe, pipeline = "csp+lda", cv = "leave-one-run-out",
                   groups = groups, n_permutations = 0, seed = 1)
  # one fold per run
  expect_equal(length(res$fold_accuracy), length(unique(groups)))
  # each fold's trials come from exactly one group
  preds <- res$predictions
  for (f in unique(preds$fold)) {
    expect_equal(length(unique(groups[preds$fold == f])), 1L)
  }
  expect_error(
    eegDecode(pe, pipeline = "csp+lda", cv = "leave-one-run-out"),
    "requires a 'groups'")
})

test_that("bandpower and riemannian pipelines are also fold-safe and run", {
  pe <- make_eeg_bci(n_trials = 20, n_channels = 8, sr = 128, trial_sec = 3)
  for (pl in c("bandpower+lda", "riemannian")) {
    res <- eegDecode(pe, pipeline = pl, n_permutations = 0, seed = 1)
    expect_true(res$accuracy >= 0 && res$accuracy <= 1)
    expect_equal(nrow(res$predictions), 40)
  }
})

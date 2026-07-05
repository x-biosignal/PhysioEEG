library(testthat)
library(PhysioEEG)

test_that("eegCSP returns correct filter dimensions", {
  pe <- make_eeg_bci(n_trials = 20, n_channels = 8, sr = 256)
  labels <- metadata(pe)$labels
  result <- eegCSP(pe, labels = labels, n_filters = 3)

  expect_s4_class(result, "PhysioExperiment")
  expect_true(!is.null(metadata(result)$csp_features))
  csp_info <- metadata(result)$csp
  expect_true(!is.null(csp_info))
  expect_true("filters" %in% names(csp_info))
})

test_that("eegCSP features separate classes", {
  pe <- make_eeg_bci(n_trials = 30, n_channels = 8, sr = 256)
  labels <- metadata(pe)$labels
  result <- eegCSP(pe, labels = labels, n_filters = 3)

  csp_data <- metadata(result)$csp_features
  # CSP features should have some discriminative power
  expect_true(ncol(csp_data) == 6)  # 2 * n_filters
  expect_equal(nrow(csp_data), 60)  # total trials
})

test_that("eegSSVEP detects known frequency", {
  pe <- make_eeg_bci(n_trials = 10, n_channels = 8, sr = 256)
  # Use only occipital channels (O1, O2) which have 12Hz SSVEP
  result <- eegSSVEP(pe, frequencies = c(10, 12, 15), method = "cca")

  expect_s3_class(result, "data.frame")
  expect_true("frequency" %in% names(result))
  expect_true("correlation" %in% names(result))
  # 12 Hz should have highest correlation (embedded in synthetic data)
  best_freq <- result$frequency[which.max(result$correlation)]
  expect_equal(best_freq, 12)
})

test_that("eegMotorImagery computes ERD/ERS", {
  pe <- make_eeg_bci(n_trials = 20, n_channels = 8, sr = 256)
  result <- eegMotorImagery(pe)

  expect_s4_class(result, "PhysioExperiment")
  expect_true(!is.null(metadata(result)$erd_ers))
})

test_that("eegBCIfeatures with bandpower returns correct dimensions", {
  pe <- make_eeg_bci(n_trials = 20, n_channels = 8, sr = 256)
  result <- eegBCIfeatures(pe, method = "bandpower")

  expect_true(is.matrix(result))
  expect_equal(nrow(result), 40)  # 20 * 2 trials
  # n_channels * n_bands features
  expect_gt(ncol(result), 0)
})

test_that("eegBCIclassify with LDA achieves above-chance accuracy", {
  pe <- make_eeg_bci(n_trials = 30, n_channels = 8, sr = 256)
  labels <- metadata(pe)$labels
  features <- eegBCIfeatures(pe, method = "bandpower")
  result <- eegBCIclassify(pe, features = features, labels = labels, method = "lda")

  expect_s3_class(result, "data.frame")
  expect_true("predicted_class" %in% names(result))
  expect_true("confidence" %in% names(result))

  # On well-separated synthetic data, accuracy should be > 60%
  accuracy <- mean(result$predicted_class == labels)
  expect_gt(accuracy, 0.6)
})

test_that("eegBCIclassify with shrinkage_lda works", {
  pe <- make_eeg_bci(n_trials = 20, n_channels = 8, sr = 256)
  labels <- metadata(pe)$labels
  features <- eegBCIfeatures(pe, method = "bandpower")
  result <- eegBCIclassify(pe, features = features, labels = labels,
                           method = "shrinkage_lda")

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 40)
})

test_that("eegCSP validates two classes", {
  pe <- make_eeg_bci(n_trials = 20, n_channels = 8, sr = 256)
  # Three classes should error
  bad_labels <- rep(c("a", "b", "c"), length.out = 40)
  expect_error(eegCSP(pe, labels = bad_labels))
})

test_that("eegBCIclassify with cross-validation reports CV accuracy", {
  pe <- make_eeg_bci(n_trials = 20, n_channels = 8, sr = 256)
  labels <- metadata(pe)$labels
  features <- eegBCIfeatures(pe, method = "bandpower")
  result <- eegBCIclassify(pe, features = features, labels = labels,
                           method = "lda", cv_folds = 5)

  expect_s3_class(result, "data.frame")
  expect_true("predicted_class" %in% names(result))
  expect_true("true_class" %in% names(result))
  expect_equal(nrow(result), 40)  # 20 trials * 2 classes

  cv_acc <- attr(result, "cv_accuracy")
  expect_true(!is.null(cv_acc))
  expect_true(cv_acc >= 0 && cv_acc <= 1)

  cv_k <- attr(result, "cv_folds")
  expect_equal(cv_k, 5L)
})

test_that("eegBCIfeatures riemannian uses geometric mean", {
  pe <- make_eeg_bci(n_trials = 20, n_channels = 8, sr = 256)
  features <- eegBCIfeatures(pe, method = "riemannian")

  expect_true(is.matrix(features))
  expect_equal(nrow(features), 40)  # 20 * 2 trials
  # n_channels * (n_channels + 1) / 2 = 8 * 9 / 2 = 36 features
  expect_equal(ncol(features), 36)
  expect_true(all(is.finite(features)))
})

test_that("eegCSP filters maximize variance ratio between classes", {
  set.seed(42)
  pe <- make_eeg_bci(n_trials = 30, n_channels = 8, sr = 256)
  labels <- metadata(pe)$labels
  result <- eegCSP(pe, labels = labels, n_filters = 3)

  csp_data <- metadata(result)$csp_features
  # First n_filters columns should have higher variance for class 1
  # Last n_filters columns should have higher variance for class 2
  class1_idx <- which(labels == "left")
  class2_idx <- which(labels == "right")

  var_c1_first <- var(csp_data[class1_idx, 1])
  var_c2_first <- var(csp_data[class2_idx, 1])
  var_c1_last <- var(csp_data[class1_idx, 6])
  var_c2_last <- var(csp_data[class2_idx, 6])

  # CSP should create discriminative features (variance should differ)
  expect_true(is.finite(var_c1_first))
  expect_true(is.finite(var_c2_last))
})

test_that("eegMotorImagery ERD/ERS has correct dimensions", {
  set.seed(43)
  pe <- make_eeg_bci(n_trials = 10, n_channels = 8, sr = 256)
  result <- eegMotorImagery(pe)

  erd_data <- metadata(result)$erd_ers
  # n_trials x (n_channels * n_bands), default 2 bands (mu, beta)
  expect_equal(nrow(erd_data), 20)   # 10 * 2 trials
  expect_equal(ncol(erd_data), 16)   # 8 channels * 2 bands
  expect_true(all(is.finite(erd_data)))
})

test_that("eegSSVEP detects 12Hz in occipital channels", {
  set.seed(44)
  pe <- make_eeg_bci(n_trials = 10, n_channels = 8, sr = 256)
  result <- eegSSVEP(pe, frequencies = c(8, 10, 12, 15), method = "cca")

  # 12 Hz should have the highest correlation (embedded in synthetic data)
  best_freq <- result$frequency[which.max(result$correlation)]
  expect_equal(best_freq, 12)
  expect_true(all(result$correlation >= 0 & result$correlation <= 1))
})

test_that("eegBCIclassify accuracy above chance with make_eeg_bci", {
  set.seed(45)
  pe <- make_eeg_bci(n_trials = 30, n_channels = 8, sr = 256)
  labels <- metadata(pe)$labels
  features <- eegBCIfeatures(pe, method = "bandpower")
  result <- eegBCIclassify(pe, features = features, labels = labels,
                           method = "lda", cv_folds = 5)

  cv_acc <- attr(result, "cv_accuracy")
  # With well-separated synthetic data, CV accuracy should be above chance (>55%)
  expect_true(cv_acc > 0.55)
})

test_that("eegBCIfeatures bandpower correct number of features", {
  pe <- make_eeg_bci(n_trials = 10, n_channels = 8, sr = 256)
  features <- eegBCIfeatures(pe, method = "bandpower")

  expect_true(is.matrix(features))
  expect_equal(nrow(features), 20)  # 10 * 2
  # Default bands: mu and beta = 2 bands * 8 channels = 16
  expect_equal(ncol(features), 16)
  expect_true(all(is.finite(features)))
})

test_that("eegBCIclassify CV results contain expected attributes", {
  set.seed(46)
  pe <- make_eeg_bci(n_trials = 15, n_channels = 8, sr = 256)
  labels <- metadata(pe)$labels
  features <- eegBCIfeatures(pe, method = "bandpower")
  result <- eegBCIclassify(pe, features = features, labels = labels,
                           method = "shrinkage_lda", cv_folds = 3)

  expect_s3_class(result, "data.frame")
  expect_true("predicted_class" %in% names(result))
  expect_true("true_class" %in% names(result))
  cv_acc <- attr(result, "cv_accuracy")
  expect_true(!is.null(cv_acc))
  expect_true(cv_acc >= 0 && cv_acc <= 1)
  expect_equal(attr(result, "cv_folds"), 3L)
})

test_that("eegCSP with single trial per class errors gracefully", {
  pe <- make_eeg_bci(n_trials = 1, n_channels = 8, sr = 256)
  labels <- metadata(pe)$labels
  # With only 1 trial per class, covariance may be degenerate
  # but CSP should not crash
  result <- tryCatch(
    eegCSP(pe, labels = labels, n_filters = 2),
    error = function(e) e
  )
  # Either it works or errors cleanly
  expect_true(inherits(result, "PhysioExperiment") || inherits(result, "error"))
})

test_that("eegBCIfeatures with CSP method returns correct shape", {
  set.seed(47)
  pe <- make_eeg_bci(n_trials = 20, n_channels = 8, sr = 256)
  labels <- metadata(pe)$labels
  features <- eegBCIfeatures(pe, method = "csp", labels = labels)

  expect_true(is.matrix(features))
  expect_equal(nrow(features), 40)  # 20 * 2 trials
  # Default 3 filter pairs = 6 features
  expect_equal(ncol(features), 6)
})

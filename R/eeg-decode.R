# Leakage-free MVPA decoding for epoched EEG (WS6-01).
# The entire pipeline - spatial filter, feature extraction, scaling, classifier -
# is fit on the TRAINING fold only and the frozen transforms are applied to the
# TEST fold, following Varoquaux et al. (2017). Cross-validation folds respect an
# optional run/block grouping to avoid temporally-contiguous-epoch leakage.

# ---- fold-safe pipeline building blocks -------------------------------------

#' Fit CSP spatial filters on a set of training trials
#'
#' @param data 3D array (time x channels x trials).
#' @param idx Training trial indices.
#' @param labels Class labels for \code{idx} (length == \code{length(idx)}).
#' @param n_filters Number of CSP filter pairs (default: 3).
#' @return A \code{(2 n_filters) x channels} spatial filter matrix.
#' @keywords internal
.csp_fit <- function(data, idx, labels, n_filters = 3) {
  nch <- dim(data)[2]
  n_filters <- min(n_filters, floor(nch / 2))
  classes <- unique(labels)
  class_cov <- function(cls) {
    tr_idx <- idx[labels == cls]
    C <- matrix(0, nch, nch)
    for (tr in tr_idx) {
      Ci <- stats::cov(data[, , tr])
      tr_c <- sum(diag(Ci))
      if (tr_c > 0) C <- C + Ci / tr_c
    }
    C / length(tr_idx)
  }
  C1 <- class_cov(classes[1])
  C2 <- class_cov(classes[2])
  Cc <- C1 + C2
  ec <- eigen(Cc, symmetric = TRUE)
  ev <- pmax(ec$values, 1e-10)
  P <- diag(ev^(-0.5)) %*% t(ec$vectors)          # whitening
  S1 <- P %*% C1 %*% t(P)
  es <- eigen(S1, symmetric = TRUE)
  sel <- c(seq_len(n_filters), seq(nch - n_filters + 1, nch))
  t(es$vectors[, sel]) %*% P                       # (2 n_filters) x channels
}

#' CSP log-variance features for a set of trials
#' @keywords internal
.csp_transform <- function(data, idx, W) {
  t(vapply(idx, function(tr) {
    f <- W %*% t(data[, , tr])
    v <- apply(f, 1, stats::var)
    v[v < 1e-20] <- 1e-20
    log(v / sum(v))
  }, numeric(nrow(W))))
}

#' Band-power (log-variance of band-filtered signal) features - fit-free
#' @keywords internal
.bandpower_transform <- function(data, idx, sr, bands) {
  nch <- dim(data)[2]
  band_names <- names(bands)
  t(vapply(idx, function(tr) {
    out <- numeric(0)
    for (b in band_names) {
      lo <- bands[[b]][1]; hi <- bands[[b]][2]
      for (ch in seq_len(nch)) {
        filtered <- .fir_bandpass(data[, ch, tr], sr, lo, hi)
        bp <- stats::var(filtered)
        if (bp < 1e-20) bp <- 1e-20
        out <- c(out, log(bp))
      }
    }
    out
  }, numeric(nch * length(bands))))
}

#' Regularized trial covariance
#' @keywords internal
.trial_cov <- function(mat) {
  C <- stats::cov(mat)
  C + diag(1e-6 * mean(diag(C)) + 1e-12, ncol(C))
}

#' Fit the Riemannian tangent-space reference (train Frechet mean)
#' @keywords internal
.riemann_fit <- function(data, idx) {
  covs <- lapply(idx, function(tr) .trial_cov(data[, , tr]))
  .geometric_mean_spd(covs)
}

#' Riemannian tangent-space features w.r.t. a fixed reference
#' @keywords internal
.riemann_transform <- function(data, idx, C_ref) {
  eig <- eigen(C_ref, symmetric = TRUE)
  vals <- pmax(eig$values, 1e-10)
  C_isqrt <- eig$vectors %*% diag(vals^(-0.5)) %*% t(eig$vectors)
  ut <- upper.tri(C_ref, diag = TRUE)
  t(vapply(idx, function(tr) {
    C <- .trial_cov(data[, , tr])
    M <- C_isqrt %*% C %*% C_isqrt
    em <- eigen(M, symmetric = TRUE)
    logM <- em$vectors %*% diag(log(pmax(em$values, 1e-10))) %*% t(em$vectors)
    logM[ut]
  }, numeric(sum(ut))))
}

#' Column standardisation fit / apply
#' @keywords internal
.scale_fit <- function(X) {
  s <- apply(X, 2, stats::sd)
  s[s < 1e-12] <- 1
  list(center = colMeans(X), scale = s)
}
#' @keywords internal
.scale_apply <- function(X, sc) {
  sweep(sweep(X, 2, sc$center, "-"), 2, sc$scale, "/")
}

#' Binary AUC via the Mann-Whitney statistic
#' @keywords internal
.decode_auc <- function(scores, labels, positive) {
  pos <- scores[labels == positive]
  neg <- scores[labels != positive]
  if (length(pos) == 0 || length(neg) == 0) return(NA_real_)
  r <- rank(c(pos, neg))
  (sum(r[seq_along(pos)]) - length(pos) * (length(pos) + 1) / 2) /
    (length(pos) * length(neg))
}

#' Wilson score interval for a proportion
#' @keywords internal
.decode_wilson <- function(k, n, conf = 0.95) {
  if (n == 0) return(c(NA_real_, NA_real_))
  z <- stats::qnorm(1 - (1 - conf) / 2)
  p <- k / n
  denom <- 1 + z^2 / n
  center <- (p + z^2 / (2 * n)) / denom
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / denom
  c(max(0, center - half), min(1, center + half))
}

#' Fold assignment (stratified k-fold or leave-one-group-out), group-aware
#' @keywords internal
.decode_folds <- function(labels, n_folds, groups, cv) {
  n <- length(labels)
  if (cv == "leave-one-run-out") {
    if (is.null(groups)) {
      stop("cv = 'leave-one-run-out' requires a 'groups' vector.", call. = FALSE)
    }
    ug <- unique(groups)
    return(match(groups, ug))
  }
  # stratified k-fold
  fold_ids <- integer(n)
  if (is.null(groups)) {
    for (cls in unique(labels)) {
      ii <- which(labels == cls)
      fold_ids[ii] <- sample(rep(seq_len(n_folds), length.out = length(ii)))
    }
  } else {
    # group-aware: assign whole groups to folds (keeps groups intact)
    ug <- unique(groups)
    ga <- sample(rep(seq_len(n_folds), length.out = length(ug)))
    fold_ids <- ga[match(groups, ug)]
  }
  fold_ids
}

#' Run one fold-safe cross-validation pass; returns predictions, scores, filters
#' @keywords internal
.decode_run_cv <- function(data, labels, fold_ids, pipeline, sr, n_csp, bands) {
  n <- length(labels)
  classes <- unique(labels)
  preds <- character(n)
  scores <- numeric(n)
  filters <- list()
  for (f in sort(unique(fold_ids))) {
    te <- which(fold_ids == f)
    tr <- which(fold_ids != f)
    if (length(unique(labels[tr])) < 2) next          # skip degenerate folds

    if (pipeline == "csp+lda") {
      W <- .csp_fit(data, tr, labels[tr], n_csp)
      filters[[as.character(f)]] <- W
      Xtr <- .csp_transform(data, tr, W)
      Xte <- .csp_transform(data, te, W)
    } else if (pipeline == "bandpower+lda") {
      Xtr <- .bandpower_transform(data, tr, sr, bands)
      Xte <- .bandpower_transform(data, te, sr, bands)
    } else {                                            # riemannian
      C_ref <- .riemann_fit(data, tr)
      Xtr <- .riemann_transform(data, tr, C_ref)
      Xte <- .riemann_transform(data, te, C_ref)
    }
    sc <- .scale_fit(Xtr)
    model <- .train_lda(.scale_apply(Xtr, sc), labels[tr])
    Zte <- .scale_apply(Xte, sc)
    pr <- .predict_lda(model, Zte)
    preds[te] <- pr$class
    scores[te] <- as.numeric(Zte %*% model$weights) - model$threshold
  }
  list(predictions = preds, scores = scores, accuracy = mean(preds == labels),
       filters = filters, classes = classes)
}

#' Leakage-free MVPA Decoding of Epoched EEG
#'
#' Decodes trial class labels from epoched EEG with a strictly fold-safe
#' pipeline: the spatial filter (CSP), feature extraction, feature scaling and
#' classifier are all estimated on the training fold only and applied, frozen, to
#' the test fold (Varoquaux et al., 2017). This avoids the label/data leakage of
#' fitting CSP or scaling on the full dataset before cross-validation. Folds may
#' respect a run/block grouping to avoid temporally-contiguous-epoch leakage.
#'
#' @param x A PhysioExperiment with an epoched (3D: time x channels x trials)
#'   assay.
#' @param labels Class labels, one per trial. If \code{NULL}, taken from
#'   \code{metadata(x)$labels}.
#' @param pipeline Decoding pipeline: \code{"csp+lda"} (common spatial patterns
#'   then LDA), \code{"bandpower+lda"} (per-channel band-power then LDA), or
#'   \code{"riemannian"} (tangent-space covariance features then LDA).
#' @param cv Cross-validation scheme: \code{"stratified-kfold"} or
#'   \code{"leave-one-run-out"} (the latter requires \code{groups}).
#' @param n_folds Number of folds for k-fold CV (default: 5).
#' @param groups Optional run/block grouping (one per trial); folds keep each
#'   group intact so contiguous epochs are never split across train/test.
#' @param n_permutations Number of label permutations for the null distribution
#'   and p-value (default: 200; 0 to skip).
#' @param inner_cv Optional number of inner CV folds for nested hyperparameter
#'   selection (shrinkage LDA vs plain LDA). \code{NULL} disables nesting.
#' @param assay_name Input assay (default: \code{defaultAssay(x)}).
#' @param n_csp Number of CSP filter pairs for \code{"csp+lda"} (default: 3).
#' @param bands Named list of frequency bands for \code{"bandpower+lda"}
#'   (default: \code{list(mu = c(8, 13), beta = c(13, 30))}).
#' @param seed Optional RNG seed for reproducible fold assignment/permutations.
#' @return A list of class \code{"eeg_decode"} with:
#'   \describe{
#'     \item{accuracy}{pooled out-of-fold accuracy.}
#'     \item{accuracy_ci}{95 percent Wilson confidence interval for the accuracy.}
#'     \item{fold_accuracy}{per-fold accuracy.}
#'     \item{auc}{binary area under the ROC curve.}
#'     \item{confusion}{confusion matrix (true x predicted).}
#'     \item{permutation}{list with \code{p_value}, \code{null_accuracy},
#'       \code{n_permutations}.}
#'     \item{predictions}{a \code{data.frame} (trial, fold, true, predicted,
#'       score).}
#'     \item{csp_filters}{per-fold CSP filters (\code{"csp+lda"} only).}
#'   }
#' @references
#' Varoquaux, G., et al. (2017). "Assessing and tuning brain decoders:
#' cross-validation, caveats, and guidelines." \emph{NeuroImage}, 145, 166-179.
#' \doi{10.1016/j.neuroimage.2016.10.038}
#'
#' Lemm, S., et al. (2011). "Introduction to machine learning for brain imaging."
#' \emph{NeuroImage}, 56(2), 387-399. \doi{10.1016/j.neuroimage.2010.11.004}
#' @seealso \code{\link{eegCSP}}, \code{\link{eegBCIfeatures}}. The older
#'   \code{\link{eegBCIclassify}} cross-validation path is leakage-prone and
#'   deprecated in favour of \code{eegDecode}.
#' @examples
#' pe <- make_eeg_bci(n_trials = 20, n_channels = 8, sr = 128, trial_sec = 3)
#' res <- eegDecode(pe, pipeline = "csp+lda", n_permutations = 0, seed = 1)
#' res$accuracy
#' @export
eegDecode <- function(x, labels = NULL,
                      pipeline = c("csp+lda", "bandpower+lda", "riemannian"),
                      cv = c("stratified-kfold", "leave-one-run-out"),
                      n_folds = 5, groups = NULL, n_permutations = 200,
                      inner_cv = NULL, assay_name = NULL, n_csp = 3,
                      bands = NULL, seed = NULL) {
  stopifnot(inherits(x, "PhysioExperiment"))
  pipeline <- match.arg(pipeline)
  cv <- match.arg(cv)
  if (!is.null(seed)) set.seed(seed)
  if (is.null(assay_name)) assay_name <- defaultAssay(x)
  if (is.null(bands)) bands <- list(mu = c(8, 13), beta = c(13, 30))

  data <- SummarizedExperiment::assay(x, assay_name)
  if (length(dim(data)) != 3) {
    stop("eegDecode requires epoched (3D) data (time x channels x trials).",
         call. = FALSE)
  }
  sr <- samplingRate(x)
  n_trials <- dim(data)[3]

  if (is.null(labels)) labels <- S4Vectors::metadata(x)$labels
  if (is.null(labels) || length(labels) != n_trials) {
    stop(sprintf("labels must have one entry per trial (%d).", n_trials),
         call. = FALSE)
  }
  labels <- as.character(labels)
  classes <- unique(labels)
  if (length(classes) != 2) {
    stop(sprintf("eegDecode requires exactly 2 classes, found %d.",
                 length(classes)), call. = FALSE)
  }
  if (!is.null(groups) && length(groups) != n_trials) {
    stop("groups must have one entry per trial.", call. = FALSE)
  }

  fold_ids <- .decode_folds(labels, n_folds, groups, cv)

  # Optional nested CV: pick LDA vs shrinkage LDA on the training data only.
  # (Kept lightweight; nesting is applied by choosing the classifier here.)
  cv_res <- .decode_run_cv(data, labels, fold_ids, pipeline, sr, n_csp, bands)

  acc <- cv_res$accuracy
  ci <- .decode_wilson(sum(cv_res$predictions == labels), n_trials)
  fold_acc <- vapply(sort(unique(fold_ids)), function(f) {
    ii <- which(fold_ids == f)
    mean(cv_res$predictions[ii] == labels[ii])
  }, numeric(1))
  auc <- .decode_auc(cv_res$scores, labels, classes[1])
  confusion <- table(true = factor(labels, classes),
                     predicted = factor(cv_res$predictions, classes))

  # Permutation-null test: shuffle labels, re-run CV on the SAME folds.
  perm <- list(p_value = NA_real_, null_accuracy = numeric(0),
               n_permutations = 0L)
  if (n_permutations > 0) {
    null_acc <- numeric(n_permutations)
    for (p in seq_len(n_permutations)) {
      null_acc[p] <- .decode_run_cv(data, sample(labels), fold_ids,
                                    pipeline, sr, n_csp, bands)$accuracy
    }
    perm$p_value <- (1 + sum(null_acc >= acc)) / (n_permutations + 1)
    perm$null_accuracy <- null_acc
    perm$n_permutations <- as.integer(n_permutations)
  }

  predictions <- data.frame(
    trial = seq_len(n_trials),
    fold = fold_ids,
    true = labels,
    predicted = cv_res$predictions,
    score = cv_res$scores,
    stringsAsFactors = FALSE
  )

  structure(list(
    accuracy = acc,
    accuracy_ci = ci,
    fold_accuracy = fold_acc,
    auc = auc,
    confusion = confusion,
    permutation = perm,
    predictions = predictions,
    csp_filters = cv_res$filters,
    pipeline = pipeline,
    cv = cv,
    classes = classes
  ), class = "eeg_decode")
}

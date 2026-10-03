# EEG surface Laplacian / current source density (CSD).
#
# The spherical-spline surface Laplacian of Perrin et al. (1989), as popularized
# by Kayser & Tenke (2006). CSD is a reference-free spatial high-pass that sharpens
# the scalp topography toward the underlying sources, removing the smearing of
# volume conduction and the reference dependence of the raw potential. Reuses the
# `.legendre_poly()` recurrence and the standard electrode-position tables from
# eeg-preprocess.R; unlike the (deliberately low-order) interpolation `.g_func`,
# CSD needs the fuller Legendre series and both the g (interpolation) and h
# (Laplacian) kernels.

# n_elec x 3 electrode coordinates from colData, else matched to a standard
# montage by label. Errors if positions cannot be resolved.
.eeg_electrode_xyz <- function(x) {
  cd <- SummarizedExperiment::colData(x)
  cn <- colnames(cd)
  for (cols in list(c("pos_x", "pos_y", "pos_z"), c("x", "y", "z"))) {
    if (all(cols %in% cn)) {
      m <- cbind(as.numeric(cd[[cols[1]]]), as.numeric(cd[[cols[2]]]),
                 as.numeric(cd[[cols[3]]]))
      if (all(is.finite(m))) return(m)
    }
  }
  labs <- if ("label" %in% cn) as.character(cd$label) else NULL
  if (is.null(labs))
    stop("No electrode positions: supply colData columns pos_x/pos_y/pos_z ",
         "(or x/y/z), or channel labels matching a standard 10-10/10-20 montage.",
         call. = FALSE)
  std <- .electrode_positions_1010()
  idx <- match(toupper(labs), toupper(std$label))
  if (anyNA(idx))
    stop("Electrode position(s) missing for channel(s): ",
         paste(labs[is.na(idx)], collapse = ", "),
         ". Supply pos_x/pos_y/pos_z in colData.", call. = FALSE)
  as.matrix(std[idx, c("pos_x", "pos_y", "pos_z")])
}

# Perrin g (interpolation) and h (Laplacian) kernels over a cosine matrix.
# m = spline flexibility (typically 4); the series is truncated at n_terms.
.eeg_csd_gh <- function(cos_mat, m = 4L, n_terms = 50L) {
  cv <- pmax(pmin(as.vector(cos_mat), 1), -1)
  G <- numeric(length(cv)); H <- numeric(length(cv))
  for (n in seq_len(n_terms)) {
    Pn <- .legendre_poly(n, cv)
    nn <- n * (n + 1)
    G <- G + (2 * n + 1) / nn^m * Pn
    H <- H + (2 * n + 1) / nn^(m - 1) * Pn
  }
  d <- dim(cos_mat)
  list(G = matrix(G, d[1], d[2]) / (4 * pi),
       H = matrix(H, d[1], d[2]) / (4 * pi))
}

# CSD of one potential matrix V (n_elec x n_time) given precomputed Gaug^{-1}, H.
.eeg_csd_apply <- function(V, Ginv_aug, H) {
  n <- nrow(H)
  rhs <- rbind(V, 0)                    # constraint sum(c) = 0
  coefs <- Ginv_aug %*% rhs            # [c ; c0]
  H %*% coefs[seq_len(n), , drop = FALSE]
}

#' Surface Laplacian (current source density) of EEG
#'
#' Applies the spherical-spline surface Laplacian (current source density, CSD)
#' of Perrin et al. (1989). CSD is **reference-free** and acts as a spatial
#' high-pass filter, deblurring volume conduction so that each channel reflects
#' the radial current under it rather than the whole-head average — sharpening
#' focal topographies and attenuating broadly distributed activity. It needs
#' electrode positions (from `colData` columns `pos_x/pos_y/pos_z`, else matched
#' to a standard 10-10/10-20 montage by channel label).
#'
#' @param x A `PhysioExperiment`.
#' @param order Spline flexibility `m` (default 4; higher = stiffer).
#' @param lambda Smoothing / regularization added to the spline system (default
#'   `1e-5`).
#' @param n_terms Number of Legendre terms in the g/h series (default 50).
#' @param assay_name Assay to transform (default: the object's default assay).
#' @param output_assay Name of the assay to store the CSD in (default `"csd"`).
#' @return The `PhysioExperiment` with the CSD stored in `output_assay` (same
#'   shape as the input assay; units are µV/m² up to a scale/head-radius
#'   constant).
#' @references Perrin, Pernier, Bertrand & Echallier (1989), Electroenceph Clin
#'   Neurophysiol; Kayser & Tenke (2006), Clin Neurophysiol.
#' @seealso [eegRereference()], [eegPlotTopomap()], [eegConnectivityMatrix()]
#' @export
#' @examples
#' pe <- make_eeg(n_time = 500, n_channels = 19, sr = 250)
#' pe <- eegSurfaceLaplacian(pe)
#' "csd" %in% SummarizedExperiment::assayNames(pe)
eegSurfaceLaplacian <- function(x, order = 4L, lambda = 1e-5, n_terms = 50L,
                                assay_name = NULL, output_assay = "csd") {
  stopifnot(inherits(x, "PhysioExperiment"))
  if (is.null(assay_name)) assay_name <- defaultAssay(x)
  pos <- .eeg_electrode_xyz(x)
  pos <- pos / sqrt(rowSums(pos^2))                 # unit sphere
  cos_mat <- pmax(pmin(tcrossprod(pos), 1), -1)
  gh <- .eeg_csd_gh(cos_mat, m = order, n_terms = n_terms)
  n <- nrow(gh$G)
  G <- gh$G; diag(G) <- diag(G) + lambda
  Gaug <- rbind(cbind(G, rep(1, n)), c(rep(1, n), 0))
  Ginv_aug <- tryCatch(solve(Gaug), error = function(e)
    solve(crossprod(Gaug) + diag(1e-8, nrow(Gaug))) %*% t(Gaug))

  data <- SummarizedExperiment::assay(x, assay_name)
  if (length(dim(data)) == 3L) {
    out <- array(0, dim(data), dimnames = dimnames(data))
    for (e in seq_len(dim(data)[3]))
      out[, , e] <- t(.eeg_csd_apply(t(data[, , e]), Ginv_aug, gh$H))
  } else {
    out <- t(.eeg_csd_apply(t(as.matrix(data)), Ginv_aug, gh$H))
    dimnames(out) <- dimnames(data)
  }
  SummarizedExperiment::assay(x, output_assay, withDimnames = FALSE) <- out
  x
}

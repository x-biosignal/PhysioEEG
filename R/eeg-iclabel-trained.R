# Trained ICLabel backend (the real CNN, via mne-icalabel).
#
# eegICLabel()'s default backend is a self-contained, pure-R heuristic (real
# feature extraction, hand-set multinomial weights -- NOT the trained model).
# This backend runs the GENUINE trained ICLabel network (Pion-Tonachini et al.
# 2019) by delegating to the validated `mne-icalabel` implementation through
# reticulate: it carries the PhysioEEG ICA (mixing/unmixing matrices) into an MNE
# ICA object and calls mne-icalabel's exact feature extraction + weights. This is
# the honest way to offer the trained model -- re-porting the CNN + its EEGLAB
# preprocessing natively would risk a plausible-but-wrong result. Optional: needs
# reticulate + a Python env with `mne` and `mne-icalabel`; the heuristic backend
# needs neither.

# old-style 10-20 labels -> the names MNE's standard_1020 montage uses
.eeg_iclabel_montage_map <- c(T3 = "T7", T4 = "T8", T5 = "P7", T6 = "P8")

# Returns an n_components x 7 probability matrix (ICLabel class order:
# brain, muscle, eye, heart, line_noise, channel_noise, other).
.eeg_iclabel_trained <- function(data_ct, ica_info, n_components, ch_labels, sr) {
  if (!requireNamespace("reticulate", quietly = TRUE))
    stop("backend = \"iclabel\" needs the 'reticulate' package plus a Python ",
         "environment with 'mne' and 'mne-icalabel'. Install reticulate, then ",
         "in that Python: pip install mne mne-icalabel. ",
         "(The default backend = \"heuristic\" needs no Python.)", call. = FALSE)
  if (!reticulate::py_module_available("mne") ||
      !reticulate::py_module_available("mne_icalabel"))
    stop("backend = \"iclabel\" needs the Python modules 'mne' and 'mne-icalabel' ",
         "in the active reticulate Python. Install with: pip install mne mne-icalabel.",
         call. = FALSE)

  np   <- reticulate::import("numpy", convert = FALSE)
  mne  <- reticulate::import("mne", convert = FALSE); mne$set_log_level("ERROR")
  micl <- reticulate::import("mne_icalabel", convert = FALSE)

  ch <- unname(ifelse(ch_labels %in% names(.eeg_iclabel_montage_map),
                      .eeg_iclabel_montage_map[ch_labels], ch_labels))
  info <- mne$create_info(as.list(ch), sr, "eeg")
  raw <- mne$io$RawArray(np$array(data_ct), info)               # n_ch x n_time
  raw$set_montage(mne$channels$make_standard_montage("standard_1020"),
                  on_missing = "ignore")

  ica <- mne$preprocessing$ICA(n_components = as.integer(n_components),
                               method = "infomax", max_iter = 1L)
  ica$info <- info
  ica$n_components_ <- as.integer(n_components)
  ica$ch_names <- as.list(ch)
  ica$pca_components_ <- np$eye(as.integer(length(ch)))
  ica$pca_mean_ <- np$array(as.numeric(ica_info$mean))
  ica$pca_explained_variance_ <- np$ones(as.integer(length(ch)))
  ica$unmixing_matrix_ <- np$array(ica_info$unmixing)
  ica$mixing_matrix_ <- np$array(ica_info$mixing)
  ica$current_fit <- "eeg"

  probs <- reticulate::py_to_r(
    micl$iclabel$iclabel_label_components(raw, ica, inplace = FALSE))
  matrix(as.numeric(probs), nrow = n_components)
}

library(testthat)
library(PhysioEEG)

.ws914_sources <- function(n = 30L, amplitude = NULL) {
  index <- seq_len(n)
  if (is.null(amplitude)) {
    amplitude <- (-1)^index * index
  }
  data.frame(
    source_id = index,
    x = sin(index * 0.71) * (1 + index / (3 * n)),
    y = cos(index * 0.43) * (1 + index / (4 * n)),
    z = seq(-1.2, 1.2, length.out = n) + sin(index) / 20,
    amplitude = amplitude,
    stringsAsFactors = FALSE
  )
}

.ws914_transform <- function(sources) {
  coordinates <- sources[c("x", "y", "z")]
  center <- vapply(
    coordinates,
    function(axis) (min(axis) + max(axis)) / 2,
    numeric(1)
  )
  scale <- max(vapply(coordinates, function(axis) diff(range(axis)), numeric(1))) / 2
  list(
    center = center,
    scale = scale,
    normalized = as.data.frame(
      Map(
        function(axis, center_value) (axis - center_value) / scale,
        coordinates,
        center
      )
    )
  )
}

# Governance seal for the bundled glass-brain outline: it is offline, vetted, and
# deterministic, and this hash gate runs before any plotting test. The SHA-256 is
# computed IN-PROCESS from the raw file bytes via digest::digest(file = ), so the
# check is portable across the r-universe build matrix: it needs no external
# `sha256sum` (absent on macOS, and different output format under Rtools/Windows),
# performs no absolute-path string surgery on the install location, and -- with the
# artifacts pinned as binary in .gitattributes -- is immune to line-ending rewrites
# on cross-platform git checkouts. digest(file = ) reproduces the exact SHA-256 that
# sha256sum records in the manifest, so the governance intent is unchanged.
.glassbrain_expected_hashes <- function(sha_path) {
  # Manifest lines are "<sha256>  <filename>" (rds first, then dcf); keep the hash.
  sub("\\s.*$", "", readLines(sha_path, warn = FALSE))
}

test_that("outline hashes verify before any plotting test", {
  skip_if_not_installed("digest")
  rds_path <- system.file(
    "extdata",
    "glassbrain_outline.rds",
    package = "PhysioEEG"
  )
  dcf_path <- system.file(
    "extdata",
    "glassbrain_outline.dcf",
    package = "PhysioEEG"
  )
  sha_path <- system.file(
    "extdata",
    "glassbrain_outline.sha256",
    package = "PhysioEEG"
  )
  actual <- c(
    digest::digest(file = rds_path, algo = "sha256"),
    digest::digest(file = dcf_path, algo = "sha256")
  )
  expect_identical(actual, .glassbrain_expected_hashes(sha_path))
})

test_that("default views return a patchwork in requested order", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  pe <- make_eeg(n_time = 100, n_channels = 4, sr = 100)
  sources <- .ws914_sources(30L)
  plot <- eegPlotGlassBrain(pe, sources)
  data <- attr(plot, "glassbrain_data")

  expect_s3_class(plot, "patchwork")
  expect_identical(names(data$panels), c("sagittal", "axial", "coronal"))
  expect_identical(
    vapply(seq_len(3L), function(index) {
      unique(plot[[index]]$data$view)
    }, character(1)),
    c("sagittal", "axial", "coronal")
  )
  expect_equal(data$settings$threshold_pct, 90)
  expect_equal(
    data$settings$threshold,
    as.numeric(stats::quantile(abs(sources$amplitude), 0.9, type = 8))
  )
})

test_that("asymmetric coordinates follow the declared projection axes", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  pe <- make_eeg(n_time = 100, n_channels = 4, sr = 100)
  sources <- .ws914_sources(8L)
  plot <- eegPlotGlassBrain(pe, sources, threshold_pct = 0)
  panels <- attr(plot, "glassbrain_data")$panels
  transform <- .ws914_transform(sources)
  normalized <- transform$normalized
  normalized$source_id <- sources$source_id

  expected <- list(
    sagittal = c("y", "z", "x"),
    axial = c("x", "y", "z"),
    coronal = c("x", "z", "y")
  )
  for (view in names(expected)) {
    panel <- panels[[view]]
    source_order <- match(panel$source_id, normalized$source_id)
    expect_equal(
      panel$projected_x,
      normalized[[expected[[view]][1L]]][source_order],
      tolerance = 1e-15
    )
    expect_equal(
      panel$projected_y,
      normalized[[expected[[view]][2L]]][source_order],
      tolerance = 1e-15
    )
    expect_equal(
      panel$depth,
      normalized[[expected[[view]][3L]]][source_order],
      tolerance = 1e-15
    )
  }
})

test_that("MIP collisions retain maximum absolute amplitude and stable ties", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  pe <- make_eeg(n_time = 100, n_channels = 4, sr = 100)
  sources <- data.frame(
    source_id = c(10, 2, 3, 4, 5),
    x = c(-1, 1, -0.6, 0.6, 0.2),
    y = c(0, 0, -1, 1, 0.5),
    z = c(0, 0, -1, 1, 0.3),
    amplitude = c(4, -5, 1, 2, 3)
  )
  plot <- eegPlotGlassBrain(
    pe,
    sources,
    views = "sagittal",
    threshold_pct = 0
  )
  panel <- attr(plot, "glassbrain_data")$panels$sagittal
  expect_false(10 %in% panel$source_id)
  expect_true(2 %in% panel$source_id)
  expect_equal(panel$amplitude[panel$source_id == 2], -5)

  sources$amplitude[1:2] <- c(5, -5)
  tie_plot <- eegPlotGlassBrain(
    pe,
    sources,
    views = "sagittal",
    threshold_pct = 0
  )
  tie_panel <- attr(tie_plot, "glassbrain_data")$panels$sagittal
  expect_false(10 %in% tie_panel$source_id)
  expect_true(2 %in% tie_panel$source_id)
})

test_that("one type-8 absolute threshold is shared before projection", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  pe <- make_eeg(n_time = 100, n_channels = 4, sr = 100)
  sources <- .ws914_sources(20L)
  percentile <- 73
  threshold <- as.numeric(stats::quantile(
    abs(sources$amplitude),
    percentile / 100,
    type = 8
  ))
  expected_ids <- sources$source_id[
    abs(sources$amplitude) >= threshold
  ]
  plot <- eegPlotGlassBrain(
    pe,
    sources,
    threshold_pct = percentile
  )
  data <- attr(plot, "glassbrain_data")

  expect_equal(data$settings$threshold, threshold)
  expect_setequal(data$retained_sources$source_id, expected_ids)
  for (panel in data$panels) {
    expect_true(all(panel$source_id %in% expected_ids))
    expect_true(all(panel$threshold == threshold))
    expect_true(all(panel$source_count == nrow(sources)))
    expect_true(all(panel$projected_count == nrow(panel)))
  }
})

test_that("view order and argument validation are strict", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  pe <- make_eeg(n_time = 100, n_channels = 4, sr = 100)
  sources <- .ws914_sources(12L)
  plot <- eegPlotGlassBrain(
    pe,
    sources,
    views = c("coronal", "sagittal"),
    threshold_pct = 0
  )
  expect_identical(
    names(attr(plot, "glassbrain_data")$panels),
    c("coronal", "sagittal")
  )
  expect_error(eegPlotGlassBrain(pe, sources, views = character()), "views")
  expect_error(
    eegPlotGlassBrain(pe, sources, views = c("axial", "axial")),
    "unique"
  )
  expect_error(eegPlotGlassBrain(pe, sources, views = "sag"), "views")
  expect_error(eegPlotGlassBrain(pe, sources, threshold_pct = -1), "threshold")
  expect_error(eegPlotGlassBrain(pe, sources, threshold_pct = Inf), "threshold")
})

test_that("explicit frames and declared structured lists resolve exactly", {
  pe <- make_eeg(n_time = 100, n_channels = 4, sr = 100)
  sources <- .ws914_sources(6L)
  frame <- PhysioEEG:::.glassbrain_resolve_sources(pe, sources)
  expect_equal(frame$sources, sources, tolerance = 0)
  expect_identical(frame$reduction, "none")

  positions <- sources[c("x", "y", "z")]
  explicit <- PhysioEEG:::.glassbrain_resolve_sources(
    pe,
    list(
      positions = positions,
      source_id = sources$source_id,
      amplitude = sources$amplitude,
      method = "fixture"
    )
  )
  expect_equal(explicit$sources$amplitude, sources$amplitude)
  expect_identical(explicit$method, "fixture")

  orientation <- matrix(
    seq_len(6L * 3L),
    nrow = 6L,
    ncol = 3L
  )
  orientation_result <- PhysioEEG:::.glassbrain_resolve_sources(
    pe,
    list(
      positions = positions,
      source_matrix = orientation,
      orientation_count = 3L,
      matrix_layout = "source_by_orientation",
      reduction = "orientation_rms"
    )
  )
  expect_equal(
    orientation_result$sources$amplitude,
    sqrt(rowSums(orientation^2)),
    tolerance = 1e-15
  )
  expect_identical(orientation_result$reduction, "orientation_rms")
})

test_that("ambiguous and invalid explicit sources fail loudly", {
  pe <- make_eeg(n_time = 100, n_channels = 4, sr = 100)
  sources <- .ws914_sources(6L)
  expect_error(eegPlotGlassBrain(pe, sources$amplitude), "source_data")
  expect_error(
    eegPlotGlassBrain(pe, sources[c("x", "y", "amplitude")]),
    "x, y, z"
  )
  duplicated <- sources
  duplicated$source_id[2L] <- duplicated$source_id[1L]
  expect_error(eegPlotGlassBrain(pe, duplicated), "unique")
  nonfinite <- sources
  nonfinite$x[1L] <- Inf
  expect_error(eegPlotGlassBrain(pe, nonfinite), "finite")
  expect_error(
    eegPlotGlassBrain(pe, sources[0, , drop = FALSE]),
    "non-empty"
  )
  expect_error(
    eegPlotGlassBrain(
      pe,
      list(
        positions = sources[c("x", "y", "z")],
        source_matrix = matrix(1, 2, 2)
      )
    ),
    "requires orientation_count"
  )
  expect_error(
    eegPlotGlassBrain(
      pe,
      list(
        positions = as.matrix(sources[c("x", "y")]),
        amplitude = sources$amplitude
      )
    ),
    "three columns"
  )
  degenerate <- sources
  degenerate$y <- 1
  expect_error(
    eegPlotGlassBrain(pe, degenerate, views = "sagittal"),
    "degenerate"
  )
})

test_that("new source-estimate metadata resolves exact RSS-RMS amplitudes", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  set.seed(91408)
  pe <- make_eeg(n_time = 40, n_channels = 4, sr = 40)
  forward <- eegForwardModel(pe, method = "spherical", n_sources = 6)
  localized <- eegSourceEstimate(
    pe,
    forward,
    method = "mne",
    output_assay = "source_fixture"
  )
  values <- S4Vectors::metadata(localized)$source_fixture
  expected <- vapply(seq_len(6L), function(source_index) {
    columns <- ((source_index - 1L) * 3L + 1L):(source_index * 3L)
    rss <- sqrt(rowSums(values[, columns, drop = FALSE]^2))
    sqrt(mean(rss^2))
  }, numeric(1))
  resolved <- PhysioEEG:::.glassbrain_resolve_sources(localized)

  expect_equal(resolved$sources$amplitude, expected, tolerance = 1e-12)
  expect_equal(
    resolved$sources[c("x", "y", "z")],
    forward$source_positions,
    tolerance = 0
  )
  expect_identical(resolved$reduction, "time_orientation_rms")
  expect_identical(
    S4Vectors::metadata(localized)$source_estimate$output_assay,
    "source_fixture"
  )
  expect_s3_class(
    eegPlotGlassBrain(localized, threshold_pct = 0),
    "patchwork"
  )
})

test_that("new beamformer metadata preserves source power and positions", {
  set.seed(91409)
  pe <- make_eeg(n_time = 80, n_channels = 4, sr = 40)
  forward <- eegForwardModel(pe, method = "spherical", n_sources = 5)
  localized <- eegBeamformer(
    pe,
    forward,
    method = "lcmv",
    output_assay = "beam_fixture"
  )
  power <- S4Vectors::metadata(localized)$beam_fixture[, 1L]
  resolved <- PhysioEEG:::.glassbrain_resolve_sources(localized)

  expect_equal(resolved$sources$amplitude, power, tolerance = 0)
  expect_equal(
    resolved$sources[c("x", "y", "z")],
    forward$source_positions,
    tolerance = 0
  )
  expect_identical(resolved$source, "beamformer")
  expect_identical(
    S4Vectors::metadata(localized)$beamformer_info$output_assay,
    "beam_fixture"
  )
})

test_that("source plotting metadata names cannot overwrite numeric outputs", {
  pe <- make_eeg(n_time = 40, n_channels = 4, sr = 40)
  forward <- eegForwardModel(pe, method = "spherical", n_sources = 5)

  expect_error(
    eegSourceEstimate(pe, forward, output_assay = "source_estimate"),
    "reserved"
  )
  expect_error(
    eegBeamformer(pe, forward, output_assay = "beamformer_info"),
    "reserved"
  )
  expect_error(
    eegSourceEstimate(pe, forward, output_assay = character()),
    "one non-empty"
  )
})

test_that("malformed source plotting metadata fail explicitly", {
  pe <- make_eeg(n_time = 40, n_channels = 4, sr = 40)
  metadata <- S4Vectors::metadata(pe)
  metadata$source_estimate <- list(
    source_positions = .ws914_sources(4L)[c("x", "y", "z")],
    output_assay = c("one", "two"),
    orientation_count = 3L
  )
  metadata$source_plot_default <- "source_estimate"
  S4Vectors::metadata(pe) <- metadata

  expect_error(eegPlotGlassBrain(pe), "incomplete or invalid")
  metadata <- S4Vectors::metadata(pe)
  metadata$source_plot_default <- c("source_estimate", "beamformer")
  S4Vectors::metadata(pe) <- metadata
  expect_error(eegPlotGlassBrain(pe), "invalid default")
})

test_that("legacy source metadata never fabricates glass-brain coordinates", {
  pe <- make_eeg(n_time = 100, n_channels = 4, sr = 100)
  metadata <- S4Vectors::metadata(pe)
  metadata$source_estimate <- list(
    method = "mne",
    lambda = 0.05,
    n_sources = 3L,
    n_source_cols = 9L
  )
  S4Vectors::metadata(pe) <- metadata
  expect_error(eegPlotGlassBrain(pe), "Legacy source metadata")
})

test_that("sparse fallback uses the first view without a second threshold", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  pe <- make_eeg(n_time = 100, n_channels = 4, sr = 100)
  sources <- .ws914_sources(10L)
  expect_warning(
    plot <- eegPlotGlassBrain(
      pe,
      sources,
      views = c("coronal", "axial"),
      threshold_pct = 100
    ),
    "sparse fallback"
  )
  data <- attr(plot, "glassbrain_data")
  expected_id <- sources$source_id[which.max(abs(sources$amplitude))]

  expect_s3_class(plot, "ggplot")
  expect_false(inherits(plot, "patchwork"))
  expect_identical(data$settings$fallback_view, "coronal")
  expect_equal(nrow(data$settings$fallback_data), 1L)
  expect_identical(data$retained_sources$source_id, expected_id)
  expect_equal(
    data$settings$fallback_data$amplitude,
    sources$amplitude[sources$source_id == expected_id]
  )

  tied_sources <- sources
  tied_sources$amplitude[9:10] <- c(-20, 20)
  expect_warning(
    tied_plot <- eegPlotGlassBrain(
      pe,
      tied_sources,
      views = "axial",
      threshold_pct = 100
    ),
    "Only 2 source rows"
  )
  tied_data <- attr(tied_plot, "glassbrain_data")
  expect_equal(nrow(tied_data$settings$fallback_data), 2L)
  expect_setequal(tied_data$retained_sources$source_id, c(9L, 10L))
})

test_that("all-zero non-empty sources render as zero activation", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  pe <- make_eeg(n_time = 100, n_channels = 4, sr = 100)
  sources <- .ws914_sources(8L, amplitude = rep(0, 8L))
  plot <- eegPlotGlassBrain(pe, sources)
  data <- attr(plot, "glassbrain_data")

  expect_s3_class(plot, "patchwork")
  expect_true(data$settings$zero_activation)
  expect_equal(data$settings$retained_count, 8L)
  expect_match(plot[[1L]]$labels$subtitle, "Zero activation")
})

test_that("outline artifact is offline, closed, and hash verified", {
  path <- system.file(
    "extdata",
    "glassbrain_outline.rds",
    package = "PhysioEEG"
  )
  dcf_path <- system.file(
    "extdata",
    "glassbrain_outline.dcf",
    package = "PhysioEEG"
  )
  sha_path <- system.file(
    "extdata",
    "glassbrain_outline.sha256",
    package = "PhysioEEG"
  )
  outline <- readRDS(path)
  metadata <- read.dcf(dcf_path)

  expect_identical(
    names(outline),
    c("view", "path_id", "vertex_order", "x", "y", "structure")
  )
  expect_setequal(unique(outline$view), c("sagittal", "axial", "coronal"))
  expect_true(all(is.finite(outline$x)))
  expect_true(all(is.finite(outline$y)))
  expect_lte(max(abs(c(outline$x, outline$y))), 1)
  for (path_data in split(outline, outline$path_id)) {
    expect_equal(
      unname(as.numeric(path_data[1L, c("x", "y")])),
      unname(as.numeric(path_data[nrow(path_data), c("x", "y")])),
      tolerance = 0
    )
    expect_identical(path_data$vertex_order, seq_len(nrow(path_data)))
  }
  expect_identical(
    unname(metadata[1L, "Anatomical_fidelity"]),
    "schematic"
  )
  expect_identical(unname(metadata[1L, "License"]), "MIT")

  skip_if_not_installed("digest")
  actual <- c(
    digest::digest(file = path, algo = "sha256"),
    digest::digest(file = dcf_path, algo = "sha256")
  )
  expect_identical(actual, .glassbrain_expected_hashes(sha_path))
})

test_that("panel scales, limits, source family, and numeric data are shared", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  pe <- make_eeg(n_time = 100, n_channels = 4, sr = 100)
  sources <- .ws914_sources(15L)
  plot <- eegPlotGlassBrain(pe, sources, threshold_pct = 0)
  data_before <- attr(plot, "glassbrain_data")$panels
  invisible(lapply(seq_len(3L), function(index) {
    ggplot2::ggplot_build(plot[[index]])
  }))
  data_after <- attr(plot, "glassbrain_data")$panels

  expect_identical(data_after, data_before)
  color_limits <- lapply(seq_len(3L), function(index) {
    plot[[index]]$scales$get_scales("colour")$limits
  })
  size_limits <- lapply(seq_len(3L), function(index) {
    plot[[index]]$scales$get_scales("size")$limits
  })
  coord_limits <- lapply(seq_len(3L), function(index) {
    c(plot[[index]]$coordinates$limits$x, plot[[index]]$coordinates$limits$y)
  })
  expect_true(all(vapply(color_limits[-1L], identical, color_limits[[1L]], FUN.VALUE = logical(1))))
  expect_true(all(vapply(size_limits[-1L], identical, size_limits[[1L]], FUN.VALUE = logical(1))))
  expect_true(all(vapply(coord_limits[-1L], identical, coord_limits[[1L]], FUN.VALUE = logical(1))))
  expected_ids <- sources$source_id
  expect_true(all(vapply(data_before, function(panel) {
    all(panel$source_id %in% expected_ids)
  }, logical(1))))
})

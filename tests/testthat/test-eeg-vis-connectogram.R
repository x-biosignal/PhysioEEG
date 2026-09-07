make_connectogram_pe <- function(labels, matrix, col_data = list(),
                                 method = "coherence", band = c(8, 13)) {
  n <- length(labels)
  cd <- S4Vectors::DataFrame(
    label = labels,
    type = rep("EEG", n)
  )
  for (name in names(col_data)) {
    cd[[name]] <- col_data[[name]]
  }
  pe <- PhysioExperiment(
    assays = list(raw = matrix(0, nrow = 20, ncol = n)),
    colData = cd,
    samplingRate = 100
  )
  S4Vectors::metadata(pe)$connectivity <- list(
    matrix = matrix,
    method = method,
    band = band
  )
  pe
}


named_connectivity <- function(values, labels, byrow = TRUE) {
  matrix <- matrix(values, nrow = length(labels), byrow = byrow)
  dimnames(matrix) <- list(labels, labels)
  matrix
}


connectogram_data <- function(plot) {
  attr(plot, "connectogram_data", exact = TRUE)
}


connectogram_edge_id_ref <- function(source, target, directed = FALSE) {
  separator <- if (directed) " -> " else " -- "
  paste0(
    nchar(enc2utf8(source), type = "bytes"),
    ":",
    enc2utf8(source),
    separator,
    nchar(enc2utf8(target), type = "bytes"),
    ":",
    enc2utf8(target)
  )
}


test_that("symmetric matrices emit one canonical edge per unordered pair", {
  skip_if_not_installed("ggplot2")
  labels <- c("Fp1", "Fp2", "Cz", "O1")
  mat <- named_connectivity(c(
    1, 0.8, 0, -0.6,
    0.8, 1, 0.3, 0,
    0, 0.3, 1, 0.2,
    -0.6, 0, 0.2, 1
  ), labels)
  pe <- make_connectogram_pe(labels, mat)

  p <- eegPlotConnectogram(pe, threshold = 0.25)
  data <- connectogram_data(p)

  expect_s3_class(p, "ggplot")
  expect_equal(nrow(data$nodes), 4L)
  expect_equal(
    data$edges$edge_id,
    c(
      connectogram_edge_id_ref("Cz", "Fp2"),
      connectogram_edge_id_ref("Fp1", "Fp2"),
      connectogram_edge_id_ref("Fp1", "O1")
    )
  )
  expect_equal(data$edges$value, c(0.3, 0.8, -0.6))
  expect_false(any(data$edges$directed))
  expect_true(data$settings$symmetric)
  expect_equal(data$settings$method, "coherence")
  expect_equal(data$settings$band, c(8, 13))
})


test_that("asymmetric matrices retain ordered source-to-target edges", {
  skip_if_not_installed("ggplot2")
  labels <- c("F3", "C3", "P3")
  mat <- named_connectivity(c(
    1, 0.2, 0,
    0.9, 1, -0.4,
    0.6, 0, 1
  ), labels)
  pe <- make_connectogram_pe(labels, mat, method = "dtf")

  p <- eegPlotConnectogram(pe, threshold = 0.15, bundle = FALSE)
  data <- connectogram_data(p)

  expect_false(data$settings$symmetric)
  expect_equal(data$settings$matrix_direction,
               "rows_target_columns_source")
  expect_setequal(
    data$edges$edge_id,
    c(
      connectogram_edge_id_ref("F3", "C3", TRUE),
      connectogram_edge_id_ref("F3", "P3", TRUE),
      connectogram_edge_id_ref("C3", "F3", TRUE),
      connectogram_edge_id_ref("P3", "C3", TRUE)
    )
  )
  expect_equal(
    data$edges$value[match(
      connectogram_edge_id_ref("F3", "C3", TRUE),
      data$edges$edge_id
    )],
    mat["C3", "F3"]
  )
  expect_equal(
    data$edges$value[match(
      connectogram_edge_id_ref("P3", "C3", TRUE),
      data$edges$edge_id
    )],
    mat["C3", "P3"]
  )
  expect_true(all(data$edges$directed))
  expect_s3_class(p$layers[[1L]]$geom_params$arrow, "arrow")
})


test_that("thresholding is strict, absolute, and excludes self-connections", {
  skip_if_not_installed("ggplot2")
  labels <- c("F3", "C3", "P3")
  mat <- named_connectivity(c(
    99, 0.5, -0.6,
    0.5, -99, 0,
    -0.6, 0, 42
  ), labels)
  pe <- make_connectogram_pe(labels, mat)

  zero <- connectogram_data(eegPlotConnectogram(pe, threshold = NULL))
  strict <- connectogram_data(eegPlotConnectogram(pe, threshold = 0.5))

  expect_equal(nrow(zero$edges), 2L)
  expect_equal(
    strict$edges$edge_id,
    connectogram_edge_id_ref("F3", "P3")
  )
  expect_equal(strict$edges$value, -0.6)
  expect_false(any(strict$edges$source_id == strict$edges$target_id))
  expect_equal(strict$settings$threshold, 0.5)
})


test_that("matrix and channel identities fail loudly when ambiguous", {
  skip_if_not_installed("ggplot2")
  labels <- c("F3", "C3")
  mat <- named_connectivity(c(1, 0.2, 0.2, 1), labels)
  pe <- make_connectogram_pe(labels, mat)

  bad <- mat
  colnames(bad) <- NULL
  S4Vectors::metadata(pe)$connectivity$matrix <- bad
  expect_error(eegPlotConnectogram(pe), "both be present")

  bad <- mat
  colnames(bad) <- rev(labels)
  S4Vectors::metadata(pe)$connectivity$matrix <- bad
  expect_error(eegPlotConnectogram(pe), "identical in the same order")

  bad <- mat
  dimnames(bad) <- list(c("F3", "F3"), c("F3", "F3"))
  S4Vectors::metadata(pe)$connectivity$matrix <- bad
  expect_error(eegPlotConnectogram(pe), "complete, unique")

  bad <- matrix(c(1, 0.2, 0.2, 1), 2)
  S4Vectors::metadata(pe)$connectivity$matrix <- bad
  expect_no_error(eegPlotConnectogram(pe))
  SummarizedExperiment::colData(pe)$label <- NULL
  expect_error(eegPlotConnectogram(pe), "requires exact unique")

  pe2 <- make_connectogram_pe(labels, mat)
  S4Vectors::metadata(pe2)$connectivity$matrix <- matrix(1:6, 2, 3)
  expect_error(eegPlotConnectogram(pe2), "square")
  S4Vectors::metadata(pe2)$connectivity$matrix <- mat
  S4Vectors::metadata(pe2)$connectivity$matrix[1, 2] <- Inf
  expect_error(eegPlotConnectogram(pe2), "finite")

  one <- named_connectivity(1, "F3")
  pe1 <- make_connectogram_pe("F3", one)
  expect_error(eegPlotConnectogram(pe1), "at least 2 x 2")
})


test_that("matrix names and colData labels reconcile by exact identity", {
  skip_if_not_installed("ggplot2")
  labels <- c("F3", "C3", "P3")
  mat <- named_connectivity(c(
    1, 0.2, 0.3,
    0.2, 1, 0.4,
    0.3, 0.4, 1
  ), labels)
  pe <- make_connectogram_pe(rev(labels), mat)
  expect_no_error(eegPlotConnectogram(pe))

  SummarizedExperiment::colData(pe)$label <- c("F3", "C3", "O1")
  expect_error(eegPlotConnectogram(pe), "same channels")

  pe_dup <- make_connectogram_pe(labels, mat)
  SummarizedExperiment::colData(pe_dup)$label <- c("F3", "F3", "P3")
  expect_error(eegPlotConnectogram(pe_dup), "complete unique")
})


test_that("recognized EEG labels map to exact hemisphere and lobe classes", {
  skip_if_not_installed("ggplot2")
  labels <- c(
    "Fp1", "AF2", "Fz", "FC3", "CP4", "Pz", "PO7",
    "O2", "I1", "FT8", "T3", "TP4", "EEG7"
  )
  mat <- diag(length(labels))
  dimnames(mat) <- list(labels, labels)
  pe <- make_connectogram_pe(labels, mat)
  nodes <- connectogram_data(eegPlotConnectogram(pe))$nodes
  nodes <- nodes[match(labels, nodes$label), ]

  expect_equal(
    nodes$hemisphere,
    c("left", "right", "midline", "left", "right", "midline",
      "left", "right", "left", "right", "left", "right", "unknown")
  )
  expect_equal(
    nodes$lobe,
    c("frontal", "frontal", "frontal", "central", "central",
      "parietal", "parietal", "occipital", "occipital", "temporal",
      "temporal", "temporal", "unknown")
  )
})


test_that("explicit anatomy overrides inference and rejects unknown values", {
  skip_if_not_installed("ggplot2")
  labels <- c("F3", "C4", "X")
  mat <- diag(3)
  dimnames(mat) <- list(labels, labels)
  pe <- make_connectogram_pe(
    labels,
    mat,
    col_data = list(
      hemi = c("RIGHT", NA, "UNKNOWN"),
      lobe = c("TEMPORAL", "CENTRAL", NA)
    )
  )
  nodes <- connectogram_data(eegPlotConnectogram(pe))$nodes
  nodes <- nodes[match(labels, nodes$label), ]
  expect_equal(nodes$hemisphere, c("right", "right", "unknown"))
  expect_equal(nodes$lobe, c("temporal", "central", "unknown"))

  SummarizedExperiment::colData(pe)$hemi[1] <- "ipsilateral"
  expect_error(eegPlotConnectogram(pe), "Unrecognized explicit hemisphere")
})


test_that("anatomical orders are invariant to common input permutation", {
  skip_if_not_installed("ggplot2")
  labels <- c("O2", "Fp1", "Cz", "P3", "X", "F4")
  mat <- outer(seq_along(labels), seq_along(labels), function(i, j) {
    ifelse(i == j, 1, (i + j) / 20)
  })
  dimnames(mat) <- list(labels, labels)
  permutation <- c(4, 1, 6, 3, 2, 5)

  for (ordering in c("hemisphere", "lobe")) {
    first <- connectogram_data(
      eegPlotConnectogram(make_connectogram_pe(labels, mat), order = ordering)
    )$nodes
    second <- connectogram_data(
      eegPlotConnectogram(
        make_connectogram_pe(
          labels[permutation],
          mat[permutation, permutation]
        ),
        order = ordering
      )
    )$nodes
    expect_equal(first$label, second$label)
    expect_equal(first[, c("plot_index", "angle", "x", "y")],
                 second[, c("plot_index", "angle", "x", "y")])
  }
})


test_that("cluster order is permutation invariant and tied matrices diagnose", {
  skip_if_not_installed("ggplot2")
  labels <- c("D", "A", "C", "B")
  mat <- named_connectivity(c(
    1, 0.9, 0.1, 0.2,
    0.9, 1, 0.2, 0.1,
    0.1, 0.2, 1, 0.8,
    0.2, 0.1, 0.8, 1
  ), labels)
  permutation <- c(3, 1, 4, 2)
  first <- connectogram_data(
    eegPlotConnectogram(make_connectogram_pe(labels, mat), order = "cluster")
  )
  second <- connectogram_data(
    eegPlotConnectogram(
      make_connectogram_pe(
        labels[permutation],
        mat[permutation, permutation]
      ),
      order = "cluster"
    )
  )
  expect_equal(first$nodes$label, second$nodes$label)

  tied <- matrix(0.5, 4, 4)
  diag(tied) <- 1
  dimnames(tied) <- list(labels, labels)
  tied_data <- connectogram_data(
    eegPlotConnectogram(make_connectogram_pe(labels, tied), order = "cluster")
  )
  expect_equal(tied_data$nodes$label, sort(labels))
  expect_true("cluster_all_tied" %in%
                tied_data$settings$diagnostics$code)
})


test_that("module declarations resolve exactly by channel label", {
  skip_if_not_installed("ggplot2")
  labels <- c("F3", "C3", "P3", "O1")
  mat <- diag(4)
  dimnames(mat) <- list(labels, labels)
  pe <- make_connectogram_pe(labels, mat)
  modules <- c(P3 = "posterior", F3 = "anterior",
               O1 = "posterior", C3 = "anterior")
  nodes <- connectogram_data(
    eegPlotConnectogram(pe, modules = modules)
  )$nodes
  expect_equal(
    stats::setNames(nodes$module, nodes$label)[labels],
    c(F3 = "anterior", C3 = "anterior",
      P3 = "posterior", O1 = "posterior")
  )

  module_df <- data.frame(
    label = rev(labels),
    module = c("posterior", "posterior", "anterior", "anterior")
  )
  df_nodes <- connectogram_data(
    eegPlotConnectogram(pe, modules = module_df)
  )$nodes
  expect_equal(
    stats::setNames(df_nodes$module, df_nodes$label)[labels],
    c(F3 = "anterior", C3 = "anterior",
      P3 = "posterior", O1 = "posterior")
  )

  expect_error(eegPlotConnectogram(pe, modules = unname(modules)),
               "uniquely named")
  expect_error(eegPlotConnectogram(pe, modules = modules[-1]),
               "uniquely named")
  expect_error(
    eegPlotConnectogram(pe, modules = c(modules, Extra = "other")),
    "uniquely named"
  )
  duplicate <- modules
  names(duplicate)[2] <- names(duplicate)[1]
  expect_error(eegPlotConnectogram(pe, modules = duplicate),
               "uniquely named")
  empty <- modules
  empty[1] <- ""
  expect_error(eegPlotConnectogram(pe, modules = empty), "non-empty")
})


test_that("non-contiguous and wrap-around module runs remain separate", {
  skip_if_not_installed("ggplot2")
  labels <- LETTERS[1:6]
  mat <- diag(6)
  dimnames(mat) <- list(labels, labels)
  pe <- make_connectogram_pe(labels, mat)

  alternating <- stats::setNames(rep(c("u", "v"), 3), labels)
  arcs <- connectogram_data(
    eegPlotConnectogram(pe, modules = alternating)
  )$module_arcs
  run_rows <- arcs[!duplicated(arcs$module_run_id), ]
  expect_equal(
    stats::setNames(as.integer(table(run_rows$module)),
                    names(table(run_rows$module))),
    c(u = 3L, v = 3L)
  )

  wrapping <- stats::setNames(c("u", "v", "v", "v", "v", "u"), labels)
  wrap_arcs <- connectogram_data(
    eegPlotConnectogram(pe, modules = wrapping)
  )$module_arcs
  wrap_runs <- wrap_arcs[!duplicated(wrap_arcs$module_run_id), ]
  expect_equal(nrow(wrap_runs), 2L)
  u <- wrap_runs[wrap_runs$module == "u", ]
  expect_equal(u$start_plot_index, 6L)
  expect_equal(u$end_plot_index, 1L)
})


test_that("straight and bundled paths preserve edges and exact endpoints", {
  skip_if_not_installed("ggplot2")
  labels <- c("F3", "C3", "P3", "O1")
  mat <- named_connectivity(c(
    1, 0.8, 0, -0.4,
    0.8, 1, 0.5, 0,
    0, 0.5, 1, 0.7,
    -0.4, 0, 0.7, 1
  ), labels)
  pe <- make_connectogram_pe(labels, mat)
  modules <- c(F3 = "front", C3 = "front", P3 = "back", O1 = "back")
  straight <- connectogram_data(
    eegPlotConnectogram(pe, modules = modules, bundle = FALSE)
  )
  bundled <- connectogram_data(
    eegPlotConnectogram(pe, modules = modules, bundle = TRUE)
  )
  repeated <- connectogram_data(
    eegPlotConnectogram(pe, modules = modules, bundle = TRUE)
  )

  expect_equal(straight$edges, bundled$edges)
  expect_equal(unique(straight$paths$edge_id),
               unique(bundled$paths$edge_id))
  expect_equal(bundled$paths, repeated$paths)
  expect_true(all(table(straight$paths$edge_id) == 2L))
  expect_true(all(table(bundled$paths$edge_id) == 40L))

  for (edge_id in straight$edges$edge_id) {
    edge <- straight$edges[straight$edges$edge_id == edge_id, ]
    path <- straight$paths[straight$paths$edge_id == edge_id, ]
    source <- straight$nodes[straight$nodes$node_id == edge$source_id, ]
    target <- straight$nodes[straight$nodes$node_id == edge$target_id, ]
    expect_equal(c(path$x[1], path$y[1]), c(source$x, source$y))
    expect_equal(c(tail(path$x, 1), tail(path$y, 1)),
                 c(target$x, target$y))
  }
})


test_that("bundled paths equal the declared cubic Bezier construction", {
  skip_if_not_installed("ggplot2")
  labels <- c("A", "B", "C", "D")
  mat <- named_connectivity(c(
    1, 0.8, 0, 0,
    0.8, 1, 0, 0,
    0, 0, 1, 0,
    0, 0, 0, 1
  ), labels)
  modules <- c(A = "left", B = "right", C = "left", D = "right")
  data <- connectogram_data(
    eegPlotConnectogram(
      make_connectogram_pe(labels, mat),
      modules = modules,
      bundle = TRUE
    )
  )
  edge <- data$edges[1, ]
  path <- data$paths[data$paths$edge_id == edge$edge_id, ]
  source <- data$nodes[data$nodes$node_id == edge$source_id, ]
  target <- data$nodes[data$nodes$node_id == edge$target_id, ]
  source_hub <- data$settings$module_hubs[
    data$settings$module_hubs$module == source$module, ]
  target_hub <- data$settings$module_hubs[
    data$settings$module_hubs$module == target$module, ]
  t <- seq(0, 1, length.out = 40)
  omt <- 1 - t
  ref_x <- omt^3 * source$x + 3 * omt^2 * t * source_hub$x +
    3 * omt * t^2 * target_hub$x + t^3 * target$x
  ref_y <- omt^3 * source$y + 3 * omt^2 * t * source_hub$y +
    3 * omt * t^2 * target_hub$y + t^3 * target$y

  expect_equal(path$x, ref_x, tolerance = 1e-15)
  expect_equal(path$y, ref_y, tolerance = 1e-15)
  expect_true("antipodal_module_hub" %in%
                data$settings$diagnostics$code)
})


test_that("resolved tables are the immutable data used by plot layers", {
  skip_if_not_installed("ggplot2")
  labels <- c("F3", "C3", "P3")
  mat <- named_connectivity(c(
    1, -0.8, 0,
    -0.8, 1, 0.4,
    0, 0.4, 1
  ), labels)
  pe <- make_connectogram_pe(labels, mat)
  p <- eegPlotConnectogram(pe, bundle = TRUE)
  data <- connectogram_data(p)

  expect_equal(p$layers[[1L]]$data, data$paths)
  expect_equal(p$layers[[2L]]$data, data$module_arcs)
  expect_equal(p$layers[[3L]]$data, data$nodes)
  expect_equal(p$layers[[4L]]$data, data$nodes)
  expect_equal(
    names(data$settings$module_colors),
    sort(unique(data$nodes$module))
  )
  expect_no_warning(built <- ggplot2::ggplot_build(p))
  expect_equal(length(built$data), 4L)
  expect_equal(p$layers[[1L]]$data$value,
               rep(data$edges$value, each = 40L))
  expect_equal(p$layers[[1L]]$data$abs_value,
               abs(p$layers[[1L]]$data$value))
})


test_that("edge identifiers remain unique for delimiter-bearing labels", {
  skip_if_not_installed("ggplot2")
  labels <- c("A", "B -- C", "A -- B", "C")
  mat <- diag(4)
  dimnames(mat) <- list(labels, labels)
  mat["A", "B -- C"] <- mat["B -- C", "A"] <- 0.8
  mat["A -- B", "C"] <- mat["C", "A -- B"] <- -0.7
  pe <- make_connectogram_pe(labels, mat)
  data <- connectogram_data(eegPlotConnectogram(pe, bundle = TRUE))

  expect_equal(nrow(data$edges), 2L)
  expect_identical(anyDuplicated(data$edges$edge_id), 0L)
  expect_equal(length(unique(data$paths$edge_id)), 2L)
  expect_setequal(
    data$edges$edge_id,
    c(
      connectogram_edge_id_ref("A", "B -- C"),
      connectogram_edge_id_ref("A -- B", "C")
    )
  )
})


test_that("argument contracts reject invalid threshold, bundle, and modules", {
  skip_if_not_installed("ggplot2")
  labels <- c("F3", "C3")
  mat <- named_connectivity(c(1, 0.2, 0.2, 1), labels)
  pe <- make_connectogram_pe(labels, mat)

  expect_error(eegPlotConnectogram(pe, threshold = -1), "non-negative")
  expect_error(eegPlotConnectogram(pe, threshold = Inf), "finite")
  expect_error(eegPlotConnectogram(pe, threshold = c(0, 1)), "scalar")
  expect_error(eegPlotConnectogram(pe, bundle = NA), "non-missing")
  expect_error(eegPlotConnectogram(pe, bundle = 1), "logical")
  expect_error(eegPlotConnectogram(pe, order = "input"))
  expect_error(eegPlotConnectogram(pe, order = "h"), "exactly one")
  expect_error(eegPlotConnectogram(pe, modules = list(F3 = "a", C3 = "b")),
               "modules must")
})

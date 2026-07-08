test_that("placeLabels produces a conflict-free layout on the example data", {
  data("exampleMascarade", package = "mascarade")
  g <- .buildTestGeom(exampleMascarade$dims, exampleMascarade$clusters)
  P <- placeLabels(g$geom, g$xlim, g$ylim, g$hw, g$hh, g$char_h)

  expect_equal(nrow(P), nrow(g$geom$poi))          # every cluster placed (feasibility)
  s <- layoutScore(P)
  expect_equal(unname(s["bb"]), 0)                 # no box-box overlap
  expect_equal(unname(s["ll"]), 0)                 # no leader-leader crossing
  expect_equal(unname(s["lb"]), 0)                 # no leader through a foreign box

  # regression gate: no worse than the recorded golden (conflict-free + no longer total)
  golden <- readRDS(test_path("fixtures", "golden_scores.rds"))$example_default
  expect_false(scoreBetter(golden, s))             # golden must not beat the current layout
})

test_that("placement stays conflict-free across box scales (view-change property)", {
  data("exampleMascarade", package = "mascarade")
  for (frac in c(0.03, 0.06, 0.09)) {              # smaller/larger boxes = zoom in/out
    g <- .buildTestGeom(exampleMascarade$dims, exampleMascarade$clusters, char_frac = frac)
    P <- placeLabels(g$geom, g$xlim, g$ylim, g$hw, g$hh, g$char_h)
    s <- layoutScore(P)
    expect_true(s["bb"] == 0 && s["ll"] == 0 && s["lb"] == 0,
                info = sprintf("char_frac=%.2f gave bb/ll/lb = %d/%d/%d", frac, s["bb"], s["ll"], s["lb"]))
  }
})

test_that("a single-label degenerate case does not error", {
  data("exampleMascarade", package = "mascarade")
  keep <- exampleMascarade$clusters == exampleMascarade$clusters[1]
  g <- .buildTestGeom(exampleMascarade$dims[keep, , drop = FALSE],
                      exampleMascarade$clusters[keep])
  P <- placeLabels(g$geom, g$xlim, g$ylim, g$hw, g$hh, g$char_h)
  expect_equal(nrow(P), 1L)
})

test_that("fancyMask renders a plot end-to-end (draw-stage placement)", {
  skip_if_not_installed("ggplot2")
  data("exampleMascarade", package = "mascarade")
  mt <- generateMask(dims = exampleMascarade$dims, clusters = exampleMascarade$clusters)
  p <- ggplot2::ggplot(data.frame(x = exampleMascarade$dims[, 1],
                                  y = exampleMascarade$dims[, 2],
                                  cluster = exampleMascarade$clusters)) +
    ggplot2::geom_point(ggplot2::aes(x, y, color = cluster)) +
    fancyMask(mt, ratio = 1, cols = "inherit")
  tmp <- tempfile(fileext = ".png")
  expect_error(ggplot2::ggsave(tmp, p, width = 8, height = 6, dpi = 90), NA)
  unlink(tmp)
})

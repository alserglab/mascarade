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

test_that("both leader styles (cl, cm) produce conflict-free layouts", {
  data("exampleMascarade", package = "mascarade")
  g <- .buildTestGeom(exampleMascarade$dims, exampleMascarade$clusters)
  for (ct in c("cl", "cm")) {
    P <- placeLabels(g$geom, g$xlim, g$ylim, g$hw, g$hh, g$char_h, con_type = ct)
    s <- layoutScore(P)
    expect_true(s["bb"] == 0 && s["ll"] == 0 && s["lb"] == 0,
                info = sprintf("con_type=%s gave bb/ll/lb = %d/%d/%d", ct, s["bb"], s["ll"], s["lb"]))
  }
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

test_that("twoMoveBnB resolves an input conflict that needs a longer candidate", {
  # Global candidate list (0-based indices used by the kernel):
  #   0 = A: label 1, SHORT leader, box overlaps B  -> the conflicting start pick
  #   1 = C: label 1, LONG  leader, box far away    -> the resolution (length-pruned!)
  #   2 = B: label 2, its only candidate
  # Leaders point away from all boxes so the only conflict is A-B box-box overlap.
  cxmin <- c(0, 10, 1); cxmax <- c(2, 12, 3)
  cymin <- c(0, 10, 1); cymax <- c(2, 12, 3)
  ex <- c(1, 11, 2);  ey <- c(-1, 9, 4)
  tx <- c(1, 11, 2);  ty <- c(-2, 8, 5)
  len <- c(1, 5, 1)                         # C (idx1) is far longer than A (idx0)
  rows <- list(c(0L, 1L), c(2L))            # label1 -> {A, C} (length-ascending), label2 -> {B}
  init <- c(0L, 2L)                         # start at A, B  -> A and B overlap (bb conflict)

  # sanity: the starting layout really is in conflict
  P0 <- data.frame(cxmin = cxmin[init + 1L], cxmax = cxmax[init + 1L],
                   cymin = cymin[init + 1L], cymax = cymax[init + 1L],
                   ex = ex[init + 1L], ey = ey[init + 1L], tx = tx[init + 1L], ty = ty[init + 1L],
                   len = len[init + 1L])
  expect_gt(unname(layoutScore(P0)["bb"]), 0)

  sel <- mascarade:::twoMoveBnB(cxmin, cxmax, cymin, cymax, ex, ey, tx, ty, len,
                                rows, init, maxpass = 50L, sq = TRUE) + 1L
  P <- data.frame(cxmin = cxmin[sel], cxmax = cxmax[sel], cymin = cymin[sel], cymax = cymax[sel],
                  ex = ex[sel], ey = ey[sel], tx = tx[sel], ty = ty[sel], len = len[sel])
  s <- layoutScore(P)
  expect_equal(unname(s["bb"]), 0)          # conflict resolved (label 1 moved A -> C)
  expect_equal(unname(s["ll"]), 0)
  expect_equal(unname(s["lb"]), 0)
})

test_that("hungarian solves a small assignment to the known optimum", {
  # each row's unique minimum sits in a distinct column -> optimal is that permutation, cost 3
  cost <- matrix(c(9, 1, 9,
                   9, 9, 1,
                   1, 9, 9), nrow = 3, byrow = TRUE)
  res <- mascarade:::hungarian(cost)          # res[i] = 0-indexed column assigned to row i
  expect_equal(res, c(1L, 2L, 0L))
  total <- sum(vapply(seq_len(3), function(i) cost[i, res[i] + 1L], 0))
  expect_equal(total, 3)
})

test_that("placeLabels produces a conflict-free layout on the example data", {
  scene <- loadPlacementScene("example")
  placement <- placeLabels(scene$geom, scene$xlim, scene$ylim,
                           scene$halfWidth, scene$halfHeight, scene$charHeight)

  expect_equal(nrow(placement), nrow(scene$geom$poi))  # every cluster placed

  score <- layoutScore(placement)
  expect_equal(unname(score["boxBox"]), 0)
  expect_equal(unname(score["leaderLeader"]), 0)
  expect_equal(unname(score["leaderBox"]), 0)
  expectNoWorseThan(score, loadGoldenScore("example_default"), tolerance = 0.05)
})

test_that("all leader styles (ledge, line, box) produce conflict-free layouts", {
  scene <- loadPlacementScene("example")
  for (connectorType in c("ledge", "line", "box")) {
    placement <- placeLabels(scene$geom, scene$xlim, scene$ylim,
                             scene$halfWidth, scene$halfHeight, scene$charHeight,
                             con_type = connectorType)
    score <- layoutScore(placement)
    expect_true(
      score["boxBox"] == 0 && score["leaderLeader"] == 0 && score["leaderBox"] == 0,
      info = sprintf("con_type=%s gave boxBox/leaderLeader/leaderBox = %d/%d/%d",
                     connectorType, score["boxBox"], score["leaderLeader"],
                     score["leaderBox"])
    )
  }
})

test_that("a single label is placed in free space through the normal pipeline", {
  scene <- loadPlacementScene("singleCluster")
  placement <- placeLabels(scene$geom, scene$xlim, scene$ylim,
                           scene$halfWidth, scene$halfHeight, scene$charHeight)

  expect_equal(nrow(placement), 1L)
  expect_true(all(is.finite(c(placement$cx, placement$cy,
                              placement$ex, placement$ey,
                              placement$bx, placement$by))))

  # The unpadded label box lands clear of its own cluster.
  box <- list(list(
    x = placement$cx + c(-1, 1, 1, -1) * placement$hw,
    y = placement$cy + c(-1, -1, 1, 1) * placement$hh
  ))
  cluster <- list(list(x = scene$geom$polysx[[1]], y = scene$geom$polysy[[1]]))
  expect_length(polyclip::polyclip(box, cluster, "intersection"), 0)
})

test_that("an empty view (no clusters) returns an empty layout", {
  geom <- list(poi = matrix(numeric(0), 0, 2),
               rtree = buildBoxFit(list(), list()),
               polysx = list(), polysy = list())
  placement <- placeLabels(geom, c(-1, 1), c(-1, 1),
                           hw = numeric(0), hh = numeric(0), char_h = 1)
  expect_equal(nrow(placement), 0L)
  expect_true(all(c("cx", "cy", "ex", "ey", "corner", "bx", "by",
                    "cxmin", "cxmax", "cymin", "cymax") %in% names(placement)))
})

test_that("fancyMask renders a plot end-to-end (draw-stage placement)", {
  skip_if_not_installed("ggplot2")
  data("exampleMascarade", package = "mascarade")
  plotData <- data.frame(x = exampleMascarade$dims[, 1],
                         y = exampleMascarade$dims[, 2],
                         cluster = exampleMascarade$clusters)
  plot <- ggplot2::ggplot(plotData) +
    ggplot2::geom_point(ggplot2::aes(x, y, color = cluster)) +
    fancyMask(exampleMask(), ratio = 1, cols = "inherit")

  outputFile <- tempfile(fileext = ".png")
  on.exit(unlink(outputFile), add = TRUE)
  expect_no_error(ggplot2::ggsave(outputFile, plot, width = 8, height = 6, dpi = 90))
})

test_that("twoMoveSweep kernel resolves an input conflict that needs a longer candidate", {
  # Global candidate list (0-based indices used by the kernel):
  #   0 = A: label 1, SHORT leader, box overlaps B  -> the conflicting start pick
  #   1 = C: label 1, LONG  leader, box far away    -> the resolution (length-pruned!)
  #   2 = B: label 2, its only candidate
  # Leaders point away from all boxes so the only conflict is A-B box-box overlap.
  cxmin <- c(0, 10, 1)
  cxmax <- c(2, 12, 3)
  cymin <- c(0, 10, 1)
  cymax <- c(2, 12, 3)
  ex <- c(1, 11, 2)
  ey <- c(-1, 9, 4)
  tx <- c(1, 11, 2)
  ty <- c(-2, 8, 5)
  len <- c(1, 5, 1)               # C (index 1) is far longer than A (index 0)
  candidatesPerLabel <- list(c(0L, 1L), c(2L))  # label 1 -> {A, C}, label 2 -> {B}
  initial <- c(0L, 2L)            # start at A, B -> A and B overlap

  asLayout <- function(selection) {
    data.frame(cxmin = cxmin[selection], cxmax = cxmax[selection],
               cymin = cymin[selection], cymax = cymax[selection],
               ex = ex[selection], ey = ey[selection],
               tx = tx[selection], ty = ty[selection],
               len = len[selection])
  }

  expect_gt(unname(layoutScore(asLayout(initial + 1L))["boxBox"]), 0)  # starts in conflict

  selection <- mascarade:::twoMoveSweepKernel(cxmin, cxmax, cymin, cymax,
                                              ex, ey, tx, ty, len,
                                              candidatesPerLabel, initial,
                                              maxpass = 50L, sq = TRUE) + 1L
  score <- layoutScore(asLayout(selection))
  expect_equal(unname(score["boxBox"]), 0)   # resolved: label 1 moved A -> C
  expect_equal(unname(score["leaderLeader"]), 0)
  expect_equal(unname(score["leaderBox"]), 0)
})

test_that("hungarian solves a small assignment to the known optimum", {
  # Each row's unique minimum sits in a distinct column -> optimum is that permutation, cost 3.
  cost <- matrix(c(9, 1, 9,
                   9, 9, 1,
                   1, 9, 9), nrow = 3, byrow = TRUE)
  assignment <- mascarade:::hungarian(cost)  # assignment[i] = 0-indexed column for row i
  expect_equal(assignment, c(1L, 2L, 0L))

  total <- sum(vapply(seq_len(3), function(i) cost[i, assignment[i] + 1L], 0))
  expect_equal(total, 3)
})

test_that(".sideColumn on a crowded column lays labels on a uniform tallest-box grid", {
  # A column with more labels than the viewport can hold (m > capacity) takes the
  # crowded path: exactly m slots, one tallest-box tall, centred on the pole span
  # and extended past the viewport. boxH = slotH = 1, viewH = 3 -> capacity = 3 < m = 6.
  set.seed(1)
  nLabels <- 6
  scene <- list(poi = cbind(runif(nLabels, -3, -1), runif(nLabels, 0, 3)),
                hh = rep(0.5, nLabels), ylim = c(0, 3), hardPad = 0)
  slotHeight <- max(2 * scene$hh)

  placement <- mascarade:::.sideColumn(scene, seq_len(nLabels), Xline = -3, side = -1)

  expect_equal(nrow(placement), nLabels)     # every label placed, one per slot
  centres <- sort(placement$cy)
  expect_equal(diff(centres), rep(slotHeight, nLabels - 1))       # uniform pitch = tallest box
  expect_gt(max(centres) - min(centres), diff(scene$ylim))        # grid extends past viewport
  expect_equal(mean(range(centres)), mean(scene$ylim), tolerance = 1e-9)  # centred on it

  # hardPad widens the pitch to the tallest padded box (box + 2*hardPad).
  scene$hardPad <- 0.25
  padded <- sort(
    mascarade:::.sideColumn(scene, seq_len(nLabels), Xline = -3, side = -1)$cy
  )
  expect_equal(diff(padded), rep(slotHeight + 2 * 0.25, nLabels - 1))
})

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

test_that("all leader styles (ledge, line, box) produce conflict-free layouts", {
  data("exampleMascarade", package = "mascarade")
  g <- .buildTestGeom(exampleMascarade$dims, exampleMascarade$clusters)
  for (ct in c("ledge", "line", "box")) {
    P <- placeLabels(g$geom, g$xlim, g$ylim, g$hw, g$hh, g$char_h, con_type = ct)
    s <- layoutScore(P)
    expect_true(s["bb"] == 0 && s["ll"] == 0 && s["lb"] == 0,
                info = sprintf("con_type=%s gave bb/ll/lb = %d/%d/%d", ct, s["bb"], s["ll"], s["lb"]))
  }
})

test_that("a single label is placed in free space through the normal pipeline", {
  data("exampleMascarade", package = "mascarade")
  keep <- exampleMascarade$clusters == exampleMascarade$clusters[1]
  g <- .buildTestGeom(exampleMascarade$dims[keep, , drop = FALSE],
                      exampleMascarade$clusters[keep])
  P <- placeLabels(g$geom, g$xlim, g$ylim, g$hw, g$hh, g$char_h)

  expect_equal(nrow(P), 1L)
  expect_true(all(is.finite(c(P$cx, P$cy, P$ex, P$ey, P$bx, P$by))))

  # feasibility: the (unpadded) label box lands clear of its own cluster, like any other label
  box <- list(list(x = c(P$cx - P$hw, P$cx + P$hw, P$cx + P$hw, P$cx - P$hw),
                   y = c(P$cy - P$hh, P$cy - P$hh, P$cy + P$hh, P$cy + P$hh)))
  cluster <- list(list(x = g$geom$polysx[[1]], y = g$geom$polysy[[1]]))
  expect_length(polyclip::polyclip(box, cluster, "intersection"), 0)
})

test_that("an empty view (no clusters) returns an empty layout", {
  geom <- list(poi = matrix(numeric(0), 0, 2),
               rtree = buildBoxFit(list(), list()),
               polysx = list(), polysy = list())
  P <- placeLabels(geom, c(-1, 1), c(-1, 1),
                   hw = numeric(0), hh = numeric(0), char_h = 1)
  expect_equal(nrow(P), 0L)
  expect_true(all(c("cx", "cy", "ex", "ey", "corner", "bx", "by",
                    "cxmin", "cxmax", "cymin", "cymax") %in% names(P)))
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
  expect_no_error(ggplot2::ggsave(tmp, p, width = 8, height = 6, dpi = 90))
  unlink(tmp)
})

test_that("twoMoveSweep kernel resolves an input conflict that needs a longer candidate", {
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

  sel <- mascarade:::twoMoveSweepKernel(cxmin, cxmax, cymin, cymax, ex, ey, tx, ty, len,
                                        rows, init, maxpass = 50L, sq = TRUE) + 1L
  P <- data.frame(cxmin = cxmin[sel], cxmax = cxmax[sel], cymin = cymin[sel], cymax = cymax[sel],
                  ex = ex[sel], ey = ey[sel], tx = tx[sel], ty = ty[sel], len = len[sel])
  s <- layoutScore(P)
  expect_equal(unname(s["bb"]), 0)          # conflict resolved (label 1 moved A -> C)
  expect_equal(unname(s["ll"]), 0)
  expect_equal(unname(s["lb"]), 0)
})

test_that("packLen stacks in order without overlap and minimises true leader length", {
  set.seed(1)
  K <- 6
  px <- runif(K, 1, 12); py <- runif(K, -8, 8); h <- rep(1, K); gap <- 0.25
  ord <- order(-py)                                   # fixed top-to-bottom order
  dx <- 0 - px[ord]; slot <- gap / 5
  cy <- mascarade:::packLen(dx, py[ord], h[ord], gap, slot, -10, 10)

  sep <- (h[ord][-K] + h[ord][-1]) / 2 + gap
  expect_true(all(diff(cy) < 0))                      # centres run strictly top-to-bottom
  expect_true(all(-diff(cy) >= sep - 1e-9))           # separation respected (no box overlap)

  # local optimality (necessary + sufficient here: separable convex cost on an ordered chain):
  # perturbing any centre by one slot, keeping feasibility, never shortens the total.
  totalLen <- function(y) sum(sqrt(dx^2 + (y - py[ord])^2))
  base <- totalLen(cy)
  for (i in seq_len(K)) {
    for (d in c(-slot, slot)) {
      y2 <- cy; y2[i] <- y2[i] + d
      if (all(-diff(y2) >= sep - 1e-9)) {
        expect_gte(totalLen(y2), base - 1e-9)
      }
    }
  }
})

test_that("packLen leaves well-separated labels on their poles but spreads crowded ones", {
  gap <- 0.25; slot <- gap / 5; h <- c(1, 1); dx <- c(2, 2)   # equal horizontal offset
  # two poles far enough apart (> one box + gap): no packing needed, centres land on the poles
  far <- mascarade:::packLen(dx, c(3, -3), h, gap, slot, -10, 10)
  expect_equal(far, c(3, -3), tolerance = slot)
  expect_true(all(-diff(far) >= (h[1] + h[2]) / 2 + gap - 1e-9))

  # two poles closer than the minimum separation: pushed symmetrically apart about their midpoint
  # (a real horizontal offset makes the length cost strictly convex, so the optimum is unique)
  near <- mascarade:::packLen(dx, c(0.4, -0.4), h, gap, slot, -10, 10)
  expect_equal(-diff(near), (h[1] + h[2]) / 2 + gap, tolerance = slot)  # exactly touching + gap
  expect_equal(mean(near), 0, tolerance = slot)                        # centred on pole midpoint
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

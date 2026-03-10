library(polyclip)

# Helper: circle approximation with n vertices (CCW)
circle_poly <- function(n = 100, r = 1) {
  theta <- seq(0, 2 * pi, length.out = n + 1)[-(n + 1)]
  list(x = r * cos(theta), y = r * sin(theta))
}

test_that("simplify_outer reduces vertex count on a non-convex polygon", {
  # Build a circle with a small inward notch so it has concave vertices
  poly <- circle_poly(100)
  # Nudge one vertex inward to create a concave dent
  poly$x[50] <- poly$x[50] * 0.5
  poly$y[50] <- poly$y[50] * 0.5
  result <- simplify_outer(poly, max_area = 0.1)
  expect_lt(length(result$x), length(poly$x))
})

poly_area <- function(p) {
  nx <- c(seq(2L, length(p$x)), 1L)
  abs(sum(p$x * p$y[nx] - p$x[nx] * p$y)) / 2
}
remainder_area <- function(rem) sum(vapply(rem, poly_area, numeric(1)))

test_that("simplified polygon encloses the original (moderate threshold)", {
  poly <- circle_poly(100)
  result <- simplify_outer(poly, max_area = 0.01)

  remainder <- polyclip(list(poly), list(result), "minus")
  # Any remainder must be floating-point noise, not real area
  expect_lt(remainder_area(remainder), 1e-6)
})

test_that("simplified polygon encloses the original (aggressive threshold)", {
  poly <- circle_poly(200)
  result <- simplify_outer(poly, max_area = 1, min_vertices = 4L)

  remainder <- polyclip(list(poly), list(result), "minus")
  expect_lt(remainder_area(remainder), 1e-6)
})

test_that("convex polygon is unchanged (no concave vertices to remove)", {
  # A convex polygon has no inward-dipping vertices; nothing can be removed
  # while guaranteeing enclosure, so the polygon is returned as-is.
  poly <- circle_poly(50)
  result <- simplify_outer(poly, max_area = Inf, min_vertices = 4L)
  expect_equal(length(result$x), length(poly$x))
})

test_that("max_area = 0 returns polygon unchanged", {
  poly <- circle_poly(20)
  result <- simplify_outer(poly, max_area = 0)
  expect_equal(result$x, poly$x)
  expect_equal(result$y, poly$y)
})

test_that("polygon at min_vertices is returned unchanged", {
  poly <- list(x = c(0, 1, 0.5), y = c(0, 0, 1))  # triangle
  result <- simplify_outer(poly, max_area = Inf, min_vertices = 3L)
  expect_equal(length(result$x), 3L)
})

test_that("scale invariance: simp_ratio gives same vertex count at any coordinate scale", {
  poly1 <- circle_poly(100, r = 1)
  poly2 <- list(x = poly1$x * 100, y = poly1$y * 100)

  simp_ratio <- 0.001

  # Threshold well below all triangle areas => nothing removed at either scale
  bbox_area1 <- diff(range(poly1$x)) * diff(range(poly1$y))
  bbox_area2 <- diff(range(poly2$x)) * diff(range(poly2$y))
  r1_low <- simplify_outer(poly1, max_area = 1e-12 * bbox_area1)
  r2_low <- simplify_outer(poly2, max_area = 1e-12 * bbox_area2)
  expect_equal(length(r1_low$x), length(r2_low$x))

  # Threshold well above all triangle areas => both reduced to min_vertices
  r1_high <- simplify_outer(poly1, max_area = 1e12 * bbox_area1)
  r2_high <- simplify_outer(poly2, max_area = 1e12 * bbox_area2)
  expect_equal(length(r1_high$x), length(r2_high$x))
})

# Test-only layout scoring (not shipped in the package).
#
# A layout is a placement table from `placeLabels()`: one row per label with the
# label box (`cxmin`/`cxmax`/`cymin`/`cymax`), the leader segment from its elbow
# (`ex`, `ey`) to its tip (`tx`, `ty`), and the leader length (`len`).

#' Test whether two line segments properly cross (vectorised)
#'
#' Touching at an endpoint or being parallel does not count as crossing.
#'
#' @param ax,ay,bx,by Start and end points of the first segments.
#' @param cx,cy,dx,dy Start and end points of the second segments.
#' @return A logical vector, `TRUE` where the two segments properly cross.
segmentsCross <- function(ax, ay, bx, by, cx, cy, dx, dy) {
  firstDx <- bx - ax
  firstDy <- by - ay
  secondDx <- dx - cx
  secondDy <- dy - cy

  denominator <- firstDx * secondDy - firstDy * secondDx
  alongFirst <- ((cx - ax) * secondDy - (cy - ay) * secondDx) / denominator
  alongSecond <- ((cx - ax) * firstDy - (cy - ay) * firstDx) / denominator

  notParallel <- abs(denominator) > 1e-12
  strictlyInside <- function(t) {
    t > 1e-9 & t < 1 - 1e-9
  }
  notParallel & strictlyInside(alongFirst) & strictlyInside(alongSecond)
}

#' Test whether a line segment meets an axis-aligned box (vectorised)
#'
#' @param ax,ay,bx,by Start and end points of the segments.
#' @param xmin,xmax,ymin,ymax Box bounds.
#' @return A logical vector, `TRUE` where segment and box intersect.
segmentCrossesBox <- function(ax, ay, bx, by, xmin, xmax, ymin, ymax) {
  endpointInside <- function(x, y) {
    x > xmin & x < xmax & y > ymin & y < ymax
  }
  crossesEdge <- function(x0, y0, x1, y1) {
    segmentsCross(ax, ay, bx, by, x0, y0, x1, y1)
  }

  endpointInside(ax, ay) |
    endpointInside(bx, by) |
    crossesEdge(xmin, ymin, xmax, ymin) |
    crossesEdge(xmax, ymin, xmax, ymax) |
    crossesEdge(xmax, ymax, xmin, ymax) |
    crossesEdge(xmin, ymax, xmin, ymin)
}

#' Score a label layout
#'
#' A feasible layout has all three conflict counts at zero.
#'
#' @param layout A placement table (see file header for the columns used).
#' @return A named numeric vector: `boxBox` (label boxes overlapping each
#'   other), `leaderLeader` (leaders crossing each other), `leaderBox` (a leader
#'   through another label's box), `totalLength` and `maxLength`.
layoutScore <- function(layout) {
  nLabels <- nrow(layout)
  if (nLabels < 2) {
    return(c(boxBox = 0, leaderLeader = 0, leaderBox = 0,
             totalLength = sum(layout$len),
             maxLength = if (nLabels > 0) max(layout$len) else 0))
  }

  pairs <- utils::combn(nLabels, 2)
  first <- pairs[1, ]
  second <- pairs[2, ]

  boxBox <- sum(
    layout$cxmin[first] < layout$cxmax[second] &
      layout$cxmin[second] < layout$cxmax[first] &
      layout$cymin[first] < layout$cymax[second] &
      layout$cymin[second] < layout$cymax[first]
  )

  # Leader-box conflicts are directional, so check both ways round.
  leaderThroughBox <- function(leader, box) {
    segmentCrossesBox(layout$ex[leader], layout$ey[leader],
                      layout$tx[leader], layout$ty[leader],
                      layout$cxmin[box], layout$cxmax[box],
                      layout$cymin[box], layout$cymax[box])
  }
  leaderBox <- sum(leaderThroughBox(first, second) |
                     leaderThroughBox(second, first))

  leaderLeader <- sum(
    segmentsCross(layout$ex[first], layout$ey[first],
                  layout$tx[first], layout$ty[first],
                  layout$ex[second], layout$ey[second],
                  layout$tx[second], layout$ty[second])
  )

  c(boxBox = boxBox, leaderLeader = leaderLeader, leaderBox = leaderBox,
    totalLength = sum(layout$len), maxLength = max(layout$len))
}

#' Assert that a layout has not regressed against a golden score
#'
#' Conflict counts are gated strictly; total leader length is allowed `tolerance`
#' fractional slack. `maxLength` is not gated.
#'
#' @param score Score of the layout under test, from `layoutScore()`.
#' @param golden Recorded golden score, in the same form.
#' @param tolerance Fractional slack allowed on `totalLength`.
#' @return Invisibly `score`.
expectNoWorseThan <- function(score, golden, tolerance) {
  for (conflict in c("boxBox", "leaderLeader", "leaderBox")) {
    expect_lte(unname(score[conflict]), unname(golden[conflict]))
  }
  expect_lte(unname(score["totalLength"]),
             unname(golden["totalLength"]) * (1 + tolerance))
  invisible(score)
}

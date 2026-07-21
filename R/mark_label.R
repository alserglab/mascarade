#' Enclosing polygon simplification
#'
#' Greedily removes small concave (inward) vertices from a polygon: a vertex is dropped when the
#' triangle it cuts off has area below `max_area`. Because only concave vertices are removed
#' the simplified polygon ENCLOSES the original, so the box-fit keep-out built from it stays
#' conservative. Used to cut vertex counts before placement.
#'
#' @param poly A list with numeric `x`, `y` (the polygon vertices).
#' @param max_area Numeric area threshold; vertices whose cut-off triangle is smaller are removed.
#' @param min_vertices Integer floor on the number of vertices kept.
#' @return A list with simplified numeric `x`, `y`.
#' @keywords internal
simplify_outer <- function(poly, max_area, min_vertices = 4L) {
  x <- poly$x
  y <- poly$y
  n <- length(x)
  if (n <= min_vertices || max_area <= 0) return(poly)

  # Circular linked list: prv[i] / nxt[i] are the neighbour indices of vertex i
  prv <- c(n, seq_len(n - 1L))
  nxt <- c(seq(2L, n), 1L)

  # Signed polygon area via shoelace: positive = CCW, negative = CW
  ccw <- sum(x * y[nxt] - x[nxt] * y) / 2 > 0

  # Vectorised initial cross products: (P_i - P_prev) x (P_next - P_prev)
  # Positive = left turn (convex in CCW), Negative = right turn (concave in CCW).
  crosses <- (x - x[prv]) * (y[nxt] - y[prv]) - (y - y[prv]) * (x[nxt] - x[prv])

  # For a CCW polygon a concave vertex (right turn, cross <= 0) dips inward.
  # Removing it fills the dent so the simplified polygon encloses the original.
  # Mark non-removable (convex) vertices NA so which.min() skips them; dead vertices are set NA
  # too as they are removed (see below).
  areas <- abs(crosses) / 2
  if (ccw) areas[crosses > 0] <- NA_real_ else areas[crosses < 0] <- NA_real_

  alive <- rep(TRUE, n)
  n_alive <- n
  repeat {
    if (n_alive <= min_vertices) break
    best_i <- which.min(areas)  # NA values are ignored by which.min
    if (length(best_i) == 0L || areas[best_i] > max_area) break

    p <- prv[best_i]; nx <- nxt[best_i]
    nxt[p] <- nx; prv[nx] <- p
    alive[best_i] <- FALSE
    areas[best_i] <- NA_real_     # exclude from future which.min (dead vertex)
    n_alive <- n_alive - 1L

    # Recompute cross products only for the two affected neighbours
    for (nb in c(p, nx)) {
      a_nb <- prv[nb]; b_nb <- nxt[nb]
      cr <- (x[nb] - x[a_nb]) * (y[b_nb] - y[a_nb]) -
            (y[nb] - y[a_nb]) * (x[b_nb] - x[a_nb])
      removable <- if (ccw) cr <= 0 else cr >= 0
      areas[nb] <- if (removable) abs(cr) / 2 else NA_real_
    }
  }

  cur <- which(alive)[1L]
  rx <- numeric(n_alive); ry <- numeric(n_alive)
  for (j in seq_len(n_alive)) { rx[j] <- x[cur]; ry[j] <- y[cur]; cur <- nxt[cur] }
  list(x = rx, y = ry)
}

#' Clean and simplify polygons for placement
#'
#' Drops non-finite vertices (draw-stage polygons can be cropped by axis limits). Degenerate
#' polygons are dropped upstream by `degeneratePolygon()` in `makeContent.shape_enc()`, so this
#' asserts at least three finite vertices remain rather than guarding the case. When
#' `simp_ratio > 0` it then removes small inward dents with `simplify_outer()` -- a big speed-up
#' (box-fit and foreignLength walk every edge), and since only concave vertices are removed each
#' polygon still ENCLOSES the original, so the box-fit keep-out stays conservative.
#'
#' @param polys List of polygons (`list(x, y)`), already subset to the drawn labels.
#' @param simp_ratio Numeric simplification fraction; `0` disables simplification.
#' @return A list of cleaned (and optionally simplified) polygons.
#' @keywords internal
#' @noRd
prepPolygons <- function(polys, simp_ratio) {
  cleaned <- lapply(polys, function(p) {
    finite <- is.finite(p$x) & is.finite(p$y)
    x <- p$x[finite]
    y <- p$y[finite]
    # degeneratePolygon() drops sub-3-vertex / zero-area polygons upstream, so the box-fit and
    # pole solvers never see a degenerate polygon here.
    stopifnot(length(x) >= 3)
    list(x = x, y = y)
  })
  if (simp_ratio > 0) {
    allX <- unlist(lapply(cleaned, `[[`, "x"))
    allY <- unlist(lapply(cleaned, `[[`, "y"))
    maxArea <- simp_ratio * diff(range(allX)) * diff(range(allY))
    cleaned <- lapply(cleaned, simplify_outer, max_area = maxArea)
  }
  cleaned
}

#' Find each polygon's pole of inaccessibility
#'
#' Returns the most-interior point of every polygon via `polylabelr::poi()`. Polygons reaching
#' here are non-degenerate (see `degeneratePolygon()`), so `poi()` is expected to succeed and
#' return a finite point; no fallback -- if it does not, that is a bug and should surface (the
#' caller's `stopifnot(all(is.finite(poles)))` catches a non-finite result).
#'
#' @param polys List of cleaned polygons (`list(x, y)`).
#' @return A two-column matrix of pole `x`, `y` coordinates, one row per polygon.
#' @keywords internal
#' @noRd
#' @importFrom polylabelr poi
polesOfInaccessibility <- function(polys) {
  t(vapply(polys, function(p) {
    pole <- polylabelr::poi(p$x, p$y)
    c(pole$x, pole$y)
  }, numeric(2)))
}

#' Warn when placed label boxes spill outside the panel
#'
#' The placer minimises viewport overflow but will still push a box off-panel (where it is
#' clipped) when there is no room left. Emits a single message naming how many labels overflowed.
#'
#' @param layout The placement table returned by `placeLabels()`.
#' @param bounds Numeric `c(width, height)` of the panel in mm.
#' @return `invisible(NULL)`; called for its warning side effect.
#' @keywords internal
#' @noRd
warnOnOverflow <- function(layout, bounds) {
  overflow <- pmax(0,
                   -(layout$cx - layout$hw), (layout$cx + layout$hw) - bounds[1],
                   -(layout$cy - layout$hh), (layout$cy + layout$hh) - bounds[2])
  # Count only real spills: 1e-3 mm (1 um) tolerance ignores floating-point noise from a box
  # that the placer parked flush against the panel edge.
  nOver <- sum(overflow > 1e-3)
  if (nOver == 0) {
    return(invisible(NULL))
  }
  cli::cli_warn(c(
    "!" = "{nOver} cluster label{?s} did not fully fit inside the plot area.",
    "i" = paste("Decrease {.arg label.fontsize} or expand the plot limits",
                "to make more room for label placement.")
  ))
}

#' Place labels at draw time (boundary-seed placer)
#'
#' Runs at draw time in the panel's millimetre space: builds the poles and the box-fit R-tree
#' from the cluster polygons (cheap, ~20 ms for ~40 clusters; the expensive mask is not
#' recomputed) and calls `placeLabels()`. The box-fit keep-out uses the dilated polygons while
#' poles, leader ends and foreign-routing use the true ones, so leaders reach the real outline.
#'
#' Note: each cluster polygon here is a single `list(x, y)` ring; any mask holes are resolved
#' upstream in `generateMask()`, so this layer treats every polygon as one simple ring.
#'
#' @param rects List of measured label box sizes `c(w, h)` in mm (a zeroed entry = not drawn).
#' @param polygons List of true cluster polygons (`list(x, y)`) in mm.
#' @param polygons_pad List of the same polygons dilated by `label.buffer` (the box keep-out).
#' @param bounds Numeric `c(width, height)` of the panel in mm.
#' @param simp_ratio Numeric polygon-simplification fraction (see `simplify_outer()`).
#' @param con_type Leader style: `"ledge"`, `"line"`, `"box"`, or `"none"`.
#' @param buffer Numeric `label.buffer` in mm; the overflow viewport is inset by it.
#' @param hardpad Numeric `label.hardpad` in mm: hard box clearance folded into every placement
#'   rectangle (seed slots, sweeps and polish alike).
#' @param softpad Numeric `label.softpad` in mm: extra target box spacing the polish aims for, on
#'   top of `hardpad`.
#' @return A list, one entry per input label: the placed centre `c(x, y)` in mm (`NULL` if not
#'   drawn), carrying `attr(., "leaders")` with `c(ex, ey, bx, by, corner)` per drawn label.
#' @keywords internal
#' @importFrom stats median
my_place_labels <- function(rects, polygons, polygons_pad, bounds,
                            simp_ratio = 0.001, con_type = "ledge", buffer = 0,
                            hardpad = 0, softpad = 0) {
  nLabels <- length(rects)

  # Assemble the return value: the centres list, carrying per-label leaders in an attribute.
  # Each leader is c(ex, ey, bx, by, corner): leader start -> visible end (mask boundary) plus
  # the ledge flag, consumed by the drawing code.
  withLeaders <- function(centres, leaders) {
    attr(centres, "leaders") <- leaders
    centres
  }

  # A label is drawn only if its measured box is non-zero.
  drawn <- which(vapply(rects, function(r) !all(r == 0), logical(1)))
  if (length(drawn) == 0) {
    return(withLeaders(vector("list", nLabels), vector("list", nLabels)))
  }

  # When placement genuinely cannot run -- a zero-size panel, or the C++ solver failing -- warn
  # once and draw no labels. (Degenerate polygons are dropped upstream, so they never reach here.)
  skipPlacement <- function(reason) {
    cli::cli_warn(c(
      "!" = "Could not place cluster labels: {reason}.",
      "i" = "No labels were drawn."
    ))
    withLeaders(vector("list", nLabels), vector("list", nLabels))
  }

  # True polygons drive poles, leader ends and foreign routing; padded polygons (dilated by
  # label.buffer) drive the box-fit keep-out, so leaders still reach the real cluster outline.
  truePolys <- prepPolygons(polygons[drawn], simp_ratio)
  paddedPolys <- prepPolygons(polygons_pad[drawn], simp_ratio)
  poles <- polesOfInaccessibility(truePolys)

  # Poles are finite by construction (>= 3-vertex polygons; polesOfInaccessibility falls back to
  # the centroid), and the panel size (convertWidth/Height of 1 npc during an active draw) is
  # always finite. A collapsed / zero-size panel is still possible, so guard that (the finiteness
  # assert also keeps the `<= 0` comparison NA-safe).
  stopifnot(all(is.finite(poles)), all(is.finite(bounds)))
  if (bounds[1] <= 0 || bounds[2] <= 0) {
    return(skipPlacement("the plot panel has zero size"))
  }

  centres <- vector("list", nLabels)
  leaders <- vector("list", nLabels)

  halfWidth <- vapply(drawn, function(i) rects[[i]][1] / 2, 0)
  halfHeight <- vapply(drawn, function(i) rects[[i]][2] / 2, 0)
  charHeight <- stats::median(2 * halfHeight)

  truePolysX <- lapply(truePolys, `[[`, "x")
  truePolysY <- lapply(truePolys, `[[`, "y")
  paddedPolysX <- lapply(paddedPolys, `[[`, "x")
  paddedPolysY <- lapply(paddedPolys, `[[`, "y")
  geom <- list(
    poi = poles,
    rtree = buildBoxFit(paddedPolysX, paddedPolysY),
    polysx = truePolysX,
    polysy = truePolysY,
    pad_xrange = range(unlist(paddedPolysX))
  )

  # Inset the overflow viewport by label.buffer, so labels keep the same gap from the panel
  # edge that they keep from clusters (overflow is measured against xhi - buffer, etc.).
  viewportInset <- min(buffer, 0.4 * min(bounds))
  layout <- tryCatch(
    placeLabels(geom,
                c(viewportInset, bounds[1] - viewportInset),
                c(viewportInset, bounds[2] - viewportInset),
                halfWidth, halfHeight, charHeight, con_type = con_type,
                hardPad = hardpad, softPad = softpad),
    error = function(e) NULL
  )
  if (is.null(layout)) {
    return(skipPlacement("the label placement solver failed"))
  }
  layout <- layout[order(layout$label)]

  warnOnOverflow(layout, bounds)

  for (j in seq_along(drawn)) {
    centres[[drawn[j]]] <- c(layout$cx[j], layout$cy[j])
    leaders[[drawn[j]]] <- c(layout$ex[j], layout$ey[j],
                             layout$bx[j], layout$by[j],
                             as.numeric(layout$corner[j]))
  }
  withLeaders(centres, leaders)
}

#' Compute a polygon's absolute area (shoelace formula)
#'
#' @param p A polygon (`list(x, y)`).
#' @return The non-negative area enclosed by the polygon.
#' @keywords internal
#' @noRd
polygonArea <- function(p) {
  abs(sum(p$x * c(p$y[-1], p$y[1]) - c(p$x[-1], p$x[1]) * p$y)) / 2
}

#' Is a polygon too degenerate to place a label against?
#'
#' A polygon is degenerate when the pole / box-fit solvers have nothing to work with: fewer than
#' three finite vertices (a point or a line, possibly after NA cropping), or three-plus vertices
#' enclosing a ~zero area (all-collinear or a repeated-vertex sliver). `generateMask()` never
#' emits such polygons, but a real cluster can still collapse this way after axis-limit cropping or
#' a negative `expand`.
#'
#' In principle these could be *supported* rather than dropped: `prepPolygons()` already
#' substitutes an eps-triangle for a sub-three-point polygon, and a point/line has a well-defined
#' anchor (its centroid / midpoint). The caller drops them for now because the supported input
#' never produces them.
#'
#' @param p A polygon (`list(x, y)`) in mm.
#' @param area_eps Minimum enclosed area (mm^2) a non-degenerate polygon must exceed. `1e-6` is a
#'   1 um square -- conservative enough to catch only genuinely collapsed polygons.
#' @return `TRUE` if the polygon is degenerate.
#' @keywords internal
#' @noRd
degeneratePolygon <- function(p, area_eps = 1e-6) {
  finite <- is.finite(p$x) & is.finite(p$y)
  x <- p$x[finite]
  y <- p$y[finite]
  if (length(x) < 3) {
    return(TRUE)
  }
  polygonArea(list(x = x, y = y)) < area_eps
}

#' Dilate each cluster polygon individually by `delta` mm (the label keep-out)
#'
#' Offsets one polygon at a time: `polyoffset()` on the whole list reorders (and can merge) its
#' output, which would break the 1:1 mapping the placer relies on (polygon `i` pairs with
#' `rects[i]`). `delta` (`label.buffer`) is non-negative, so dilating a single connected polygon
#' always yields exactly one connected piece (it can neither vanish nor split -- those need a
#' negative offset), which the assertion enforces. Used only for the box-fit keep-out -- poles,
#' leader ends and foreign routing use the true polygons.
#'
#' @param polygons List of cluster polygons (`list(x, y)`) in mm.
#' @param delta Dilation distance in mm (`>= 0`).
#' @return A list of dilated polygons, aligned 1:1 with `polygons`.
#' @keywords internal
#' @noRd
#' @importFrom polyclip polyoffset
dilatePolygons <- function(polygons, delta) {
  lapply(polygons, function(p) {
    dilated <- polyoffset(list(p), delta)
    stopifnot(length(dilated) == 1)
    dilated[[1]]
  })
}

#' Build the connector grob for the drawn labels
#'
#' Each drawn label contributes zero or more 2-point line segments: the leader (box anchor ->
#' visible mask-boundary end, trimmed to leave a `con_cap` gap at the cluster), the horizontal
#' ledge along the box edge for `con_type == "ledge"`, and the four box-outline edges for
#' `con_type == "box"`. All segments are drawn as one `polylineGrob`, one id per segment, with
#' the connector `gp` subset to each segment's label.
#'
#' @param labels_drawn Integer indices (into the full label list) of the drawn labels.
#' @param leaders Per-label `c(ex, ey, bx, by, corner)` leader descriptors.
#' @param labelpos Per-label placed centres `c(x, y)`.
#' @param dims List of measured label box sizes `c(w, h)` in mm.
#' @param con_type Leader style: `"ledge"`, `"line"`, `"box"`, or `"none"`.
#' @param con_cap Numeric gap (mm) left between the leader end and the cluster.
#' @param con_gp A `gpar` for the connectors (per drawn label).
#' @param arrow Optional `grid::arrow` for the connectors.
#' @return A single `polylineGrob`, or `nullGrob()` when there is nothing to draw.
#' @keywords internal
#' @noRd
#' @importFrom grid nullGrob polylineGrob
buildConnectorGrob <- function(labels_drawn, leaders, labelpos, dims,
                               con_type, con_cap, con_gp, arrow) {
  if (con_type == "none") {
    return(nullGrob())
  }

  # The 2-point segments contributed by one drawn label (index `i` into `labels_drawn`).
  segmentsFor <- function(i) {
    idx <- labels_drawn[i]
    leader <- leaders[[idx]]
    centre <- labelpos[[idx]]
    halfWidth <- dims[[idx]][1] / 2
    halfHeight <- dims[[idx]][2] / 2
    segs <- list()

    # Leader: box anchor c(ex, ey) -> visible end c(bx, by). A degenerate zero-length leader
    # (start == end) falls below the 1e-4 mm floor and draws nothing.
    x0 <- leader[1]
    y0 <- leader[2]
    x1 <- leader[3]
    y1 <- leader[4]
    leaderLen <- sqrt((x1 - x0)^2 + (y1 - y0)^2)
    if (leaderLen > 1e-4) {
      if (con_cap > 0 && leaderLen > con_cap) {
        # Stop the leader con_cap mm short of the cluster.
        shrink <- (leaderLen - con_cap) / leaderLen
        x1 <- x0 + (x1 - x0) * shrink
        y1 <- y0 + (y1 - y0) * shrink
      }
      segs <- c(segs, list(list(x = c(x0, x1), y = c(y0, y1), label = i)))
    }

    # Ledge: a short horizontal line along the box edge at the anchor's y.
    if (con_type == "ledge" && length(leader) >= 5 && leader[5] == 1) {
      segs <- c(segs, list(list(
        x = c(centre[1] - halfWidth, centre[1] + halfWidth),
        y = c(y0, y0), label = i
      )))
    }

    # Box: outline the label's bounding box as four edges (corners in a closed loop).
    if (con_type == "box") {
      cornerX <- centre[1] + c(-halfWidth, halfWidth, halfWidth, -halfWidth, -halfWidth)
      cornerY <- centre[2] + c(-halfHeight, -halfHeight, halfHeight, halfHeight, -halfHeight)
      for (e in seq_len(4L)) {
        segs <- c(segs, list(list(
          x = cornerX[c(e, e + 1L)], y = cornerY[c(e, e + 1L)], label = i
        )))
      }
    }
    segs
  }

  segments <- do.call(c, lapply(seq_along(labels_drawn), segmentsFor))
  if (length(segments) == 0) {
    return(nullGrob())
  }
  if (!is.null(arrow)) {
    arrow$ends <- 2L
  }
  labelIndex <- vapply(segments, `[[`, integer(1), "label")
  gp <- subset_gp(subset_gp(con_gp, labels_drawn), labelIndex)
  polylineGrob(
    x = unlist(lapply(segments, `[[`, "x")),
    y = unlist(lapply(segments, `[[`, "y")),
    id = rep(seq_along(segments), each = 2L),
    default.units = "mm", gp = gp, arrow = arrow
  )
}

#' Build the label + leader grobs for a mark
#'
#' Draw-time worker for `makeContent.shape_enc()`: dilates the cluster polygons by `buffer`
#' (the box keep-out), calls `my_place_labels()` for the placement, positions the label box
#' grobs and builds the leader polylines (anchor -> visible mask-boundary end, plus the
#' horizontal ledge for `con_type == "ledge"`).
#'
#' @param labels List of label-box grobs (one per mark part).
#' @param dims List of measured label box sizes `c(w, h)` in mm.
#' @param polygons List of cluster polygons (`list(x, y)`) in mm.
#' @param ghosts Points to avoid (currently unused by the placer).
#' @param buffer Grid unit: the `label.buffer` polygon padding / box keep-out.
#' @param con_type Leader style: `"ledge"`, `"line"`, `"box"`, or `"none"`.
#' @param con_cap Numeric gap (mm) left between the leader end and the cluster.
#' @param con_gp A `gpar` for the connectors (per drawn label).
#' @param arrow Optional `grid::arrow` for the connectors.
#' @param simp_ratio Numeric polygon-simplification fraction (see `simplify_outer()`).
#' @param hardpad Grid unit: the `label.hardpad` hard box clearance.
#' @param softpad Grid unit: the `label.softpad` extra polish-only target box spacing.
#' @return A `gList`-ready list: the positioned label grobs followed by the connector grob.
#' @keywords internal
#' @importFrom grid convertWidth convertHeight nullGrob unit
my_make_label <- function(labels, dims, polygons, ghosts, buffer, con_type,
                          con_cap, con_gp, arrow, simp_ratio = 0.001,
                          hardpad = unit(0, "pt"), softpad = unit(0, "pt")) {
  # `label.buffer` keep-out: dilate each cluster by `buffer` (mm). Used only for the box-fit
  # keep-out; poles, leader ends and foreign routing use the true `polygons`, so leaders reach
  # the actual cluster outline. label.buffer is a padding, so it must be non-negative -- a
  # negative offset would erode (and could split or empty) the clusters.
  delta <- convertWidth(buffer, "mm", TRUE)
  stopifnot(delta >= 0)
  polygons_pad <- dilatePolygons(polygons, delta)

  # Hard/soft box clearances (mm). Both are non-negative padding amounts.
  hardpad_mm <- convertWidth(hardpad, "mm", TRUE)
  softpad_mm <- convertWidth(softpad, "mm", TRUE)
  stopifnot(hardpad_mm >= 0, softpad_mm >= 0)

  panel <- c(
    convertWidth(unit(1, "npc"), "mm", TRUE),
    convertHeight(unit(1, "npc"), "mm", TRUE)
  )
  labelpos <- my_place_labels(dims, polygons, polygons_pad, panel,
                              simp_ratio = simp_ratio, con_type = con_type,
                              buffer = delta,
                              hardpad = hardpad_mm, softpad = softpad_mm)
  if (all(lengths(labelpos) == 0)) {
    return(list(nullGrob()))
  }

  labels_drawn <- which(!vapply(labelpos, is.null, logical(1)))
  leaders <- attr(labelpos, "leaders")

  # Position each label-box grob at its placed centre; undrawn entries become nullGrobs.
  labels <- Map(function(lab, pos) {
    if (is.null(pos) || inherits(lab, "null")) {
      return(nullGrob())
    }
    lab$vp$x <- unit(pos[1], "mm")
    lab$vp$y <- unit(pos[2], "mm")
    lab
  }, lab = labels, pos = labelpos)

  connect <- buildConnectorGrob(labels_drawn, leaders, labelpos, dims,
                                con_type, con_cap, con_gp, arrow)
  c(labels, list(connect))
}


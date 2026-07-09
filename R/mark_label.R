#' Enclosing polygon simplification
#'
#' Greedily removes small concave (inward) vertices from a ring: a vertex is dropped when the
#' triangle it cuts off has area below `max_area`. Because only concave vertices are removed
#' the simplified ring ENCLOSES the original, so the box-fit keep-out built from it stays
#' conservative. Used to cut vertex counts before placement.
#'
#' @param poly A list with numeric `x`, `y` (the ring vertices).
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
  # Mark non-removable vertices with Inf so which.min() skips them.
  areas <- abs(crosses) / 2
  # Mark non-removable vertices with NA so which.min() skips them;
  # dead vertices (removed) are marked Inf.
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

#' Draw-time label placement (boundary-seed placer)
#'
#' Runs at draw time in the panel's millimetre space: builds the poles and the box-fit R-tree
#' from the cluster rings (cheap, ~20 ms for ~40 clusters; the expensive mask is not
#' recomputed) and calls `placeLabels()`. The box-fit keep-out uses the dilated polygons while
#' poles, leader ends and foreign-routing use the true ones, so leaders reach the real outline.
#'
#' @param rects List of measured label box sizes `c(w, h)` in mm (a zeroed entry = not drawn).
#' @param polygons List of true cluster rings (`list(x, y)`) in mm.
#' @param polygons_pad List of the same rings dilated by `label.buffer` (the box keep-out).
#' @param bounds Numeric `c(width, height)` of the panel in mm.
#' @param anchors List of fallback anchor points, used only for degenerate input.
#' @param simp_ratio Numeric polygon-simplification fraction (see `simplify_outer()`).
#' @param con_type Leader style: `"cl"`, `"cm"`, or `"none"`.
#' @param buffer Numeric `label.buffer` in mm; the overflow viewport is inset by it.
#' @return A list, one entry per input label: the placed centre `c(x, y)` in mm (`NULL` if not
#'   drawn), carrying `attr(., "leaders")` with `c(ex, ey, bx, by, corner)` per drawn label.
#' @keywords internal
#' @importFrom polylabelr poi
#' @importFrom stats median
my_place_labels <- function(rects, polygons, polygons_pad, bounds, anchors,
                            simp_ratio = 0.001, con_type = "cl", buffer = 0) {
  res <- vector('list', length(rects))    # label centres (mm)
  lead <- vector('list', length(rects))   # per label c(ex, ey, bx, by, corner): leader start ->
                                          # visible end (mask boundary) + ledge flag, for drawing
  withLeaders <- function() { attr(res, "leaders") <- lead; res }
  active <- which(vapply(rects, function(r) !all(r == 0), logical(1)))
  if (length(active) == 0) return(withLeaders())

  # Clean + simplify a polygon list (aligned to `active`). Clean: drop non-finite vertices and
  # guarantee >= 3 points so the C++ box-fit / pole solver never see a degenerate ring (draw-
  # stage polygons can be cropped by axis limits, dropped by expansion, or reduced to a point).
  # Simplify: drop small inward dents -- a big speed-up (box-fit and foreignLength walk every
  # edge), and since simplify_outer only removes concave vertices the ring ENCLOSES the
  # original, so the box-fit keep-out stays conservative.
  prep <- function(polys) {
    pa <- lapply(polys[active], function(p) {
      ok <- is.finite(p$x) & is.finite(p$y); x <- p$x[ok]; y <- p$y[ok]
      if (length(x) < 3) {
        cx <- if (length(x)) mean(x) else 0; cy <- if (length(y)) mean(y) else 0
        x <- cx + c(-1e-3, 1e-3, 0); y <- cy + c(-1e-3, -1e-3, 1e-3)
      }
      list(x = x, y = y)
    })
    if (simp_ratio > 0) {
      ax <- unlist(lapply(pa, `[[`, 'x')); ay <- unlist(lapply(pa, `[[`, 'y'))
      pa <- lapply(pa, simplify_outer, max_area = simp_ratio * diff(range(ax)) * diff(range(ay)))
    }
    pa
  }
  polys_t <- prep(polygons)               # true clusters: poles, leader ends, foreign routing
  polys_d <- prep(polygons_pad)           # dilated by label.buffer: box-fit keep-out only
  poimat <- t(vapply(polys_t, function(p) {
    pl <- tryCatch(polylabelr::poi(p$x, p$y), error = function(e) NULL)
    if (is.null(pl) || !is.finite(pl$x) || !is.finite(pl$y)) c(mean(p$x), mean(p$y)) else c(pl$x, pl$y)
  }, numeric(2)))

  # anchor fallback if anything is still degenerate (start == end => no visible leader drawn)
  anchorFallback <- function() {
    for (j in seq_along(active)) {
      a <- as.numeric(anchors[[active[j]]])
      res[[active[j]]] <<- a; lead[[active[j]]] <<- c(a, a, 0)
    }
    withLeaders()
  }
  if (any(!is.finite(poimat)) || any(!is.finite(bounds)) || bounds[1] <= 0 || bounds[2] <= 0)
    return(anchorFallback())
  if (length(active) == 1) {
    res[[active]] <- poimat[1, ]; lead[[active]] <- c(poimat[1, ], poimat[1, ], 0)
    return(withLeaders())
  }

  hw <- vapply(active, function(i) rects[[i]][1] / 2, 0)
  hh <- vapply(active, function(i) rects[[i]][2] / 2, 0)
  char_h <- stats::median(2 * hh)
  pxt <- lapply(polys_t, `[[`, 'x'); pyt <- lapply(polys_t, `[[`, 'y')
  pxd <- lapply(polys_d, `[[`, 'x'); pyd <- lapply(polys_d, `[[`, 'y')
  # box-fit keep-out from the dilated polygons; poles/foreign/leader-ends from the true ones,
  # so leaders are drawn to the actual cluster outline. Seed columns use the dilated extent.
  geom <- list(poi = poimat, rtree = buildBoxFit(pxd, pyd), polysx = pxt, polysy = pyt,
               pad_xrange = range(unlist(pxd)))

  # Inset the overflow viewport by label.buffer, so labels keep the same gap from the panel
  # edge that they keep from clusters (overflow is measured against xhi - buffer, etc.).
  vb <- min(buffer, 0.4 * min(bounds))
  lay <- tryCatch(placeLabels(geom, c(vb, bounds[1] - vb), c(vb, bounds[2] - vb), hw, hh, char_h,
                              con_type = con_type),
                  error = function(e) NULL)
  if (is.null(lay)) return(anchorFallback())
  lay <- lay[order(lay$label)]
  for (j in seq_along(active)) {
    res[[active[j]]] <- c(lay$cx[j], lay$cy[j])
    lead[[active[j]]] <- c(lay$ex[j], lay$ey[j], lay$bx[j], lay$by[j], as.numeric(lay$corner[j]))
  }
  withLeaders()
}
#' Build the label + leader grobs for a mark
#'
#' Draw-time worker for `makeContent.shape_enc()`: dilates the cluster polygons by `buffer`
#' (the box keep-out), calls `my_place_labels()` for the placement, positions the label box
#' grobs and builds the leader polylines (anchor -> visible mask-boundary end, plus the
#' horizontal ledge for `con_type == "cl"`).
#'
#' @param labels List of label-box grobs (one per mark part).
#' @param dims List of measured label box sizes `c(w, h)` in mm.
#' @param polygons List of cluster rings (`list(x, y)`) in mm.
#' @param ghosts Points to avoid (currently unused by the placer).
#' @param buffer Grid unit: the `label.buffer` polygon padding / box keep-out.
#' @param con_type Leader style: `"cl"`, `"cm"`, or `"none"`.
#' @param con_cap Numeric gap (mm) left between the leader end and the cluster.
#' @param con_gp A `gpar` for the connectors (per drawn label).
#' @param anchor_x,anchor_y Optional per-label anchor overrides.
#' @param arrow Optional `grid::arrow` for the connectors.
#' @param simp_ratio Numeric polygon-simplification fraction (see `simplify_outer()`).
#' @return A `gList`-ready list: the positioned label grobs followed by the connector grob.
#' @keywords internal
#' @importFrom polyclip polyoffset
#' @importFrom grid convertWidth convertHeight nullGrob polylineGrob
#' @importFrom stats runif
my_make_label <- function(labels, dims, polygons, ghosts, buffer, con_type,
                       con_cap, con_gp, anchor_x,
                       anchor_y, arrow, simp_ratio = 0.001) {
  polygons <- lapply(polygons, function(p) {
    if (length(p$x) == 1 & length(p$y) == 1) {
      list(
        x = runif(200, p$x-0.00005, p$x+0.00005),
        y = runif(200, p$y-0.00005, p$y+0.00005)
      )
    } else {
      list(
        x = p$x,
        y = p$y
      )
    }
  })

  anchors <- lapply(seq_along(polygons), function(i) {
    x <- mean(range(polygons[[i]]$x))
    if (length(anchor_x) == length(polygons) && !is.na(anchor_x[i])) x <- anchor_x[i]
    y <- mean(range(polygons[[i]]$y))
    if (length(anchor_y) == length(polygons) && !is.na(anchor_y[i])) y <- anchor_y[i]
    c(x, y)
  })
  # `label.buffer` polygon padding: dilate each cluster by `buffer` to form the label
  # keep-out zone. Offset each polygon INDIVIDUALLY -- polyoffset() on the whole list reorders
  # (and can merge) its output, breaking the 1:1 mapping the placer relies on (polygon i is
  # paired with rects[i]/anchors[i]). The dilated set is used only for the box-fit keep-out;
  # poles, leader ends and foreign-routing use the true `polygons` so leaders reach the
  # actual cluster outline.
  delta <- convertWidth(buffer, 'mm', TRUE)
  p_big <- lapply(polygons, function(p) {
    off <- polyoffset(list(p), delta)
    if (length(off) == 0) return(p)          # degenerate offset: keep the original ring
    if (length(off) == 1) return(off[[1]])
    # a ring can split into pieces; keep the largest so exactly one ring maps to this label
    a <- vapply(off, function(q)
      abs(sum(q$x * c(q$y[-1], q$y[1]) - c(q$x[-1], q$x[1]) * q$y)) / 2, numeric(1))
    off[[which.max(a)]]
  })

  area <- c(
    convertWidth(unit(1, 'npc'), 'mm', TRUE),
    convertHeight(unit(1, 'npc'), 'mm', TRUE)
  )
  labelpos <- my_place_labels(dims, polygons, p_big, area, anchors, simp_ratio = simp_ratio,
                              con_type = con_type, buffer = delta)
  if (all(lengths(labelpos) == 0)) {
    return(list(nullGrob()))
  }
  labels_drawn <- which(!vapply(labelpos, is.null, logical(1)))
  leaders <- attr(labelpos, "leaders")
  labels <- Map(function(lab, pos) {
    if (is.null(pos) || inherits(lab, 'null')) return(nullGrob())
    lab$vp$x <- unit(pos[1], 'mm')
    lab$vp$y <- unit(pos[2], 'mm')
    lab
  }, lab = labels, pos = labelpos)
  # Draw each leader as the placer scored it: from the box anchor c(ex,ey) to the visible end
  # c(bx,by) = the first mask boundary along anchor->pole (the part inside the cluster is
  # hidden). For "cl" also draw the horizontal ledge along the box edge at the anchor's y.
  # Each drawn line is one polyline id; `gi` maps it back to its label for the connector gp.
  if (con_type == 'none') {
    connect <- nullGrob()
  } else {
    xs <- list(); ys <- list(); gi <- integer(0); k <- 0L
    for (i in seq_along(labels_drawn)) {
      idx <- labels_drawn[i]; l <- leaders[[idx]]; ctr <- labelpos[[idx]]
      x0 <- l[1]; y0 <- l[2]; x1 <- l[3]; y1 <- l[4]
      d <- sqrt((x1 - x0)^2 + (y1 - y0)^2)
      if (d > 1e-4) {                                          # visible leader present
        if (con_cap > 0 && d > con_cap) {                      # leave con.cap gap at the cluster
          f <- (d - con_cap) / d; x1 <- x0 + (x1 - x0) * f; y1 <- y0 + (y1 - y0) * f
        }
        k <- k + 1L; xs[[k]] <- c(x0, x1); ys[[k]] <- c(y0, y1); gi[k] <- i
      }
      if (con_type == 'cl' && length(l) >= 5 && l[5] == 1) {    # ledge = box edge at anchor y
        hw_i <- dims[[idx]][1] / 2
        k <- k + 1L; xs[[k]] <- c(ctr[1] - hw_i, ctr[1] + hw_i); ys[[k]] <- c(l[2], l[2]); gi[k] <- i
      }
    }
    if (k == 0) {
      connect <- nullGrob()
    } else {
      if (!is.null(arrow)) arrow$ends <- 2L
      gp <- subset_gp(subset_gp(con_gp, labels_drawn), gi)
      connect <- polylineGrob(unlist(xs), unlist(ys), id = rep(seq_len(k), each = 2),
                              default.units = 'mm', gp = gp, arrow = arrow)
    }
  }
  c(labels, list(connect))
}


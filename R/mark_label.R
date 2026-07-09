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

# Label placement (mascarade boundary-seed placer). Runs at draw time in the panel's
# millimetre space, where `rects` are the measured label box sizes (w, h in mm), `polygons`
# are the cluster rings in mm, and `bounds` is the panel size in mm. Returns, per label, its
# placed centre c(x, y) in mm (NULL for labels not drawn). Poles and the box-fit R-tree are
# recomputed here each draw (cheap: ~20 ms for ~40 clusters); the expensive mask is not.
# `anchors` is used only as the degenerate-input fallback; non-labelled polygon parts are
# treated as the placed set (points to avoid / ghosts are not handled yet). `simp_ratio`
# simplifies the polygons used for placement (box-fit + foreign-crossing) — see below.
#' @importFrom polylabelr poi
#' @importFrom stats median
my_place_labels <- function(rects, polygons, bounds, anchors, simp_ratio = 0.001) {
  res <- vector('list', length(rects))    # label centres (mm)
  pol <- vector('list', length(rects))    # matching poles = leader targets (mm), for the drawn connector
  withPoles <- function() { attr(res, "poles") <- pol; res }
  active <- which(vapply(rects, function(r) !all(r == 0), logical(1)))
  if (length(active) == 0) return(withPoles())

  # clean each polygon: drop non-finite vertices; guarantee >= 3 points so the C++ box-fit
  # and pole solver never receive an empty/degenerate ring (draw-stage polygons can be
  # cropped by axis limits, dropped by expansion, or reduced to a point).
  polys_a <- lapply(polygons[active], function(p) {
    ok <- is.finite(p$x) & is.finite(p$y); x <- p$x[ok]; y <- p$y[ok]
    if (length(x) < 3) {
      cx <- if (length(x)) mean(x) else 0; cy <- if (length(y)) mean(y) else 0
      x <- cx + c(-1e-3, 1e-3, 0); y <- cy + c(-1e-3, -1e-3, 1e-3)
    }
    list(x = x, y = y)
  })

  # Simplify the placement polygons (drop small inward dents). This is a big speed-up:
  # box-fit `intersects(box, polygon)` and `foreignLength` both walk every polygon edge,
  # and mask rings can have hundreds of vertices. Because `simplify_outer` only removes
  # concave vertices the simplified ring ENCLOSES the original, so box-fit stays
  # conservative (never lets a box onto the real cluster). Placement only; the drawn
  # connectors still project onto the full-resolution polygons in my_make_label().
  if (simp_ratio > 0 && length(active) > 0) {
    ax <- unlist(lapply(polys_a, `[[`, 'x')); ay <- unlist(lapply(polys_a, `[[`, 'y'))
    max_area <- simp_ratio * diff(range(ax)) * diff(range(ay))
    polys_a <- lapply(polys_a, simplify_outer, max_area = max_area)
  }
  poimat <- t(vapply(polys_a, function(p) {
    pl <- tryCatch(polylabelr::poi(p$x, p$y), error = function(e) NULL)
    if (is.null(pl) || !is.finite(pl$x) || !is.finite(pl$y)) c(mean(p$x), mean(p$y)) else c(pl$x, pl$y)
  }, numeric(2)))

  # anchor fallback if anything is still degenerate (pole = centre => no visible leader)
  anchorFallback <- function() {
    for (j in seq_along(active)) {
      a <- as.numeric(anchors[[active[j]]]); res[[active[j]]] <<- a; pol[[active[j]]] <<- a
    }
    withPoles()
  }
  if (any(!is.finite(poimat)) || any(!is.finite(bounds)) || bounds[1] <= 0 || bounds[2] <= 0)
    return(anchorFallback())
  if (length(active) == 1) { res[[active]] <- poimat[1, ]; pol[[active]] <- poimat[1, ]; return(withPoles()) }

  hw <- vapply(active, function(i) rects[[i]][1] / 2, 0)
  hh <- vapply(active, function(i) rects[[i]][2] / 2, 0)
  char_h <- stats::median(2 * hh)
  px <- lapply(polys_a, `[[`, 'x'); py <- lapply(polys_a, `[[`, 'y')
  geom <- list(poi = poimat, rtree = buildBoxFit(px, py), polysx = px, polysy = py)

  lay <- tryCatch(placeLabels(geom, c(0, bounds[1]), c(0, bounds[2]), hw, hh, char_h),
                  error = function(e) NULL)
  if (is.null(lay)) return(anchorFallback())
  lay <- lay[order(lay$label)]
  for (j in seq_along(active)) {
    res[[active[j]]] <- c(lay$cx[j], lay$cy[j]); pol[[active[j]]] <- poimat[j, ]
  }
  withPoles()
}
#' @importFrom polyclip polyoffset
#' @importFrom grid convertWidth convertHeight nullGrob polylineGrob
#' @importFrom stats runif
my_make_label <- function(labels, dims, polygons, ghosts, buffer, con_type,
                       con_border, con_cap, con_gp, anchor_mod, anchor_x,
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
  p_big <- polyoffset(polygons, convertWidth(buffer, 'mm', TRUE))

  area <- c(
    convertWidth(unit(1, 'npc'), 'mm', TRUE),
    convertHeight(unit(1, 'npc'), 'mm', TRUE)
  )
  labelpos <- my_place_labels(dims, p_big, area, anchors, simp_ratio = simp_ratio)
  if (all(lengths(labelpos) == 0)) {
    return(list(nullGrob()))
  }
  labels_drawn <- which(!vapply(labelpos, is.null, logical(1)))
  labels <- Map(function(lab, pos) {
    if (is.null(pos) || inherits(lab, 'null')) return(nullGrob())
    lab$vp$x <- unit(pos[1], 'mm')
    lab$vp$y <- unit(pos[2], 'mm')
    lab
  }, lab = labels, pos = labelpos)
  # Leader target = the cluster POLE the placer optimised each leader against, so the drawn
  # connector (box-edge -> pole, straight() clips it at the box) matches the placement
  # algorithm exactly. (Previously this projected onto the nearest polygon boundary point,
  # which is a different endpoint and could look inconsistent with the optimised layout.)
  # NOTE (future): the placement algorithm assumes STRAIGHT leaders; if we ever want the
  # connectors drawn as elbows to also be conflict-free, the elbow geometry would have to be
  # modelled inside the placer (box -> bend -> pole) -- worth exposing as a parameter then.
  connect <- rlang::inject(rbind(!!!attr(labelpos, "poles")[lengths(labelpos) != 0]))
  labeldims <- rlang::inject(rbind(!!!dims[lengths(labelpos) != 0])) / 2
  labelpos <- rlang::inject(rbind(!!!labelpos))
  if (con_type == 'none' || !con_type %in% c('elbow', 'straight')) {
    connect <- nullGrob()
  } else {
    con_fun <- switch(con_type, elbow = elbow, straight = straight)
    connect <- con_fun(
      labelpos[, 1] - labeldims[, 1], labelpos[, 1] + labeldims[, 1],
      labelpos[, 2] - labeldims[, 2], labelpos[, 2] + labeldims[, 2],
      connect[, 1], connect[, 2]
    )
    if (con_border == 'one') {
      connect <- with_borderline(
        labelpos[, 1] - labeldims[, 1],
        labelpos[, 1] + labeldims[, 1], connect
      )
    }
    connect <- end_cap(connect, con_cap)
    connect <- zip_points(connect)
    if (!is.null(arrow)) arrow$ends <- 2L
    con_gp <- subset_gp(con_gp, labels_drawn)
    connect <- polylineGrob(connect$x, connect$y,
      id = connect$id,
      default.units = 'mm', gp = con_gp, arrow = arrow
    )
  }
  c(labels, list(connect))
}


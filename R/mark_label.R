# Label placement (mascarade boundary-seed placer). Runs at draw time in the panel's
# millimetre space, where `rects` are the measured label box sizes (w, h in mm), `polygons`
# are the cluster rings in mm, and `bounds` is the panel size in mm. Returns, per label, its
# placed centre c(x, y) in mm (NULL for labels not drawn). Poles and the box-fit R-tree are
# recomputed here each draw (cheap: ~20 ms for ~40 clusters); the expensive mask is not.
# `anchors` is used only as the degenerate-input fallback; non-labelled polygon parts are
# treated as the placed set (points to avoid / ghosts are not handled yet).
#' @importFrom polylabelr poi
#' @importFrom stats median
my_place_labels <- function(rects, polygons, bounds, anchors) {
  res <- vector('list', length(rects))
  active <- which(vapply(rects, function(r) !all(r == 0), logical(1)))
  if (length(active) == 0) return(res)

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
  poimat <- t(vapply(polys_a, function(p) {
    pl <- tryCatch(polylabelr::poi(p$x, p$y), error = function(e) NULL)
    if (is.null(pl) || !is.finite(pl$x) || !is.finite(pl$y)) c(mean(p$x), mean(p$y)) else c(pl$x, pl$y)
  }, numeric(2)))

  # anchor fallback if anything is still degenerate
  anchorFallback <- function() {
    for (j in seq_along(active)) res[[active[j]]] <<- as.numeric(anchors[[active[j]]])
    res
  }
  if (any(!is.finite(poimat)) || any(!is.finite(bounds)) || bounds[1] <= 0 || bounds[2] <= 0)
    return(anchorFallback())
  if (length(active) == 1) { res[[active]] <- poimat[1, ]; return(res) }

  hw <- vapply(active, function(i) rects[[i]][1] / 2, 0)
  hh <- vapply(active, function(i) rects[[i]][2] / 2, 0)
  char_h <- stats::median(2 * hh)
  px <- lapply(polys_a, `[[`, 'x'); py <- lapply(polys_a, `[[`, 'y')
  geom <- list(poi = poimat, rtree = buildBoxFit(px, py), polysx = px, polysy = py)

  lay <- tryCatch(placeLabels(geom, c(0, bounds[1]), c(0, bounds[2]), hw, hh, char_h),
                  error = function(e) NULL)
  if (is.null(lay)) return(anchorFallback())
  lay <- lay[order(lay$label)]
  for (j in seq_along(active)) res[[active[j]]] <- c(lay$cx[j], lay$cy[j])
  res
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
  labelpos <- my_place_labels(dims, p_big, area, anchors)
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
  connect <- rlang::inject(rbind(!!!Map(function(pol, pos, dim) {
    if (is.null(pos)) return(NULL)
    dim <- dim / anchor_mod
    pos <- cbind(
      c(pos[1] - dim[1], pos[1] + dim[1], pos[1] + dim[1], pos[1] - dim[1]),
      c(pos[2] - dim[2], pos[2] - dim[2], pos[2] + dim[2], pos[2] + dim[2])
    )
    pos <- points_to_path(pos, list(cbind(pol$x, pol$y)), TRUE)
    pos$projection[which.min(pos$distance), ]
  }, pol = polygons, pos = labelpos, dim = dims)))
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


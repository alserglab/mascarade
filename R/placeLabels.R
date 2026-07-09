# Anchor-leader label placement driver.
#
# Pipeline: boundary seed -> radial candidates (+ seed fallback, visible-length filtered)
# -> one-move sweep -> all-pairs two-move (length B&B, lexicographic) -> continuous force
# polish (off the candidate grid). The leader START on the label box is decoupled from the
# box centre (see .anchorPoint); the C++ kernels in src/ score conflicts on the actual
# anchor->pole segment (forcePolish mirrors .anchorPoint). This file is the R orchestration.
#
# `placeLabels()` is pure given the per-label box half-sizes (hw, hh) and line height
# (char_h); the draw-stage hook supplies those from text metrics in the panel's mm space.
# `geom` is the box-fit structure the caller builds from the cluster polygons:
#   list(poi = K x 2 pole matrix, rtree = XPtr<BoxFit>, polysx, polysy = per-cluster rings)
#
# `con_type` selects the leader style: "cl" corners+ledge (one sign-quadrant corner),
# "cm" corners+midpoints (8-point rule, no ledge), "none" (placed like cm, no connector).
# The keep-out padding between labels and clusters is baked into `geom` itself: the caller
# passes polygons already dilated by `label.buffer`, so box-fit forbids boxes within that
# distance of the true cluster (see my_make_label).

# pool-adjacent-violators (non-increasing L2 fit)
.pavaDec <- function(y) {
  vals <- numeric(0); wts <- numeric(0)
  for (i in seq_along(y)) {
    cv <- y[i]; cw <- 1
    while (length(vals) > 0 && vals[length(vals)] < cv) {
      k <- length(vals); cv <- (vals[k] * wts[k] + cv * cw) / (wts[k] + cw)
      cw <- wts[k] + cw; vals <- vals[-k]; wts <- wts[-k]
    }
    vals <- c(vals, cv); wts <- c(wts, cw)
  }
  rep(vals, wts)
}

# isotonic 1D placement: centres closest to targets t with min separation (h_i+h_j)/2 + gap
.place1d <- function(t, h, gap) {
  n <- length(t); if (n == 1) return(t)
  d <- (h[-n] + h[-1]) / 2 + gap; S <- c(0, cumsum(d))
  e <- .pavaDec(t + S); e - S
}

# The point on the box border where the leader starts, aimed at the pole (tx,ty). xin/yin
# mark whether the pole coordinate falls between the box edges (ggforce get_end_points rule).
#   "cl": always the sign(dx),sign(dy) quadrant corner.
#   else: corner only when the pole is fully to the side in BOTH axes; if the pole's x is
#         between the vertical edges -> top/bottom edge centre; if its y is between the
#         horizontal edges -> left/right edge centre.
# Returns list(ex, ey, corner) (vectorised).
.anchorPoint <- function(cx, cy, hw, hh, tx, ty, con_type) {
  dx <- tx - cx; dy <- ty - cy
  sX <- ifelse(dx >= 0, 1, -1); sY <- ifelse(dy >= 0, 1, -1)
  if (con_type == "cl")
    return(list(ex = cx + sX * hw, ey = cy + sY * hh, corner = rep(TRUE, length(cx))))
  xin <- abs(dx) < hw; yin <- abs(dy) < hh
  list(ex = ifelse(xin & !yin, cx, cx + sX * hw),
       ey = ifelse(!xin & yin, cy, cy + sY * hh),
       corner = !xin & !yin)
}

# First intersection of each leader (anchor_i -> pole) with a single polygon ring, as a
# fraction t in (0,1] of the leader; visible length = t * |anchor - pole|. Vectorised over
# anchors; the pole (tx,ty) is inside the ring. Returns list(vis, hx, hy).
.firstHit <- function(ax, ay, tx, ty, px, py) {
  n <- length(px); dx <- tx - ax; dy <- ty - ay
  best <- rep(Inf, length(ax))
  for (e in seq_len(n)) {
    x1 <- px[e]; y1 <- py[e]; x2 <- px[e %% n + 1L]; y2 <- py[e %% n + 1L]
    ex <- x2 - x1; ey <- y2 - y1
    den <- dx * ey - dy * ex
    ok <- abs(den) > 1e-12
    t <- ((x1 - ax) * ey - (y1 - ay) * ex) / den   # along leader
    u <- ((x1 - ax) * dy - (y1 - ay) * dx) / den   # along edge
    hit <- ok & t > 1e-6 & t <= 1 + 1e-9 & u >= -1e-9 & u <= 1 + 1e-9
    best <- ifelse(hit & t < best, t, best)
  }
  best[!is.finite(best)] <- 1                      # no crossing (degenerate): whole segment
  d <- sqrt(dx^2 + dy^2)
  list(vis = best * d, hx = ax + dx * best, hy = ay + dy * best)
}

# derive box + leader geometry columns for a (label, cx, cy) table. `len` is the ranking
# length (centre->pole distance); ex,ey,corner come from the anchor rule.
.geoCols <- function(dt, hw, hh, poi, pad, con_type) {
  L <- dt$label; h_w <- hw[L]; h_h <- hh[L]; tx <- poi[L, 1]; ty <- poi[L, 2]
  a <- .anchorPoint(dt$cx, dt$cy, h_w, h_h, tx, ty, con_type)
  data.table::data.table(
    label = L, cx = dt$cx, cy = dt$cy, hw = h_w, hh = h_h, tx = tx, ty = ty,
    cxmin = dt$cx - h_w - pad, cxmax = dt$cx + h_w + pad,
    cymin = dt$cy - h_h - pad, cymax = dt$cy + h_h + pad,
    ex = a$ex, ey = a$ey, corner = a$corner, len = sqrt((dt$cx - tx)^2 + (dt$cy - ty)^2))
}

# Conflict-free boundary seed (guaranteed for equal box heights). Two columns hang off the
# cluster cloud on the vertical lines X = min/max polygon x (the polygons are already dilated
# by label.buffer, so this sits a keep-out margin outside the true clusters); each box's
# padded near edge sits on its line and the leader starts on the line (bottom corner for
# "cl", edge centre otherwise), so it never clips a box. Labels are matched to stacked slots
# by leader length via the Hungarian assignment (min total length => crossing-free).
.reorderBase <- function(poi, hw, hh, polyxlim, gap, pad, con_type) {
  K <- nrow(poi); fh <- 2 * hh
  ox <- order(poi[, 1]); cs <- cumsum(fh[ox]); k <- which(cs >= sum(fh) / 2)[1]
  if (is.na(k)) k <- K; k <- max(1L, min(k, K - 1L))
  Lset <- ox[seq_len(k)]; Rset <- ox[(k + 1L):K]
  XL <- polyxlim[1]; XR <- polyxlim[2]
  col <- function(set, Xline, side) {              # side: -1 left column, +1 right column
    m <- length(set)
    if (m == 1) return(data.table::data.table(label = set, cy = poi[set, 2], Xline = Xline, side = side))
    lab <- set[order(-poi[set, 2])]
    for (it in seq_len(10L)) {
      sy <- .place1d(poi[lab, 2], fh[lab], gap)
      Cm <- sqrt(outer(poi[set, 1], sy, function(a, b) (Xline - a)^2) +   # leader length^2 from
                 outer(poi[set, 2], sy, function(a, b) (a - b)^2))        # (Xline, slot y) to pole
      asg <- hungarian(Cm) + 1L
      newlab <- integer(m); newlab[asg] <- set
      if (identical(newlab, lab)) break
      lab <- newlab
    }
    data.table::data.table(label = lab, cy = .place1d(poi[lab, 2], fh[lab], gap),
                           Xline = Xline, side = side)
  }
  s <- rbind(col(Lset, XL, -1), col(Rset, XR, +1))
  # box centre so the padded near edge sits on the line; leader starts at the true (unpadded)
  # near edge so a "cl" ledge meets it. It stays a padded-margin width inside the line, still
  # clear of every other column box (they are separated in y by gap > pad).
  s[, `:=`(cx = Xline + side * (hw[label] + pad), ex = Xline + side * pad,
           ey = if (con_type == "cl") cy - hh[label] else cy, corner = con_type == "cl")]
  s[]
}

# Place labels for one view. Returns a data.table (one row per cluster) with cx, cy, the box
# columns, the leader anchor (ex, ey), its `corner` flag, and the visible leader end (bx, by
# = first mask-boundary hit). Conflict-free by construction given a feasible pool.
placeLabels <- function(geom, xlim, ylim, hw, hh, char_h, con_type = "cl",
                        MU = 55, iters = 120L) {
  poi <- geom$poi; K <- nrow(poi)
  pad <- 0.05 * char_h                                             # hard box clearance
  gap <- 0.25 * char_h                                             # seed column spacing
  if (K == 1) {
    r <- .geoCols(data.table::data.table(label = 1L, cx = poi[1, 1], cy = poi[1, 2]),
                  hw, hh, poi, pad, con_type)
    return(r[, `:=`(bx = poi[1, 1], by = poi[1, 2])][])            # single label sits on its pole
  }

  ndir <- 48L; radStep <- 0.3 * char_h; radStart <- 0.2 * char_h
  radReach <- 16 * char_h; radFill <- 1.2 * char_h; dedup <- 0.3 * char_h
  polyxlim <- geom$pad_xrange                                       # dilated extent (keep-out)
  if (is.null(polyxlim)) polyxlim <- range(unlist(geom$polysx))     # undilated geom (tests)

  cand <- data.table::as.data.table(radialCandidates(
    geom$rtree, poi, hw, hh, pad, xlim[1], xlim[2], ylim[1], ylim[2],
    ndir, radStep, radStart, radReach, radFill, dedup))
  pool <- .geoCols(cand, hw, hh, poi, pad, con_type)
  pool$isb <- FALSE

  # boundary seed = guaranteed-clean fallback slot / label
  sd <- .reorderBase(poi, hw, hh, polyxlim, gap, pad, con_type)
  seed <- data.table::data.table(
    label = sd$label, cx = sd$cx, cy = sd$cy, hw = hw[sd$label], hh = hh[sd$label],
    tx = poi[sd$label, 1], ty = poi[sd$label, 2],
    cxmin = sd$cx - hw[sd$label] - pad, cxmax = sd$cx + hw[sd$label] + pad,
    cymin = sd$cy - hh[sd$label] - pad, cymax = sd$cy + hh[sd$label] + pad,
    ex = sd$ex, ey = sd$ey, corner = sd$corner,
    len = sqrt((sd$cx - poi[sd$label, 1])^2 + (sd$cy - poi[sd$label, 2])^2), isb = TRUE)
  pool <- rbind(pool, seed)

  pool <- pool[order(pool$label, pool$len)]
  pool$idx <- seq_len(nrow(pool)) - 1L                             # 0-based for the C++ kernels
  rows <- lapply(seq_len(K), function(i) pool$idx[pool$label == i])
  init <- vapply(seq_len(K), function(i) {
    sel <- pool$label == i; pool$idx[sel][pool$isb[sel]][1] }, 0L)

  # effective length = leader length + length inside a foreign cluster (routes leaders around)
  elen <- pool$len + foreignLength(pool$ex, pool$ey, pool$tx, pool$ty,
                                   as.integer(pool$label), geom$polysx, geom$polysy)
  geomArgs <- list(cxmin = pool$cxmin, cxmax = pool$cxmax, cymin = pool$cymin, cymax = pool$cymax,
                   ex = pool$ex, ey = pool$ey, tx = pool$tx, ty = pool$ty, len = elen, rows = rows)
  rs <- do.call(oneMoveSweep, c(geomArgs, list(init = as.integer(init), maxpass = 100L)))
  tw <- do.call(twoMoveBnB, c(geomArgs, list(init = as.integer(rs), maxpass = 50L, sq = TRUE)))

  # continuous force-directed polish off the candidate grid (anchor-aware; preserves the
  # discrete solution's conflict-freeness). con_type 0 = "cl" corner, else "cm"/"none".
  two <- pool[tw + 1L][order(label)]
  r <- forcePolish(geom$rtree, two$cx, two$cy, hw, hh, poi[, 1], poi[, 2], pad,
                   xlim[1], xlim[2], ylim[1], ylim[2], as.integer(iters), 0.4 * char_h, MU,
                   0.6 * char_h, 0.03 * char_h, if (con_type == "cl") 0L else 1L)
  res <- .geoCols(data.table::data.table(label = seq_len(K), cx = r$cx, cy = r$cy),
                  hw, hh, poi, pad, con_type)[order(label)]
  # visible leader end (first mask-boundary hit) for drawing
  hit <- do.call(rbind, lapply(seq_len(nrow(res)), function(i) {
    L <- res$label[i]
    h <- .firstHit(res$ex[i], res$ey[i], poi[L, 1], poi[L, 2], geom$polysx[[L]], geom$polysy[[L]])
    c(h$hx, h$hy) }))
  res[, `:=`(bx = hit[, 1], by = hit[, 2])][]
}

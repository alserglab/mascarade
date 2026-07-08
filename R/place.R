# Boundary-seed label placement driver.
#
# Pipeline: min-cost seed -> radial candidates (+ boundary-slot fallback) -> one-move
# sweep -> all-pairs two-move (length B&B, lexicographic) -> force-directed polish.
# The C++ kernels live in src/; this file is the R orchestration ("prep" glue: candidate
# pool, per-label row lists, index maps, and the derived box/leader geometry columns).
#
# `placeLabels()` is pure given the per-label box half-sizes (hw, hh) and line height
# (char_h) in data units; `labelBoxes()`/`fancyMask()` supply those from text metrics at
# draw time. `geom` is the compute-once structure built by `fancyMask()`:
#   list(poi = K x 2 pole matrix, rtree = XPtr<BoxFit>, polysx, polysy = per-cluster rings)

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

# box-edge point of the leader (centre -> pole), vectorized
.leaderEdges <- function(cx, cy, hw, hh, tx, ty) {
  dx <- tx - cx; dy <- ty - cy; adx <- abs(dx); ady <- abs(dy)
  tt <- pmin(ifelse(adx > 1e-9, hw / adx, Inf), ifelse(ady > 1e-9, hh / ady, Inf))
  tt[!is.finite(tt)] <- 0
  list(ex = cx + dx * tt, ey = cy + dy * tt)
}

# derive box + leader geometry columns for a (label, cx, cy) table
.geoCols <- function(dt, hw, hh, poi, pad) {
  L <- dt$label
  h_w <- hw[L]; h_h <- hh[L]; tx <- poi[L, 1]; ty <- poi[L, 2]
  e <- .leaderEdges(dt$cx, dt$cy, h_w, h_h, tx, ty)
  data.table::data.table(
    label = L, cx = dt$cx, cy = dt$cy, hw = h_w, hh = h_h, tx = tx, ty = ty,
    cxmin = dt$cx - h_w - pad, cxmax = dt$cx + h_w + pad,
    cymin = dt$cy - h_h - pad, cymax = dt$cy + h_h + pad,
    ex = e$ex, ey = e$ey, len = sqrt((dt$cx - tx)^2 + (dt$cy - ty)^2))
}

# min-cost boundary seed: height-balanced x-cut split, per-column Hungarian matching
# (crossing-free); returns a (label, cx, cy) table. hw/hh are half-sizes.
.reorderBase <- function(poi, hw, hh, xlim, gap) {
  K <- nrow(poi); fh <- 2 * hh
  ox <- order(poi[, 1]); cs <- cumsum(fh[ox]); k <- which(cs >= sum(fh) / 2)[1]
  if (is.na(k)) k <- K; k <- max(1L, min(k, K - 1L))
  Lset <- ox[seq_len(k)]; Rset <- ox[(k + 1L):K]
  XcL <- xlim[1] + max(hw[Lset]); XcR <- xlim[2] - max(hw[Rset])
  col <- function(set, Xc) {
    m <- length(set)
    if (m == 1) return(data.table::data.table(label = set, cx = Xc, cy = poi[set, 2]))
    lab <- set[order(-poi[set, 2])]
    for (it in seq_len(10L)) {
      sy <- .place1d(poi[lab, 2], fh[lab], gap)
      Cm <- sqrt(outer(poi[set, 1], sy, function(a, b) (Xc - a)^2) +
                 outer(poi[set, 2], sy, function(a, b) (a - b)^2))
      asg <- hungarian(Cm) + 1L
      newlab <- integer(m); newlab[asg] <- set
      if (identical(newlab, lab)) break
      lab <- newlab
    }
    data.table::data.table(label = lab, cx = Xc, cy = .place1d(poi[lab, 2], fh[lab], gap))
  }
  rbind(col(Lset, XcL), col(Rset, XcR))
}

# Place labels for one view. Returns a data.table (one row per cluster) with cx, cy and
# the derived box/leader columns. Conflict-free by construction given a feasible pool.
placeLabels <- function(geom, xlim, ylim, hw, hh, char_h,
                        MU = 55, iters = 120L, hard_ll = TRUE) {
  poi <- geom$poi; K <- nrow(poi)
  pad <- 0.05 * char_h; gap <- 0.25 * char_h
  if (K == 1) return(.geoCols(data.table::data.table(label = 1L, cx = poi[1, 1], cy = poi[1, 2]),
                              hw, hh, poi, pad))                    # single label sits on its pole

  seed <- .reorderBase(poi, hw, hh, xlim, gap)                      # guaranteed-clean fallback slot / label
  cand <- data.table::as.data.table(radialCandidates(
    geom$rtree, poi, hw, hh, pad, xlim[1], xlim[2], ylim[1], ylim[2],
    48L, 0.3 * char_h, 0.2 * char_h, 16 * char_h, 1.2 * char_h, 0.3 * char_h))

  cand$isb <- FALSE
  sb <- data.table::data.table(label = seed$label, cx = seed$cx, cy = seed$cy, isb = TRUE)
  allc <- rbind(cand, sb)                                          # boundary slot = feasibility fallback
  pool <- .geoCols(allc, hw, hh, poi, pad)
  pool$isb <- allc$isb
  ord <- order(pool$label, pool$len); pool <- pool[ord]
  pool$idx <- seq_len(nrow(pool)) - 1L
  rows <- lapply(seq_len(K), function(i) pool$idx[pool$label == i])
  init <- vapply(seq_len(K), function(i) {
    sel <- pool$label == i; pool$idx[sel][pool$isb[sel]][1] }, 0L)

  elen <- pool$len + foreignLength(pool$ex, pool$ey, pool$tx, pool$ty,
                                   as.integer(pool$label), geom$polysx, geom$polysy)
  a <- list(pool$cxmin, pool$cxmax, pool$cymin, pool$cymax,
            pool$ex, pool$ey, pool$tx, pool$ty, elen, rows)
  rs <- do.call(oneMoveSweep, c(a, list(as.integer(init), 100L)))
  tw <- do.call(twoMoveBnB, c(a, list(as.integer(rs), 50L, TRUE)))

  two <- .geoCols(pool[tw + 1L], hw, hh, poi, pad)
  two <- two[order(two$label)]
  r <- forcePolish(geom$rtree, two$cx, two$cy, hw, hh, poi[, 1], poi[, 2], pad,
                   xlim[1], xlim[2], ylim[1], ylim[2],
                   as.integer(iters), 0.4 * char_h, MU, 0.6 * char_h, 0.03 * char_h,
                   hard_ll, TRUE)
  .geoCols(data.table::data.table(label = seq_len(K), cx = r$cx, cy = r$cy), hw, hh, poi, pad)
}

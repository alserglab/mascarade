# Anchor-leader label placement driver.
#
# Pipeline (see `placeLabels()`): candidatePool() -> oneMoveSweep() -> twoMoveSweep() ->
# polishLayout() -> withLeaderEnds(). Each stage takes the placement structure and returns it,
# so the driver reads as a straight chain. The leader START on the label box is decoupled from
# the box centre (see `.anchorPoint`); the C++ kernels in src/ score conflicts on the actual
# anchor->pole segment (forcePolish mirrors `.anchorPoint`). This file is the R orchestration.
#
# Two structures flow through the pipeline:
#   * Layout (see `.layout`) -- a one-row-per-label data.table carrying the box centre, padded
#     box extents, leader anchor and ranking length: the columns the C++ kernels and the test
#     scorer read. The final placement is a Layout.
#   * Pool (see `candidatePool`) -- the discrete search state: a set of candidate Layout rows
#     (many per label), the per-label candidate index lists, and the currently selected index
#     per label. oneMoveSweep()/twoMoveSweep() take a Pool and return a Pool with the selection
#     improved; polishLayout() collapses the selection to a single polished Layout.
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

#' Leader start anchor on a label box
#'
#' The point on the box border where the leader begins, aimed at the pole (`tx`, `ty`). For
#' `con_type == "cl"` this is always the sign(dx),sign(dy) quadrant corner. Otherwise (the
#' ggforce `get_end_points` rule) it is a corner only when the pole is fully to the side in
#' both axes; if the pole's x lies between the vertical edges it is the top/bottom edge
#' centre, and if its y lies between the horizontal edges it is the left/right edge centre.
#'
#' @param cx,cy Numeric box-centre coordinates (vectorised).
#' @param hw,hh Numeric box half-width / half-height.
#' @param tx,ty Numeric pole (leader target) coordinates.
#' @param con_type Leader style: `"cl"`, `"cm"`, or `"none"`.
#' @return A list with numeric `ex`, `ey` (the anchor) and logical `corner` (is it a corner?).
#' @keywords internal
#' @noRd
.anchorPoint <- function(cx, cy, hw, hh, tx, ty, con_type) {
  dx <- tx - cx
  dy <- ty - cy
  signX <- ifelse(dx >= 0, 1, -1)
  signY <- ifelse(dy >= 0, 1, -1)
  if (con_type == "cl") {
    return(list(ex = cx + signX * hw, ey = cy + signY * hh,
                corner = rep(TRUE, length(cx))))
  }
  xInside <- abs(dx) < hw
  yInside <- abs(dy) < hh
  list(ex = ifelse(xInside & !yInside, cx, cx + signX * hw),
       ey = ifelse(!xInside & yInside, cy, cy + signY * hh),
       corner = !xInside & !yInside)
}

#' Build a Layout from explicit box centres and leader anchors
#'
#' The Layout structure shared across the placement pipeline: one row per label with the box
#' centre, padded box extents, pole, leader anchor and ranking length (centre-to-pole
#' distance). This is the low-level constructor -- the caller supplies the leader anchor
#' directly (the boundary seed needs anchors that do not follow the `.anchorPoint()` rule).
#' Use `.layoutFromCentres()` when the anchor should follow the leader style.
#'
#' @param label Integer 1-indexed cluster of each row.
#' @param cx,cy Numeric box-centre coordinates.
#' @param ex,ey Numeric leader-start (anchor) coordinates.
#' @param corner Logical: is the anchor a box corner?
#' @param hw,hh Numeric per-label box half-sizes (indexed by `label`).
#' @param poi K x 2 matrix of cluster poles.
#' @param pad Numeric hard box clearance added around each box.
#' @return A data.table with one row per input, carrying `label`, `cx`, `cy`, `hw`, `hh`,
#'   `tx`, `ty`, the padded extents `cxmin`/`cxmax`/`cymin`/`cymax`, `ex`, `ey`, `corner`
#'   and `len`.
#' @keywords internal
#' @noRd
.layout <- function(label, cx, cy, ex, ey, corner, hw, hh, poi, pad) {
  halfW <- hw[label]
  halfH <- hh[label]
  tx <- poi[label, 1]
  ty <- poi[label, 2]
  data.table::data.table(
    label = label, cx = cx, cy = cy, hw = halfW, hh = halfH, tx = tx, ty = ty,
    cxmin = cx - halfW - pad, cxmax = cx + halfW + pad,
    cymin = cy - halfH - pad, cymax = cy + halfH + pad,
    ex = ex, ey = ey, corner = corner,
    len = sqrt((cx - tx)^2 + (cy - ty)^2))
}

#' Build a Layout from box centres, deriving leader anchors from the leader style
#'
#' Convenience wrapper over `.layout()`: derives each leader anchor from the box centre and
#' pole via `.anchorPoint()`, then assembles the Layout.
#'
#' @param dt A data.table/data.frame with columns `label`, `cx`, `cy`.
#' @param hw,hh Numeric per-label box half-sizes (indexed by `label`).
#' @param poi K x 2 matrix of cluster poles.
#' @param pad Numeric hard box clearance added around each box.
#' @param con_type Leader style passed to `.anchorPoint()`.
#' @return A Layout data.table (see `.layout()`).
#' @keywords internal
#' @noRd
.layoutFromCentres <- function(dt, hw, hh, poi, pad, con_type) {
  anchor <- .anchorPoint(dt$cx, dt$cy, hw[dt$label], hh[dt$label],
                         poi[dt$label, 1], poi[dt$label, 2], con_type)
  .layout(dt$label, dt$cx, dt$cy, anchor$ex, anchor$ey, anchor$corner,
          hw, hh, poi, pad)
}

#' One boundary-seed column of stacked label slots
#'
#' Places the labels assigned to one side (`side` -1 left, +1 right) as a vertical stack on the
#' line `x = Xline`. Each side chooses its slot height once:
#'
#' * If the labels all fit the viewport in slots as tall as the tallest box (`tallest + gap`), a
#'   single Hungarian assigns labels to those slots -- one box per slot, so any assignment is
#'   overlap-free and no packing is needed (single- and multi-line labels alike). The slots cover
#'   the pole span (a bounded window clamped into the viewport) so labels sit near their poles.
#' * Otherwise the column is too crowded for full-height slots: the Hungarian only fixes the
#'   top-to-bottom order and `packLen()` packs the mixed-height boxes tightly on a one-line grid
#'   (`one line height + gap`), extended past the viewport when they still cannot fit.
#'
#' @param set Integer labels assigned to this column.
#' @param Xline Numeric x-coordinate of the column line.
#' @param side Integer -1 (left column) or +1 (right column).
#' @param poi K x 2 matrix of cluster poles.
#' @param boxH Numeric per-label full box heights.
#' @param gap Numeric seed column spacing.
#' @param char_h Numeric line height (the fine packing grid uses `char_h + gap`).
#' @param ylim Numeric length-2 viewport y-bounds (slots span at least this).
#' @return A data.table with `label`, `cy`, `Xline`, `side`.
#' @keywords internal
#' @noRd
.seedColumn <- function(set, Xline, side, poi, boxH, gap, char_h, ylim) {
  m <- length(set)
  poleY <- poi[set, 2]
  if (m == 1) {
    return(data.table::data.table(label = set, cy = poleY, Xline = Xline, side = side))
  }
  dx2 <- (Xline - poi[set, 1])^2               # squared horizontal pole-to-column distance
  viewH <- ylim[2] - ylim[1]
  coarseSlot <- max(boxH[set]) + gap           # uniform slot as tall as the tallest box
  capacity <- floor(viewH / coarseSlot)        # such slots the viewport holds

  if (m <= capacity) {
    # room for every label to get its own tallest-box slot: the Hungarian assignment is itself
    # overlap-free, so no packing DP is needed. Slots cover the pole span (a bounded window,
    # clamped into the viewport) so labels sit near their poles.
    spanSlots <- ceiling((max(poleY) - min(poleY)) / coarseSlot)
    nSlot <- max(m, min(capacity, spanSlots + 2L * m))
    lo <- mean(range(poleY)) - nSlot * coarseSlot / 2
    lo <- max(ylim[1], min(lo, ylim[2] - nSlot * coarseSlot))
    slotY <- lo + (seq_len(nSlot) - 0.5) * coarseSlot
    cost <- matrix(0, nSlot, nSlot)            # square: m real rows + dummy zero-cost rows
    for (i in seq_len(m)) {
      cost[i, ] <- sqrt(dx2[i] + (poleY[i] - slotY)^2)
    }
    assignment <- hungarian(cost) + 1L
    return(data.table::data.table(label = set, cy = slotY[assignment[seq_len(m)]],
                                  Xline = Xline, side = side))
  }

  # too crowded for full-height slots: the Hungarian fixes the top-to-bottom order and packLen()
  # packs the mixed-height boxes tightly on a one-line grid (extended past the viewport if needed).
  orderY <- mean(range(poleY)) + (seq_len(m) - (m + 1) / 2) * coarseSlot
  cost <- sqrt(outer(dx2, rep(1, m)) +         # leader length from label i to order-slot j
               outer(poleY, orderY, function(a, b) (a - b)^2))
  assignment <- hungarian(cost) + 1L
  ordered <- set[order(-orderY[assignment])]   # top-to-bottom label order
  data.table::data.table(
    label = ordered,
    cy = packLen(Xline - poi[ordered, 1], poi[ordered, 2], boxH[ordered],
                 gap, char_h + gap, ylim[1], ylim[2]),
    Xline = Xline, side = side)
}

#' Conflict-free boundary seed
#'
#' Guaranteed-clean fallback placement: two columns hang off the cluster cloud on the vertical
#' lines `x = min/max` polygon x. Since the polygons are already dilated by `label.buffer`, those
#' lines sit a keep-out margin outside the true clusters; each box's padded near edge sits on its
#' line and the leader starts on the line (bottom corner for `"cl"`, edge centre otherwise), so it
#' never clips a box. Labels are split left/right by pole x (balancing total height) and each
#' column is stacked by `.seedColumn()`.
#'
#' @param poi K x 2 matrix of cluster poles.
#' @param hw,hh Numeric per-label box half-sizes.
#' @param polyxlim Numeric length-2 x-range of the (dilated) cluster polygons.
#' @param ylim Numeric length-2 viewport y-bounds (slots span at least this).
#' @param char_h Numeric line height (the fine packing grid uses `char_h + gap`).
#' @param gap Numeric seed column spacing.
#' @param pad Numeric hard box clearance.
#' @param con_type Leader style (`"cl"` anchors at the bottom corner, else the edge centre).
#' @return A data.table with `label`, `cx`, `cy`, `ex`, `ey`, `corner`.
#' @keywords internal
#' @noRd
.reorderBase <- function(poi, hw, hh, polyxlim, ylim, char_h, gap, pad, con_type) {
  K <- nrow(poi)
  boxH <- 2 * hh                               # full box height per label

  # split labels into a left / right column by pole x, balancing total box height
  byX <- order(poi[, 1])
  cumH <- cumsum(boxH[byX])
  split <- which(cumH >= sum(boxH) / 2)[1]
  if (is.na(split)) {
    split <- K
  }
  split <- max(1L, min(split, K - 1L))
  leftSet <- byX[seq_len(split)]
  rightSet <- byX[(split + 1L):K]

  slots <- rbind(
    .seedColumn(leftSet, polyxlim[1], -1, poi, boxH, gap, char_h, ylim),
    .seedColumn(rightSet, polyxlim[2], +1, poi, boxH, gap, char_h, ylim))

  # box centre so the padded near edge sits on the column line; the leader starts on the true
  # (unpadded) near edge so a "cl" ledge meets it. It stays a padded-margin width inside the
  # line, still clear of every other column box (they are separated in y by gap > pad).
  slots[, `:=`(
    cx = Xline + side * (hw[label] + pad),
    ex = Xline + side * pad,
    ey = if (con_type == "cl") cy - hh[label] else cy,
    corner = con_type == "cl")]
  slots[]
}

#' Candidate pool for one view: radial free-space candidates + boundary-seed fallback
#'
#' Builds the Pool the discrete refinement searches over. It gathers many radial free-space
#' candidate Layout rows per label (from `radialCandidates()`) plus one guaranteed-clean
#' boundary-seed row per label (from `.reorderBase()`), ranks every candidate by its effective
#' length (leader length + arc inside foreign clusters + viewport overflow, from
#' `effectiveLength()`), orders the rows by `(label, len)` and gives each a 0-based `idx` for the
#' C++ kernels.
#'
#' @param geom Box-fit structure (see `placeLabels()`).
#' @param xlim,ylim Numeric length-2 viewport bounds.
#' @param hw,hh Numeric per-label box half-sizes.
#' @param char_h Numeric line height.
#' @param pad Numeric hard box clearance.
#' @param gap Numeric seed column spacing.
#' @param con_type Leader style.
#' @return A Pool: a list with `cand` (the candidate Layout, plus `eff` effective length, `isb`
#'   seed flag and 0-based `idx`), `rows` (`rows[[i]]` = label i's candidate indices) and `sel`
#'   (the currently selected `idx` per label, initialised to each label's boundary seed -- the
#'   conflict-free starting pick).
#' @keywords internal
#' @noRd
candidatePool <- function(geom, xlim, ylim, hw, hh, char_h, pad, gap, con_type) {
  poi <- geom$poi
  K <- nrow(poi)

  # radial free-space candidates around every pole (many per label)
  ndir <- 48L
  radStep <- 0.3 * char_h
  radStart <- 0.2 * char_h
  radReach <- 16 * char_h
  radFill <- 1.2 * char_h
  dedup <- 0.3 * char_h
  cand <- data.table::as.data.table(radialCandidates(
    geom$rtree, poi, hw, hh, pad, xlim[1], xlim[2], ylim[1], ylim[2],
    ndir, radStep, radStart, radReach, radFill, dedup))
  candidates <- .layoutFromCentres(cand, hw, hh, poi, pad, con_type)
  candidates$isb <- FALSE

  # guaranteed-clean boundary seed: one fallback slot per label off the cluster cloud
  polyxlim <- geom$pad_xrange                  # dilated extent (keep-out)
  if (is.null(polyxlim)) {
    polyxlim <- range(unlist(geom$polysx))     # undilated geom (tests)
  }
  seedSlots <- .reorderBase(poi, hw, hh, polyxlim, ylim, char_h, gap, pad, con_type)
  seed <- .layout(seedSlots$label, seedSlots$cx, seedSlots$cy,
                  seedSlots$ex, seedSlots$ey, seedSlots$corner, hw, hh, poi, pad)
  seed$isb <- TRUE

  cand <- rbind(candidates, seed)
  cand <- cand[order(cand$label, cand$len)]
  cand$idx <- seq_len(nrow(cand)) - 1L         # 0-based row id for the C++ kernels

  # effective length = leader length + arc inside foreign clusters (routes leaders around them)
  # + how far the box overflows the viewport (steers labels in-bounds); the kernels rank by it
  cand$eff <- effectiveLength(
    cand$len, cand$ex, cand$ey, cand$tx, cand$ty, as.integer(cand$label),
    geom$polysx, geom$polysy, cand$cxmin, cand$cxmax, cand$cymin, cand$cymax,
    xlim[1], xlim[2], ylim[1], ylim[2])

  rows <- lapply(seq_len(K), function(i) cand$idx[cand$label == i])
  sel <- vapply(seq_len(K), function(i) {
    cand$idx[cand$isb & cand$label == i][1]    # each label has exactly one seed row
  }, 0L)
  list(cand = cand, rows = rows, sel = sel)
}

#' One-move conflict sweep over a Pool
#'
#' Wraps `oneMoveSweepKernel()`: moves each label to its lexicographically best candidate versus
#' the current others (including staying put), ranking by effective length. Takes a Pool and
#' returns it with the selection improved.
#'
#' @param pool A Pool (see `candidatePool()`).
#' @param maxpass Integer cap on the number of sweeps.
#' @return The Pool with `sel` updated.
#' @keywords internal
#' @noRd
oneMoveSweep <- function(pool, maxpass = 100L) {
  cand <- pool$cand
  pool$sel <- oneMoveSweepKernel(
    cand$cxmin, cand$cxmax, cand$cymin, cand$cymax, cand$ex, cand$ey, cand$tx, cand$ty,
    cand$eff, pool$rows, as.integer(pool$sel), maxpass)
  pool
}

#' Two-move conflict / length refinement over a Pool
#'
#' Wraps `twoMoveSweepKernel()`: per label, applies the best two-move (length branch-and-bound
#' when the label is conflict-free, all-pairs search when it is conflicted). Takes a Pool and
#' returns it with the selection improved.
#'
#' @param pool A Pool (see `candidatePool()`).
#' @param maxpass Integer cap on the number of passes.
#' @param sq Logical; if `TRUE` the length objective uses the squared length.
#' @return The Pool with `sel` updated.
#' @keywords internal
#' @noRd
twoMoveSweep <- function(pool, maxpass = 50L, sq = TRUE) {
  cand <- pool$cand
  pool$sel <- twoMoveSweepKernel(
    cand$cxmin, cand$cxmax, cand$cymin, cand$cymax, cand$ex, cand$ey, cand$tx, cand$ty,
    cand$eff, pool$rows, as.integer(pool$sel), maxpass, sq)
  pool
}

#' Continuous force-directed polish of a Pool's selection
#'
#' Collapses the Pool to its currently selected Layout (one row per label) and runs
#' `forcePolish()`: pattern-search descent off the candidate grid that preserves the discrete
#' solution's conflict-freeness. Returns the polished Layout. `con_type` maps to the kernel's
#' integer leader style (0 = "cl" corner, else "cm"/"none").
#'
#' @param pool A Pool (see `candidatePool()`).
#' @param geom Box-fit structure (see `placeLabels()`).
#' @param xlim,ylim Numeric length-2 viewport bounds.
#' @param hw,hh Numeric per-label box half-sizes.
#' @param char_h Numeric line height (scales the polish step sizes).
#' @param pad Numeric hard box clearance.
#' @param con_type Leader style.
#' @return A Layout (one row per label) at the polished centres.
#' @keywords internal
#' @noRd
polishLayout <- function(pool, geom, xlim, ylim, hw, hh, char_h, pad, con_type) {
  poi <- geom$poi
  chosen <- pool$cand[pool$sel + 1L][order(label)]   # selected candidate per label
  MU <- 55                                     # box-spacing weight (fixed polish tuning)
  iters <- 120L                                # polish iteration count
  polished <- forcePolish(
    geom$rtree, chosen$cx, chosen$cy, hw, hh, poi[, 1], poi[, 2],
    geom$polysx, geom$polysy, pad, xlim[1], xlim[2], ylim[1], ylim[2],
    as.integer(iters), 0.4 * char_h, MU, 0.6 * char_h, 0.03 * char_h,
    if (con_type == "cl") 0L else 1L)
  .layoutFromCentres(
    data.table::data.table(label = seq_len(nrow(poi)), cx = polished$cx, cy = polished$cy),
    hw, hh, poi, pad, con_type)[order(label)]
}

#' Add the visible leader end to a Layout
#'
#' Fills the `bx`, `by` columns: where each leader (anchor -> pole) first meets the label's own
#' cluster ring, the point the drawn leader stops at (`firstLeaderHit()`).
#'
#' @param layout A Layout (one row per label; see `.layout()`).
#' @param geom Box-fit structure (for the cluster rings).
#' @return The Layout with `bx`, `by` added.
#' @keywords internal
#' @noRd
withLeaderEnds <- function(layout, geom) {
  poi <- geom$poi
  hit <- firstLeaderHit(layout$ex, layout$ey,
                        poi[layout$label, 1], poi[layout$label, 2],
                        geom$polysx[layout$label], geom$polysy[layout$label])
  layout[, `:=`(bx = hit$bx, by = hit$by)][]
}

#' Place cluster labels for one view
#'
#' Boundary-seed placement driver, run as a straight pipeline over the placement structures:
#' `candidatePool()` -> `oneMoveSweep()` -> `twoMoveSweep()` -> `polishLayout()` ->
#' `withLeaderEnds()`. Conflict-free by construction given a feasible pool. Pure given the
#' per-label box half-sizes and line height (the draw hook supplies those from text metrics in
#' the panel's mm space).
#'
#' @param geom Box-fit structure: a list with `poi` (K x 2 poles), `rtree` (`XPtr<BoxFit>`),
#'   `polysx`/`polysy` (per-cluster rings) and optionally `pad_xrange` (dilated x-extent).
#' @param xlim,ylim Numeric length-2 viewport bounds (already inset by `label.buffer`).
#' @param hw,hh Numeric per-label box half-sizes (mm).
#' @param char_h Numeric line height (mm) used to scale the internal spacing constants.
#' @param con_type Leader style: `"cl"`, `"cm"`, or `"none"`.
#' @return A Layout data.table (one row per cluster) with `cx`, `cy`, the box columns, the
#'   leader anchor `ex`, `ey`, its `corner` flag and the visible leader end `bx`, `by`.
#' @keywords internal
#' @noRd
placeLabels <- function(geom, xlim, ylim, hw, hh, char_h, con_type = "cl") {
  poi <- geom$poi
  K <- nrow(poi)
  pad <- 0.05 * char_h                         # hard box clearance
  gap <- 0.25 * char_h                         # seed column spacing

  if (K == 1) {
    # single label sits on its pole; the leader is degenerate (ends at the pole)
    only <- .layoutFromCentres(
      data.table::data.table(label = 1L, cx = poi[1, 1], cy = poi[1, 2]),
      hw, hh, poi, pad, con_type)
    return(only[, `:=`(bx = poi[1, 1], by = poi[1, 2])][])
  }

  pool <- candidatePool(geom, xlim, ylim, hw, hh, char_h, pad, gap, con_type)
  pool <- oneMoveSweep(pool)
  pool <- twoMoveSweep(pool)
  layout <- polishLayout(pool, geom, xlim, ylim, hw, hh, char_h, pad, con_type)
  withLeaderEnds(layout, geom)
}

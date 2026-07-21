# Anchor-leader label placement driver.
#
# The pipeline threads a single Layout object through a chain of transforms (see
# `placeLabels()`): placementScene() -> seedLayout() -> addRadialCandidates() -> oneMoveSweep()
# -> twoMoveSweep() -> polishLayout() -> withLeaderEnds(). The leader START on the label box is
# decoupled from the box centre (see `.anchorPoint`); the C++ kernels in src/ score conflicts on
# the actual anchor->pole segment (forcePolish mirrors `.anchorPoint`). This file is the R
# orchestration.
#
# Structures:
#   * scene (see `placementScene`) -- the fixed placement problem, built ONCE and then carried
#     inside every Layout so no stage takes the geometry as a separate argument: cluster
#     geometry (poles, box-fit R-tree, polygons, dilated x-range), the viewport, the per-label box
#     sizes and leader style, and the derived spacing constants.
#   * Layout (see `seedLayout`) -- the evolving placement state: `scene`, a `place` table of
#     placement rows tagged by `label` (one row per label to start; addRadialCandidates() appends
#     candidate alternatives, so a label may then own several rows), and `sel`, the 0-based row of
#     the currently selected placement per label. Its solution is `place[sel]`. Every stage after
#     the scene is a Layout -> Layout (or Layout -> result) transform.
#
# The final result is a one-row-per-label data.table (the selected placements plus the visible
# leader ends `bx`, `by`) -- the columns `mark_label.R` and the test scorer read.
#
# `con_type` selects the leader style: "ledge" (corners + a short horizontal ledge along the box
# edge, one sign-quadrant corner), "line" (corners + edge-midpoints, 8-point rule, no ledge),
# "box" (like "line", plus the label's bounding box outline), "none" (placed like "line", no
# connector drawn). Only "ledge" differs geometrically; the rest share the line anchor rule.
# The keep-out padding between labels and clusters is baked into the scene's box-fit R-tree: the
# caller passes polygons already dilated by `label.buffer` (see my_make_label).

# Canonical Layout placement columns (one placement per row). `.assemble()` attaches an `eff`
# (effective length) ranking column on top; `currentPlacements()` selects just these.
layoutCols <- c("label", "cx", "cy", "hw", "hh", "tx", "ty",
                "cxmin", "cxmax", "cymin", "cymax", "ex", "ey", "corner", "len")

#' Assemble the fixed placement scene
#'
#' Bundles everything that defines the placement problem and never changes during optimization,
#' so it is built once and carried inside every Layout. The derived spacing constant (`pad`, the
#' hard box clearance) is computed here from the line height.
#'
#' @param geom Box-fit structure: a list with `poi` (K x 2 poles), `rtree` (`XPtr<BoxFit>`),
#'   `polysx`/`polysy` (per-cluster polygons) and optionally `pad_xrange` (dilated x-extent).
#' @param xlim,ylim Numeric length-2 viewport bounds (already inset by `label.buffer`).
#' @param hw,hh Numeric per-label box half-sizes (mm).
#' @param char_h Numeric line height (mm) used to scale the internal spacing constants.
#' @param con_type Leader style: `"ledge"`, `"line"`, `"box"`, or `"none"`.
#' @return A scene list.
#' @keywords internal
#' @noRd
placementScene <- function(geom, xlim, ylim, hw, hh, char_h, con_type) {
  polyxlim <- geom$pad_xrange                   # dilated extent (keep-out)
  if (is.null(polyxlim)) {
    xs <- unlist(geom$polysx)                   # undilated geom (tests); guard the empty view
    polyxlim <- if (length(xs) > 0) range(xs) else c(0, 0)
  }
  list(
    poi = geom$poi, K = nrow(geom$poi),
    rtree = geom$rtree, polysx = geom$polysx, polysy = geom$polysy, polyxlim = polyxlim,
    xlim = xlim, ylim = ylim,
    hw = hw, hh = hh, char_h = char_h, con_type = con_type,
    pad = 0.05 * char_h)
}

#' Find the leader's start anchor on a label box
#'
#' The point on the box border where the leader begins, aimed at the pole (`tx`, `ty`). For
#' `con_type == "ledge"` this is always the sign(dx),sign(dy) quadrant corner. Otherwise (the
#' ggforce `get_end_points` rule) it is a corner only when the pole is fully to the side in
#' both axes; if the pole's x lies between the vertical edges it is the top/bottom edge
#' centre, and if its y lies between the horizontal edges it is the left/right edge centre.
#'
#' @param cx,cy Numeric box-centre coordinates (vectorised).
#' @param hw,hh Numeric box half-width / half-height.
#' @param tx,ty Numeric pole (leader target) coordinates.
#' @param con_type Leader style: `"ledge"`, `"line"`, `"box"`, or `"none"`.
#' @return A list with numeric `ex`, `ey` (the anchor) and logical `corner` (is it a corner?).
#' @keywords internal
#' @noRd
.anchorPoint <- function(cx, cy, hw, hh, tx, ty, con_type) {
  dx <- tx - cx
  dy <- ty - cy
  signX <- ifelse(dx >= 0, 1, -1)
  signY <- ifelse(dy >= 0, 1, -1)
  if (con_type == "ledge") {
    return(list(ex = cx + signX * hw, ey = cy + signY * hh,
                corner = rep(TRUE, length(cx))))
  }
  xInside <- abs(dx) < hw
  yInside <- abs(dy) < hh
  list(ex = ifelse(xInside & !yInside, cx, cx + signX * hw),
       ey = ifelse(!xInside & yInside, cy, cy + signY * hh),
       corner = !xInside & !yInside)
}

#' Build placement rows from explicit box centres and leader anchors
#'
#' The low-level `place`-row constructor: one row per label with the box centre, padded box
#' extents, pole, leader anchor and ranking length (centre-to-pole distance). The caller supplies
#' the leader anchor directly (the boundary seed needs anchors that do not follow the
#' `.anchorPoint()` rule). Use `.layoutFromCentres()` when the anchor should follow the style.
#'
#' @param scene The placement scene (for `hw`, `hh`, `poi`, `pad`).
#' @param label Integer 1-indexed cluster of each row.
#' @param cx,cy Numeric box-centre coordinates.
#' @param ex,ey Numeric leader-start (anchor) coordinates.
#' @param corner Logical: is the anchor a box corner?
#' @return A data.table of placement rows with the `layoutCols` columns.
#' @keywords internal
#' @noRd
.layout <- function(scene, label, cx, cy, ex, ey, corner) {
  halfW <- scene$hw[label]
  halfH <- scene$hh[label]
  tx <- scene$poi[label, 1]
  ty <- scene$poi[label, 2]
  pad <- scene$pad
  data.table::data.table(
    label = label, cx = cx, cy = cy, hw = halfW, hh = halfH, tx = tx, ty = ty,
    cxmin = cx - halfW - pad, cxmax = cx + halfW + pad,
    cymin = cy - halfH - pad, cymax = cy + halfH + pad,
    ex = ex, ey = ey, corner = corner,
    len = sqrt((cx - tx)^2 + (cy - ty)^2))
}

#' Build placement rows from box centres, deriving leader anchors from the leader style
#'
#' Convenience wrapper over `.layout()`: derives each leader anchor from the box centre and
#' pole via `.anchorPoint()`, then assembles the placement rows.
#'
#' @param scene The placement scene.
#' @param dt A data.table/data.frame with columns `label`, `cx`, `cy`.
#' @return A data.table of placement rows (see `.layout()`).
#' @keywords internal
#' @noRd
.layoutFromCentres <- function(scene, dt) {
  anchor <- .anchorPoint(dt$cx, dt$cy, scene$hw[dt$label], scene$hh[dt$label],
                         scene$poi[dt$label, 1], scene$poi[dt$label, 2], scene$con_type)
  .layout(scene, dt$label, dt$cx, dt$cy, anchor$ex, anchor$ey, anchor$corner)
}

#' Wrap placement rows into a ranked Layout
#'
#' Computes each row's effective length (`eff`, the quantity the sweeps minimise; see
#' `effectiveLength()`), orders the rows by `(label, eff)` so each label's candidate group is
#' contiguous and ascending in that ranking -- the two-move branch-and-bound prunes on it and so
#' requires the ordering to match -- and records the current selection: `sel[i]` is the 0-based
#' row of label i's active placement (the one true `active` row per label).
#'
#' @param scene The placement scene.
#' @param place A data.table of placement rows (see `.layout()`).
#' @param active Logical, one per row: the currently selected placement of each label.
#' @return A Layout list `(scene, place, sel)`, with `place` ranked and carrying `eff`.
#' @keywords internal
#' @noRd
.assemble <- function(scene, place, active) {
  place$eff <- effectiveLength(
    place$len, place$ex, place$ey, place$tx, place$ty, as.integer(place$label),
    scene$polysx, scene$polysy, place$cxmin, place$cxmax, place$cymin, place$cymax,
    scene$xlim[1], scene$xlim[2], scene$ylim[1], scene$ylim[2])
  ord <- order(place$label, place$eff)          # rank by the metric the kernels prune on
  place <- place[ord]
  active <- active[ord]
  sel <- vapply(seq_len(scene$K), function(i) {
    which(active & place$label == i)[1] - 1L    # 0-based row of label i's active placement
  }, 0L)
  list(scene = scene, place = place, sel = sel)
}

#' Get a Layout's current one-row-per-label solution
#'
#' The selected placement of every label, in label order, as the canonical `layoutCols` columns
#' (dropping any transient optimizer columns such as `eff`).
#'
#' @param layout A Layout.
#' @return A data.table with one row per label.
#' @keywords internal
#' @noRd
currentPlacements <- function(layout) {
  layout$place[layout$sel + 1L, layoutCols, with = FALSE]
}

#' Stack one boundary side column of label slots
#'
#' Places the labels assigned to one side (`side` -1 left, +1 right) as a vertical stack on the
#' line `x = Xline`. A column may be empty (0 labels) or hold a single label, in which case the
#' label just sits at its pole y. Otherwise the labels are placed on a uniform grid of slots as
#' tall as the tallest box and a single Hungarian assigns each label to a slot -- one box per
#' slot, so any assignment leaves the boxes at most touching (single- and multi-line alike).
#'
#' * If the whole column fits the viewport, the grid has a few slack slots and is clamped into
#'   the viewport, so labels sit near their poles.
#' * Otherwise the column is too crowded for slack: the grid is exactly `m` slots centred on the
#'   pole span, extended past the viewport when the labels still cannot all fit.
#'
#' @param scene The placement scene (for `poi`, `hh`, `ylim`).
#' @param set Integer labels assigned to this column.
#' @param Xline Numeric x-coordinate of the column line.
#' @param side Integer -1 (left column) or +1 (right column).
#' @return A data.table with `label`, `cy`, `Xline`, `side`.
#' @keywords internal
#' @noRd
.sideColumn <- function(scene, set, Xline, side) {
  poi <- scene$poi
  boxH <- 2 * scene$hh                          # full box height per label
  ylim <- scene$ylim
  m <- length(set)
  poleY <- poi[set, 2]
  if (m <= 1) {
    # empty column (0), or a lone label that just sits at its pole y (nothing to stack)
    return(data.table::data.table(label = set, cy = poleY,
                                  Xline = rep_len(Xline, m), side = rep_len(side, m)))
  }
  dx2 <- (Xline - poi[set, 1])^2                # squared horizontal pole-to-column distance
  viewH <- ylim[2] - ylim[1]
  slotH <- max(boxH[set])                       # uniform slot as tall as the tallest box (no pad)
  capacity <- floor(viewH / slotH)              # such slots the viewport holds

  # Lay a grid of uniform tallest-box slots centred on the pole span. When the whole column
  # fits, use a few slack slots and clamp the grid into the viewport so labels sit near their
  # poles; when it is too crowded to fit, use exactly `m` slots and let the grid extend past
  # the viewport. `min(capacity, ...)` collapses to `m` in the crowded case (capacity < m).
  spanSlots <- ceiling((max(poleY) - min(poleY)) / slotH)
  nSlot <- max(m, min(capacity, spanSlots + 2L * m))
  lo <- mean(range(poleY)) - nSlot * slotH / 2
  if (nSlot * slotH <= viewH) {                 # clamp into the viewport only when it fits
    lo <- max(ylim[1], min(lo, ylim[2] - nSlot * slotH))
  }
  slotY <- lo + (seq_len(nSlot) - 0.5) * slotH

  # square cost: m real label rows (leader length to each slot) + dummy zero-cost rows. One box
  # per slot, so the Hungarian assignment leaves the boxes at most touching.
  cost <- matrix(0, nSlot, nSlot)
  for (i in seq_len(m)) {
    cost[i, ] <- sqrt(dx2[i] + (poleY[i] - slotY)^2)
  }
  assignment <- hungarian(cost) + 1L
  data.table::data.table(label = set, cy = slotY[assignment[seq_len(m)]],
                         Xline = Xline, side = side)
}

#' Build the conflict-free boundary side slots
#'
#' The starting placement: two columns hang off the cluster cloud on the vertical lines
#' `x = min/max` polygon x. Because the polygons are dilated by `label.buffer`, those lines sit
#' outside the true clusters, so a box resting on its line never clips a cluster. Labels are
#' split left/right by pole x (balancing the label count) and each column is stacked by
#' `.sideColumn()`.
#'
#' @param scene The placement scene.
#' @return A data.table with `label`, `cx`, `cy`, `ex`, `ey`, `corner`.
#' @keywords internal
#' @noRd
.sideSlots <- function(scene) {
  poi <- scene$poi
  K <- scene$K
  hw <- scene$hw
  hh <- scene$hh

  # split labels into a left / right column by pole x, balancing the label count
  byX <- order(poi[, 1])
  split <- max(1L, min(K %/% 2L, K - 1L))
  leftSet <- byX[seq_len(split)]
  rightSet <- byX[-seq_len(split)]             # the rest (empty when split == K, i.e. K == 1)

  slots <- rbind(
    .sideColumn(scene, leftSet, scene$polyxlim[1], -1),
    .sideColumn(scene, rightSet, scene$polyxlim[2], +1))

  # box centre so the padded near edge sits on the column line; the leader starts on the true
  # (unpadded) near edge so a "ledge" connector's ledge meets it. The box stays a padded-margin
  # width inside the line; column neighbours sit one tallest-box apart, so they at most touch.
  slots[, `:=`(
    cx = Xline + side * (hw[label] + scene$pad),
    ex = Xline + side * scene$pad,
    ey = if (scene$con_type == "ledge") cy - hh[label] else cy,
    corner = scene$con_type == "ledge")]
  slots[]
}

#' Build the initial feasible layout (the seed)
#'
#' Wraps the conflict-free boundary side slots (`.sideSlots()`) into a Layout -- one placement
#' per label, all selected. This is the starting state the rest of the pipeline improves.
#'
#' @param scene The placement scene.
#' @return A Layout.
#' @keywords internal
#' @noRd
seedLayout <- function(scene) {
  slots <- .sideSlots(scene)
  place <- .layout(scene, slots$label, slots$cx, slots$cy,
                   slots$ex, slots$ey, slots$corner)
  .assemble(scene, place, active = rep(TRUE, nrow(place)))
}

#' Generate radial free-space slot candidates
#'
#' From each pole, `radialCandidates()` marches rays outward and emits a box centre wherever a
#' padded box is cluster-free; this turns those centres into placement rows.
#'
#' @param scene The placement scene.
#' @return A data.table of placement rows (many per label).
#' @keywords internal
#' @noRd
.radialSlots <- function(scene) {
  ndir <- 48L
  radStep <- 0.3 * scene$char_h
  radStart <- 0.2 * scene$char_h
  radReach <- 16 * scene$char_h
  radFill <- 1.2 * scene$char_h
  dedup <- 0.3 * scene$char_h
  cand <- data.table::as.data.table(radialCandidates(
    scene$rtree, scene$poi, scene$hw, scene$hh, scene$pad,
    ndir, radStep, radStart, radReach, radFill, dedup))
  .layoutFromCentres(scene, cand)
}

#' Append radial candidate placements to a Layout
#'
#' Adds the radial free-space candidates (`.radialSlots()`) to the layout's current
#' placements as candidate alternatives, keeping the current selection active, and re-ranks
#' (`.assemble()`). After this a label may own several rows; the sweeps pick among them.
#'
#' @param layout A Layout.
#' @return The Layout with candidate rows appended and ranked.
#' @keywords internal
#' @noRd
addRadialCandidates <- function(layout) {
  scene <- layout$scene
  current <- currentPlacements(layout)          # the selected placement of each label
  radial <- .radialSlots(scene)
  place <- rbind(current, radial)
  active <- c(rep(TRUE, nrow(current)), rep(FALSE, nrow(radial)))
  .assemble(scene, place, active)
}

#' List each label's 0-based candidate rows
#'
#' `rows[[i]]` lists the 0-based `place` rows available to label i -- the candidate index lists
#' the C++ sweep kernels search over.
#'
#' @param layout A Layout.
#' @return A list of integer vectors, one per label.
#' @keywords internal
#' @noRd
.candidateRows <- function(layout) {
  unname(split(seq_len(nrow(layout$place)) - 1L, layout$place$label))
}

#' Run the one-move conflict sweep
#'
#' Wraps `oneMoveSweepKernel()`: moves each label to its lexicographically best candidate versus
#' the current others (including staying put), ranking by effective length. Updates the Layout's
#' selection.
#'
#' @param layout A Layout with ranked candidates (see `addRadialCandidates()`).
#' @param maxpass Integer cap on the number of sweeps.
#' @return The Layout with `sel` updated.
#' @keywords internal
#' @noRd
oneMoveSweep <- function(layout, maxpass = 100L) {
  p <- layout$place
  layout$sel <- oneMoveSweepKernel(
    p$cxmin, p$cxmax, p$cymin, p$cymax, p$ex, p$ey, p$tx, p$ty, p$eff,
    .candidateRows(layout), as.integer(layout$sel), maxpass)
  layout
}

#' Run the two-move conflict / length refinement
#'
#' Wraps `twoMoveSweepKernel()`: per label, applies the best two-move (length branch-and-bound
#' when the label is conflict-free, all-pairs search when it is conflicted). Updates the Layout's
#' selection.
#'
#' @param layout A Layout with ranked candidates (see `addRadialCandidates()`).
#' @param maxpass Integer cap on the number of passes.
#' @param sq Logical; if `TRUE` the length objective uses the squared length.
#' @return The Layout with `sel` updated.
#' @keywords internal
#' @noRd
twoMoveSweep <- function(layout, maxpass = 50L, sq = TRUE) {
  p <- layout$place
  layout$sel <- twoMoveSweepKernel(
    p$cxmin, p$cxmax, p$cymin, p$cymax, p$ex, p$ey, p$tx, p$ty, p$eff,
    .candidateRows(layout), as.integer(layout$sel), maxpass, sq)
  layout
}

#' Run the continuous force-directed polish
#'
#' Runs `forcePolish()` on the layout's selected placements: pattern-search descent off the
#' candidate grid that preserves the discrete solution's conflict-freeness. Returns a Layout
#' whose placements are the polished centres. `con_type` maps to the kernel's integer leader
#' style (0 = "ledge" corner, else "line"/"box"/"none" -- all share the line anchor rule).
#'
#' @param layout A Layout.
#' @return A Layout at the polished centres (one placement per label).
#' @keywords internal
#' @noRd
polishLayout <- function(layout) {
  scene <- layout$scene
  chosen <- currentPlacements(layout)
  MU <- 55                                      # box-spacing weight (fixed polish tuning)
  iters <- 120L                                 # polish iteration count
  polished <- forcePolish(
    scene$rtree, chosen$cx, chosen$cy, scene$hw, scene$hh, scene$poi[, 1], scene$poi[, 2],
    scene$polysx, scene$polysy, scene$pad,
    scene$xlim[1], scene$xlim[2], scene$ylim[1], scene$ylim[2],
    as.integer(iters), 0.4 * scene$char_h, MU, 0.6 * scene$char_h, 0.03 * scene$char_h,
    if (scene$con_type == "ledge") 0L else 1L)
  place <- .layoutFromCentres(
    scene, data.table::data.table(label = seq_len(scene$K),
                                  cx = polished$cx, cy = polished$cy))
  .assemble(scene, place, active = rep(TRUE, scene$K))
}

#' Add the visible leader ends and return the final placement table
#'
#' Collapses the Layout to its one-row-per-label solution and fills `bx`, `by`: where each leader
#' (anchor -> pole) first meets the label's own cluster polygon, the point the drawn leader stops at
#' (`firstLeaderHit()`). This is the deliverable the draw stage consumes.
#'
#' @param layout A Layout.
#' @return A data.table (one row per label) with the placement columns plus `bx`, `by`.
#' @keywords internal
#' @noRd
withLeaderEnds <- function(layout) {
  scene <- layout$scene
  sol <- currentPlacements(layout)
  hit <- firstLeaderHit(sol$ex, sol$ey,
                        scene$poi[sol$label, 1], scene$poi[sol$label, 2],
                        scene$polysx[sol$label], scene$polysy[sol$label])
  sol$bx <- hit$bx
  sol$by <- hit$by
  sol[]
}

#' Build the empty placement for a view with no clusters
#'
#' @param scene The placement scene (with `K == 0`).
#' @return A 0-row placement data.table with the full column set.
#' @keywords internal
#' @noRd
emptyPlacement <- function(scene) {
  place <- .layoutFromCentres(
    scene, data.table::data.table(label = integer(0), cx = numeric(0), cy = numeric(0)))
  place$bx <- numeric(0)
  place$by <- numeric(0)
  place[]
}

#' Place cluster labels for one view
#'
#' Boundary-seed placement driver, run as a straight pipeline over a single Layout object:
#' `seedLayout()` -> `addRadialCandidates()` -> `oneMoveSweep()` -> `twoMoveSweep()` ->
#' `polishLayout()` -> `withLeaderEnds()`. Conflict-free by construction given a feasible seed.
#' Pure given the per-label box half-sizes and line height (the draw hook supplies those from text
#' metrics in the panel's mm space). It returns one placement per cluster, so an empty view
#' (K == 0) returns an empty layout.
#'
#' @param geom Box-fit structure (see `placementScene()`).
#' @param xlim,ylim Numeric length-2 viewport bounds (already inset by `label.buffer`).
#' @param hw,hh Numeric per-label box half-sizes (mm).
#' @param char_h Numeric line height (mm) used to scale the internal spacing constants.
#' @param con_type Leader style: `"ledge"`, `"line"`, `"box"`, or `"none"`.
#' @return A data.table (one row per cluster) with `cx`, `cy`, the box columns, the leader anchor
#'   `ex`, `ey`, its `corner` flag and the visible leader end `bx`, `by`.
#' @keywords internal
#' @noRd
placeLabels <- function(geom, xlim, ylim, hw, hh, char_h, con_type = "ledge") {
  scene <- placementScene(geom, xlim, ylim, hw, hh, char_h, con_type)
  # No clusters -> no placements. The early return only spares the pipeline an empty seed
  # (`seedLayout()` assumes at least one label); the contract is the same one-row-per-cluster.
  if (scene$K == 0) {
    return(emptyPlacement(scene))
  }

  layout <- seedLayout(scene)
  layout <- addRadialCandidates(layout)
  layout <- oneMoveSweep(layout)
  layout <- twoMoveSweep(layout)
  layout <- polishLayout(layout)
  withLeaderEnds(layout)
}

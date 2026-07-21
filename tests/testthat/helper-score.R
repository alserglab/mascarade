# Test-only layout scoring (not shipped in the package). Computes the conflict counts
# (bb box-box, ll leader-leader, lb leader-box) plus leader length from a placement table
# with the geometry columns produced by placeLabels(). scoreBetter() is the golden gate.

.vsegcross <- function(ax, ay, bx, by, cx, cy, dx, dy) {
  d1x <- bx - ax; d1y <- by - ay; d2x <- dx - cx; d2y <- dy - cy
  den <- d1x * d2y - d1y * d2x
  t <- ((cx - ax) * d2y - (cy - ay) * d2x) / den
  u <- ((cx - ax) * d1y - (cy - ay) * d1x) / den
  ok <- abs(den) > 1e-12
  ok & t > 1e-9 & t < 1 - 1e-9 & u > 1e-9 & u < 1 - 1e-9
}

.vsegbox <- function(ax, ay, bx, by, x0, x1, y0, y1) {
  (ax > x0 & ax < x1 & ay > y0 & ay < y1) | (bx > x0 & bx < x1 & by > y0 & by < y1) |
    .vsegcross(ax, ay, bx, by, x0, y0, x1, y0) | .vsegcross(ax, ay, bx, by, x1, y0, x1, y1) |
    .vsegcross(ax, ay, bx, by, x1, y1, x0, y1) | .vsegcross(ax, ay, bx, by, x0, y1, x0, y0)
}

layoutScore <- function(P) {
  n <- nrow(P)
  if (n < 2) return(c(bb = 0, ll = 0, lb = 0, len = sum(P$len), maxlen = if (n) max(P$len) else 0))
  cp <- utils::combn(n, 2); a <- cp[1, ]; b <- cp[2, ]
  bb <- sum(P$cxmin[a] < P$cxmax[b] & P$cxmin[b] < P$cxmax[a] & P$cymin[a] < P$cymax[b] & P$cymin[b] < P$cymax[a])
  lb <- sum(.vsegbox(P$ex[a], P$ey[a], P$tx[a], P$ty[a], P$cxmin[b], P$cxmax[b], P$cymin[b], P$cymax[b]) |
            .vsegbox(P$ex[b], P$ey[b], P$tx[b], P$ty[b], P$cxmin[a], P$cxmax[a], P$cymin[a], P$cymax[a]))
  ll <- sum(.vsegcross(P$ex[a], P$ey[a], P$tx[a], P$ty[a], P$ex[b], P$ey[b], P$tx[b], P$ty[b]))
  c(bb = bb, ll = ll, lb = lb, len = sum(P$len), maxlen = max(P$len))
}

# strictly better? feasibility-first, then length (golden reference gate)
scoreBetter <- function(a, b) {
  for (k in c("bb", "ll", "lb", "len")) if (a[k] != b[k]) return(a[k] < b[k])
  FALSE
}

# Build a placement geom (poles + box-fit R-tree) from data + clusters, plus data-unit box
# sizes, for a headless placeLabels() call. char_h scales the boxes relative to the data.
.buildTestGeom <- function(dims, clusters, char_frac = 0.045) {
  mt <- as.data.frame(generateMask(dims = dims, clusters = clusters))
  xv <- colnames(mt)[1]; yv <- colnames(mt)[2]
  cl <- if (is.factor(mt$cluster)) levels(mt$cluster) else unique(mt$cluster)
  cl <- cl[vapply(cl, function(c) any(mt$cluster == c & grepl("#1$", mt$part)), logical(1))]
  polys <- lapply(cl, function(c) { s <- mt[mt$cluster == c & grepl("#1$", mt$part), ]; list(x = s[[xv]], y = s[[yv]]) })
  px <- lapply(polys, `[[`, "x"); py <- lapply(polys, `[[`, "y")
  poi <- t(vapply(polys, function(p) { pl <- polylabelr::poi(p$x, p$y); c(pl$x, pl$y) }, numeric(2)))
  xr <- range(unlist(px)); yr <- range(unlist(py))
  char_h <- char_frac * diff(yr)
  hh <- rep(char_h / 2, length(cl))
  hw <- nchar(cl) * (0.55 * char_h) / 2
  list(geom = list(poi = poi, rtree = buildBoxFit(px, py), polysx = px, polysy = py),
       hw = hw, hh = hh, char_h = char_h,
       xlim = xr + c(-1, 1) * 0.45 * diff(xr),
       ylim = yr + c(-1, 1) * 0.10 * diff(yr))
}

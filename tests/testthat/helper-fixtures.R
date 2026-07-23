# Test-only fixture construction and loading (not shipped in the package).
#
# The placement geometry is frozen into `fixtures/placement_geom.rds` (by
# `fixtures/make-placement-fixtures.R`) so the placement tests do not depend on
# `generateMask()`, spatstat or polylabelr. The box-fit R-tree is an external
# pointer that cannot be serialised, so the fixture stores plain data and
# `restorePlacementGeom()` rebuilds the tree on load.

#' Extract each cluster's largest (`#1`) part from a mask table
#'
#' @param maskTable A mask table as returned by `generateMask()`.
#' @return A named list of `list(x =, y =)` rings, one per cluster.
largestPartPolygons <- function(maskTable) {
  maskTable <- as.data.frame(maskTable)
  xColumn <- colnames(maskTable)[1]
  yColumn <- colnames(maskTable)[2]

  isLargestPart <- grepl("#1$", maskTable$part)
  clusters <- if (is.factor(maskTable$cluster)) {
    levels(maskTable$cluster)
  } else {
    unique(maskTable$cluster)
  }
  clusters <- clusters[vapply(clusters, function(cluster) {
    any(maskTable$cluster == cluster & isLargestPart)
  }, logical(1))]

  polygons <- lapply(clusters, function(cluster) {
    ring <- maskTable[maskTable$cluster == cluster & isLargestPart, ]
    list(x = ring[[xColumn]], y = ring[[yColumn]])
  })
  names(polygons) <- clusters
  polygons
}

#' Build a serialisable placement scene from a mask table
#'
#' Everything `placeLabels()` needs except the R-tree (see file header). Label
#' boxes are sized from the character count, not real font metrics, keeping the
#' fixture independent of the installed fonts.
#'
#' @param maskTable A mask table as returned by `generateMask()`.
#' @param charFrac Label character height, as a fraction of the data y-range.
#' @return A list of plain vectors and lists, safe to `saveRDS()`.
buildPlacementGeom <- function(maskTable, charFrac = 0.045) {
  polygons <- largestPartPolygons(maskTable)
  labels <- names(polygons)

  polysx <- lapply(polygons, `[[`, "x")
  polysy <- lapply(polygons, `[[`, "y")

  poles <- t(vapply(polygons, function(polygon) {
    pole <- polylabelr::poi(polygon$x, polygon$y)
    c(pole$x, pole$y)
  }, numeric(2)))
  dimnames(poles) <- NULL

  xRange <- range(unlist(polysx))
  yRange <- range(unlist(polysy))
  charHeight <- charFrac * diff(yRange)

  list(
    labels = labels,
    polysx = unname(polysx),
    polysy = unname(polysy),
    poi = poles,
    charHeight = charHeight,
    halfHeight = rep(charHeight / 2, length(labels)),
    halfWidth = nchar(labels) * (0.55 * charHeight) / 2,
    xlim = xRange + c(-1, 1) * 0.45 * diff(xRange),
    ylim = yRange + c(-1, 1) * 0.10 * diff(yRange)
  )
}

#' Rebuild the runtime placement scene from a frozen one
#'
#' Adds back the box-fit R-tree in the shape `placeLabels()` expects.
#'
#' @param frozen A scene from `buildPlacementGeom()`.
#' @return `frozen` with a `geom` element added, ready for `placeLabels()`.
restorePlacementGeom <- function(frozen) {
  frozen$geom <- list(
    poi = frozen$poi,
    rtree = buildBoxFit(frozen$polysx, frozen$polysy),
    polysx = frozen$polysx,
    polysy = frozen$polysy
  )
  frozen
}

#' Load a frozen placement scene by name (`"example"` or `"singleCluster"`)
#'
#' @param name Scene name in `fixtures/placement_geom.rds`.
#' @return A scene ready to pass to `placeLabels()`.
loadPlacementScene <- function(name = "example") {
  scenes <- readRDS(test_path("fixtures", "placement_geom.rds"))
  restorePlacementGeom(scenes[[name]])
}

#' Load a recorded golden layout score by name
#'
#' @param name Entry name in `fixtures/golden_scores.rds`.
#' @return A named score vector, as produced by `layoutScore()`.
loadGoldenScore <- function(name = "example_default") {
  readRDS(test_path("fixtures", "golden_scores.rds"))[[name]]
}

#' The shipped example mask table, optionally restricted to some clusters
#'
#' @param clusters Cluster names to keep; `NULL` (default) keeps all. Unused
#'   factor levels are dropped.
#' @return A mask table.
exampleMask <- function(clusters = NULL) {
  data("exampleMaskTable", package = "mascarade", envir = environment())
  maskTable <- exampleMaskTable
  if (!is.null(clusters)) {
    maskTable <- maskTable[maskTable$cluster %in% clusters, ]
    maskTable$cluster <- droplevels(maskTable$cluster)
  }
  maskTable
}

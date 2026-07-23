# Regenerates the placement test fixtures. Not run during `devtools::test()` --
# run it by hand from the package root when the placement geometry or the golden
# layout scores need re-recording:
#
#   source("tests/testthat/fixtures/make-placement-fixtures.R")
#
# Writes, next to this script:
#   placement_geom.rds  frozen placement geometry derived from `exampleMaskTable`
#   golden_scores.rds   layout scores of `placeLabels()` on that geometry

devtools::load_all(quiet = TRUE)

fixtureDir <- "tests/testthat/fixtures"
stopifnot(dir.exists(fixtureDir))

# Helpers are not on the search path outside a test run; source them so the
# fixture is built by the same code that loads it.
source("tests/testthat/helper-score.R")
source("tests/testthat/helper-fixtures.R")

data("exampleMaskTable", package = "mascarade")

firstCluster <- levels(exampleMaskTable$cluster)[1]
placementGeom <- list(
  example = buildPlacementGeom(exampleMaskTable),
  singleCluster = buildPlacementGeom(
    exampleMaskTable[exampleMaskTable$cluster == firstCluster, ]
  )
)
saveRDS(placementGeom, file.path(fixtureDir, "placement_geom.rds"))

scene <- restorePlacementGeom(placementGeom$example)
placement <- placeLabels(scene$geom, scene$xlim, scene$ylim,
                         scene$halfWidth, scene$halfHeight, scene$charHeight)
goldenScores <- list(example_default = layoutScore(placement))
saveRDS(goldenScores, file.path(fixtureDir, "golden_scores.rds"))

print(goldenScores)

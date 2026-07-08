# mascarade 0.4.0
* New leader-based label placement: cluster labels are placed by a boundary-seed
  optimizer (min-cost seed, one-/two-move refinement, force-directed polish) that
  keeps label boxes and leaders free of overlaps and crossings, replacing the
  previous label positioning in `geom_mark_shape()`. Placement runs at draw time, so
  labels re-place when the view or device size changes. Adds compiled code
  (`LinkingTo: Rcpp, BH`) and depends on `polylabelr` for cluster poles.

# mascarade 0.3.5
* bug fixes

# mascarade 0.3.4
* bug fixes

# mascarade 0.3.3
* option in fancyMask to label only the largest part of a cluster (default)
* bug fixes

# mascarade 0.3.1
* bug fixes

# mascarade 0.3.1
* fancyMask support setting colors manually or to inherit from the main layer.

# mascarade 0.3.0

* Significant performance improvements
* Introduced `generateMaskSeurat()` and `fancyMask()` helper functions.

# mascarade 0.2.0

* Major rewrite of generateMask function, with a change of interface

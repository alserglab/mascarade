# Build the label + leader grobs for a mark

Draw-time worker for `makeContent.shape_enc()`: dilates the cluster
polygons by `buffer` (the box keep-out), calls
[`my_place_labels()`](https://alserglab.github.io/mascarade/reference/my_place_labels.md)
for the placement, positions the label box grobs and builds the leader
polylines (anchor -\> visible mask-boundary end, plus the horizontal
ledge for `con_type == "ledge"`).

## Usage

``` r
my_make_label(
  labels,
  dims,
  polygons,
  ghosts,
  buffer,
  con_type,
  con_cap,
  con_gp,
  arrow,
  simp_ratio = 0.001
)
```

## Arguments

- labels:

  List of label-box grobs (one per mark part).

- dims:

  List of measured label box sizes `c(w, h)` in mm.

- polygons:

  List of cluster polygons (`list(x, y)`) in mm.

- ghosts:

  Points to avoid (currently unused by the placer).

- buffer:

  Grid unit: the `label.buffer` polygon padding / box keep-out.

- con_type:

  Leader style: `"ledge"`, `"line"`, `"box"`, or `"none"`.

- con_cap:

  Numeric gap (mm) left between the leader end and the cluster.

- con_gp:

  A `gpar` for the connectors (per drawn label).

- arrow:

  Optional [`grid::arrow`](https://rdrr.io/r/grid/arrow.html) for the
  connectors.

- simp_ratio:

  Numeric polygon-simplification fraction (see
  [`simplify_outer()`](https://alserglab.github.io/mascarade/reference/simplify_outer.md)).

## Value

A `gList`-ready list: the positioned label grobs followed by the
connector grob.

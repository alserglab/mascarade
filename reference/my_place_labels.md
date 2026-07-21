# Place labels at draw time (boundary-seed placer)

Runs at draw time in the panel's millimetre space: builds the poles and
the box-fit R-tree from the cluster polygons (cheap, ~20 ms for ~40
clusters; the expensive mask is not recomputed) and calls
`placeLabels()`. The box-fit keep-out uses the dilated polygons while
poles, leader ends and foreign-routing use the true ones, so leaders
reach the real outline.

## Usage

``` r
my_place_labels(
  rects,
  polygons,
  polygons_pad,
  bounds,
  simp_ratio = 0.001,
  con_type = "ledge",
  buffer = 0,
  hardpad = 0,
  softpad = 0
)
```

## Arguments

- rects:

  List of measured label box sizes `c(w, h)` in mm (a zeroed entry = not
  drawn).

- polygons:

  List of true cluster polygons (`list(x, y)`) in mm.

- polygons_pad:

  List of the same polygons dilated by `label.buffer` (the box
  keep-out).

- bounds:

  Numeric `c(width, height)` of the panel in mm.

- simp_ratio:

  Numeric polygon-simplification fraction (see
  [`simplify_outer()`](https://alserglab.github.io/mascarade/reference/simplify_outer.md)).

- con_type:

  Leader style: `"ledge"`, `"line"`, `"box"`, or `"none"`.

- buffer:

  Numeric `label.buffer` in mm; the overflow viewport is inset by it.

- hardpad:

  Numeric `label.hardpad` in mm: hard box clearance folded into every
  placement rectangle (seed slots, sweeps and polish alike).

- softpad:

  Numeric `label.softpad` in mm: extra target box spacing the polish
  aims for, on top of `hardpad`.

## Value

A list, one entry per input label: the placed centre `c(x, y)` in mm
(`NULL` if not drawn), carrying `attr(., "leaders")` with
`c(ex, ey, bx, by, corner)` per drawn label.

## Details

Note: each cluster polygon here is a single `list(x, y)` ring; any mask
holes are resolved upstream in
[`generateMask()`](https://alserglab.github.io/mascarade/reference/generateMask.md),
so this layer treats every polygon as one simple ring.

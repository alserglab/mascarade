# Generate ggplot2 layers for a labeled cluster mask

Convenience helper that returns a list of ggplot2 components that draws
polygon-like outlines and places cluster labels. The plotting limits are
expanded (via `limits.expand`) to provide extra room for labels.

## Usage

``` r
fancyMask(
  maskTable,
  ratio = NULL,
  limits.expand = ifelse(label, 0.1, 0.05),
  linewidth = 1,
  shape.expand = linewidth * unit(-1, "pt"),
  cols = "auto",
  label = TRUE,
  label.largest = TRUE,
  label.fontsize = 10,
  label.buffer = unit(2, "mm"),
  label.fontface = "plain",
  label.margin = margin(2, 2, 2, 2, "pt"),
  label.width = NULL,
  simp_ratio = 0.001,
  con.type = "ledge"
)
```

## Arguments

- maskTable:

  A data.frame of mask coordinates. The first two columns are
  interpreted as x/y coordinates (in that order). Must contain at least
  the columns `cluster` (a factor) and `group` (grouping identifier
  passed to
  [`geom_mark_shape()`](https://alserglab.github.io/mascarade/reference/geom_mark_shape.md)).

- ratio:

  Optional aspect ratio passed to
  [`ggplot2::coord_cartesian()`](https://ggplot2.tidyverse.org/reference/coord_cartesian.html).
  Use `1` for equal scaling. Default is `NULL` (no fixed ratio).

- limits.expand:

  Numeric scalar giving the fraction of the x/y range to expand on both
  sides when setting plot limits. Default is `0.1` with labels and 0.05
  with no labels.

- linewidth:

  Line width passed to
  [`geom_mark_shape()`](https://alserglab.github.io/mascarade/reference/geom_mark_shape.md)
  for the outline. Default is `1`.

- shape.expand:

  Expansion or contraction applied to the marked shapes, passed to
  `geom_mark_shape(expand = ...)`. Default is `unit(-linewidth, "pt")`.

- cols:

  Color specification for cluster outlines (and labels). One of:

  - `"auto"` (default) — inspects the plot at the time `fancyMask()` is
    added with `+`. If a layer maps `colour` to a discrete (non-numeric)
    variable, the mask joins that scale via `aes(colour = cluster)` so
    colours stay in sync regardless of `scale_color_*()` order.
    Otherwise (continuous colour, constant colour, or no colour
    aesthetic) explicit colours from
    [`scales::hue_pal()`](https://scales.r-lib.org/reference/pal_hue.html)
    are baked in and the plot's scale system is left untouched.

  - `"inherit"` — always maps `colour` as an aesthetic
    (`aes(colour = cluster)`), unconditionally joining whatever colour
    scale is present. Useful when you want to force scale sharing; will
    error if the existing scale is continuous.

  - A palette function that accepts a single integer `n` and returns `n`
    colors (e.g.,
    [`scales::hue_pal()`](https://scales.r-lib.org/reference/pal_hue.html),
    `rainbow`).

  - A single color string — applied to every cluster.

  - An unnamed character vector of length equal to the number of
    clusters — colors are assigned to clusters in factor-level order.

  - A named character vector — names must match cluster levels; order
    does not matter.

- label:

  Boolean flag whether the labels should be displayed.

- label.largest:

  Boolean flag. When `TRUE` (default), only the largest part of each
  cluster is labelled; smaller disconnected parts are drawn but not
  labelled. When `FALSE`, all parts are labelled. Ignored when
  `label = FALSE`.

- label.fontsize:

  Label font size passed to
  [`geom_mark_shape()`](https://alserglab.github.io/mascarade/reference/geom_mark_shape.md).
  Default is `10`.

- label.buffer:

  Polygon padding passed to
  [`geom_mark_shape()`](https://alserglab.github.io/mascarade/reference/geom_mark_shape.md):
  cluster polygons are dilated by this distance and labels are kept out
  of the dilated zone, so each label keeps a gap from its cluster
  outline. Default `unit(2, "mm")`; `unit(0, "mm")` disables.

- label.fontface:

  Label font face passed to
  [`geom_mark_shape()`](https://alserglab.github.io/mascarade/reference/geom_mark_shape.md).
  Default is `"plain"`.

- label.margin:

  Label margin passed to
  [`geom_mark_shape()`](https://alserglab.github.io/mascarade/reference/geom_mark_shape.md).
  Default is `margin(2, 2, 2, 2, "pt")`.

- label.width:

  Soft target width for wrapping labels, passed to
  [`geom_mark_shape()`](https://alserglab.github.io/mascarade/reference/geom_mark_shape.md).
  A grid unit (e.g. `unit(30, "mm")`); labels are balanced across lines
  to keep line widths even and close to this width, without a short
  dangling line. `NULL` (default) leaves labels unwrapped. See
  [`geom_mark_shape()`](https://alserglab.github.io/mascarade/reference/geom_mark_shape.md)
  for details.

- simp_ratio:

  Fraction of the polygon bounding-box area used to simplify cluster
  polygons before label placement: small inward (concave) vertices whose
  cut-off area is below `simp_ratio * bbox_area` are removed, which
  speeds up the box-fit and leader-routing geometry (both scale with
  vertex count). The simplified polygon encloses the original, so labels
  never overlap the real cluster. Larger values simplify more; set to
  `0` to disable. Default `0.001`.

- con.type:

  Leader / label-mark style passed to
  [`geom_mark_shape()`](https://alserglab.github.io/mascarade/reference/geom_mark_shape.md):
  one of `"ledge"`, `"line"`, `"box"`, or `"none"` (see the
  [`geom_mark_shape()`](https://alserglab.github.io/mascarade/reference/geom_mark_shape.md)
  Details). Default `"ledge"`.

## Value

A list of ggplot2 components suitable for adding to a plot with `+`,
containing a
[`ggplot2::coord_cartesian()`](https://ggplot2.tidyverse.org/reference/coord_cartesian.html)
specification and a
[`geom_mark_shape()`](https://alserglab.github.io/mascarade/reference/geom_mark_shape.md)
layer.

## Details

The first two columns of `maskTable` are used as x/y coordinates.
Cluster labels are taken from `maskTable$cluster`. Shapes are grouped by
`maskTable$group`.

## See also

- [`geom_mark_shape()`](https://alserglab.github.io/mascarade/reference/geom_mark_shape.md)

## Examples

``` r
data("exampleMascarade")
maskTable <- generateMask(dims=exampleMascarade$dims,
                          clusters=exampleMascarade$clusters)
library(ggplot2)
basePlot <- ggplot(do.call(cbind, exampleMascarade)) +
    geom_point(aes(x=UMAP_1, y=UMAP_2, color=GNLY)) +
    scale_color_gradient2(low = "#404040", high="red") +
    theme_classic()

basePlot + fancyMask(maskTable, ratio=1, cols=scales::hue_pal())

```

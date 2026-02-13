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
  cols = "inherit",
  label = TRUE,
  label.fontsize = 10,
  label.buffer = unit(0, "cm"),
  label.fontface = "plain",
  label.margin = margin(2, 2, 2, 2, "pt")
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

  - `"inherit"` (default) — inherits colors from the discrete color
    scale of the plot that `fancyMask()` is added to (e.g., from
    [`scale_color_manual()`](https://ggplot2.tidyverse.org/reference/scale_manual.html)).
    Falls back to black if no discrete color scale is found.

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

  Boolean flag wheter the labels should be displayed.

- label.fontsize:

  Label font size passed to
  [`geom_mark_shape()`](https://alserglab.github.io/mascarade/reference/geom_mark_shape.md).
  Default is `10`.

- label.buffer:

  Label buffer distance passed to
  [`geom_mark_shape()`](https://alserglab.github.io/mascarade/reference/geom_mark_shape.md).
  Default is `unit(0, "cm")`.

- label.fontface:

  Label font face passed to
  [`geom_mark_shape()`](https://alserglab.github.io/mascarade/reference/geom_mark_shape.md).
  Default is `"plain"`.

- label.margin:

  Label margin passed to
  [`geom_mark_shape()`](https://alserglab.github.io/mascarade/reference/geom_mark_shape.md).
  Default is `margin(2, 2, 2, 2, "pt")`.

## Value

A list of ggplot2 components suitable for adding to a plot with `+`,
containing a
[`ggplot2::coord_cartesian()`](https://ggplot2.tidyverse.org/reference/coord_cartesian.html)
specification and a
[`geom_mark_shape()`](https://alserglab.github.io/mascarade/reference/geom_mark_shape.md)
layer. When `cols = "inherit"`, returns an opaque object whose colors
are resolved when added to a plot.

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

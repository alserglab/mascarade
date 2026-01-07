# Annotate areas with polygonal shapes

This geom lets you annotate sets of points via polygonal shapes. Unlike
other `ggforce::geom_mark_*` functions, `geom_mark_shape` should be
explicitly provided with the shape coordinates. As in
[`ggforce::geom_shape`](https://ggforce.data-imaginist.com/reference/geom_shape.html),
the polygon can be expanded/contracted and corners can be rounded, which
is controlled by `expand` and `radius` parameters.

## Usage

``` r
geom_mark_shape(
  mapping = NULL,
  data = NULL,
  stat = "identity",
  position = "identity",
  expand = 0,
  radius = 0,
  label.margin = margin(2, 2, 2, 2, "mm"),
  label.width = NULL,
  label.minwidth = unit(50, "mm"),
  label.hjust = 0,
  label.fontsize = 12,
  label.family = "",
  label.lineheight = 1,
  label.fontface = c("bold", "plain"),
  label.fill = "white",
  label.colour = "black",
  label.buffer = unit(10, "mm"),
  con.colour = "black",
  con.size = 0.5,
  con.type = "elbow",
  con.linetype = 1,
  con.border = "one",
  con.cap = unit(3, "mm"),
  con.arrow = NULL,
  ...,
  na.rm = FALSE,
  show.legend = NA,
  inherit.aes = TRUE
)
```

## Arguments

- mapping:

  Set of aesthetic mappings created by
  [`aes()`](https://ggplot2.tidyverse.org/reference/aes.html). If
  specified and `inherit.aes = TRUE` (the default), it is combined with
  the default mapping at the top level of the plot. You must supply
  `mapping` if there is no plot mapping.

- data:

  The data to be displayed in this layer. There are three options:

  If `NULL`, the default, the data is inherited from the plot data as
  specified in the call to
  [`ggplot()`](https://ggplot2.tidyverse.org/reference/ggplot.html).

  A `data.frame`, or other object, will override the plot data. All
  objects will be fortified to produce a data frame. See
  [`fortify()`](https://ggplot2.tidyverse.org/reference/fortify.html)
  for which variables will be created.

  A `function` will be called with a single argument, the plot data. The
  return value must be a `data.frame`, and will be used as the layer
  data. A `function` can be created from a `formula` (e.g.
  `~ head(.x, 10)`).

- stat:

  The statistical transformation to use on the data for this layer. When
  using a `geom_*()` function to construct a layer, the `stat` argument
  can be used the override the default coupling between geoms and stats.
  The `stat` argument accepts the following:

  - A `Stat` ggproto subclass, for example `StatCount`.

  - A string naming the stat. To give the stat as a string, strip the
    function name of the `stat_` prefix. For example, to use
    [`stat_count()`](https://ggplot2.tidyverse.org/reference/geom_bar.html),
    give the stat as `"count"`.

  - For more information and other ways to specify the stat, see the
    [layer
    stat](https://ggplot2.tidyverse.org/reference/layer_stats.html)
    documentation.

- position:

  A position adjustment to use on the data for this layer. This can be
  used in various ways, including to prevent overplotting and improving
  the display. The `position` argument accepts the following:

  - The result of calling a position function, such as
    [`position_jitter()`](https://ggplot2.tidyverse.org/reference/position_jitter.html).
    This method allows for passing extra arguments to the position.

  - A string naming the position adjustment. To give the position as a
    string, strip the function name of the `position_` prefix. For
    example, to use
    [`position_jitter()`](https://ggplot2.tidyverse.org/reference/position_jitter.html),
    give the position as `"jitter"`.

  - For more information and other ways to specify the position, see the
    [layer
    position](https://ggplot2.tidyverse.org/reference/layer_positions.html)
    documentation.

- expand:

  A numeric or unit vector of length one, specifying the expansion
  amount. Negative values will result in contraction instead. If the
  value is given as a numeric it will be understood as a proportion of
  the plot area width.

- radius:

  As `expand` but specifying the corner radius.

- label.margin:

  The margin around the annotation boxes, given by a call to
  [`ggplot2::margin()`](https://ggplot2.tidyverse.org/reference/element.html).

- label.width:

  A fixed width for the label. Set to `NULL` to let the text or
  `label.minwidth` decide.

- label.minwidth:

  The minimum width to provide for the description. If the size of the
  label exceeds this, the description is allowed to fill as much as the
  label.

- label.hjust:

  The horizontal justification for the annotation. If it contains two
  elements the first will be used for the label and the second for the
  description.

- label.fontsize:

  The size of the text for the annotation. If it contains two elements
  the first will be used for the label and the second for the
  description.

- label.family:

  The font family used for the annotation. If it contains two elements
  the first will be used for the label and the second for the
  description.

- label.lineheight:

  The height of a line as a multipler of the fontsize. If it contains
  two elements the first will be used for the label and the second for
  the description.

- label.fontface:

  The font face used for the annotation. If it contains two elements the
  first will be used for the label and the second for the description.

- label.fill:

  The fill colour for the annotation box. Use `"inherit"` to use the
  fill from the enclosure or `"inherit_col"` to use the border colour of
  the enclosure.

- label.colour:

  The text colour for the annotation. If it contains two elements the
  first will be used for the label and the second for the description.
  Use `"inherit"` to use the border colour of the enclosure or
  `"inherit_fill"` to use the fill colour from the enclosure.

- label.buffer:

  The size of the region around the mark where labels cannot be placed.

- con.colour:

  The colour for the line connecting the annotation to the mark. Use
  `"inherit"` to use the border colour of the enclosure or
  `"inherit_fill"` to use the fill colour from the enclosure.

- con.size:

  The width of the connector. Use `"inherit"` to use the border width of
  the enclosure.

- con.type:

  The type of the connector. Either `"elbow"`, `"straight"`, or
  `"none"`.

- con.linetype:

  The linetype of the connector. Use `"inherit"` to use the border
  linetype of the enclosure.

- con.border:

  The bordertype of the connector. Either `"one"` (to draw a line on the
  horizontal side closest to the mark), `"all"` (to draw a border on all
  sides), or `"none"` (not going to explain that one).

- con.cap:

  The distance before the mark that the line should stop at.

- con.arrow:

  An arrow specification for the connection using
  [`grid::arrow()`](https://rdrr.io/r/grid/arrow.html) for the end
  pointing towards the mark.

- ...:

  Other arguments passed on to
  [`layer()`](https://ggplot2.tidyverse.org/reference/layer.html)'s
  `params` argument. These arguments broadly fall into one of 4
  categories below. Notably, further arguments to the `position`
  argument, or aesthetics that are required can *not* be passed through
  `...`. Unknown arguments that are not part of the 4 categories below
  are ignored.

  - Static aesthetics that are not mapped to a scale, but are at a fixed
    value and apply to the layer as a whole. For example,
    `colour = "red"` or `linewidth = 3`. The geom's documentation has an
    **Aesthetics** section that lists the available options. The
    'required' aesthetics cannot be passed on to the `params`. Please
    note that while passing unmapped aesthetics as vectors is
    technically possible, the order and required length is not
    guaranteed to be parallel to the input data.

  - When constructing a layer using a `stat_*()` function, the `...`
    argument can be used to pass on parameters to the `geom` part of the
    layer. An example of this is
    `stat_density(geom = "area", outline.type = "both")`. The geom's
    documentation lists which parameters it can accept.

  - Inversely, when constructing a layer using a `geom_*()` function,
    the `...` argument can be used to pass on parameters to the `stat`
    part of the layer. An example of this is
    `geom_area(stat = "density", adjust = 0.5)`. The stat's
    documentation lists which parameters it can accept.

  - The `key_glyph` argument of
    [`layer()`](https://ggplot2.tidyverse.org/reference/layer.html) may
    also be passed on through `...`. This can be one of the functions
    described as [key
    glyphs](https://ggplot2.tidyverse.org/reference/draw_key.html), to
    change the display of the layer in the legend.

- na.rm:

  If `FALSE`, the default, missing values are removed with a warning. If
  `TRUE`, missing values are silently removed.

- show.legend:

  logical. Should this layer be included in the legends? `NA`, the
  default, includes if any aesthetics are mapped. `FALSE` never
  includes, and `TRUE` always includes. It can also be a named logical
  vector to finely select the aesthetics to display.

- inherit.aes:

  If `FALSE`, overrides the default aesthetics, rather than combining
  with them. This is most useful for helper functions that define both
  data and aesthetics and shouldn't inherit behaviour from the default
  plot specification, e.g.
  [`borders()`](https://ggplot2.tidyverse.org/reference/annotation_borders.html).

## Value

A ggplot2 layer
([`ggplot2::layer`](https://ggplot2.tidyverse.org/reference/layer.html))
that adds polygonal shape annotations to a plot.

## Aesthetics

`geom_mark_shape` understand the following aesthetics (required
aesthetics are in bold):

- **x**

- **y**

- x0 *(used to anchor the label)*

- y0 *(used to anchor the label)*

- filter

- label

- description

- color

- fill

- group

- size

- linetype

- alpha

## Annotation

All `geom_mark_*` allow you to put descriptive textboxes connected to
the mark on the plot, using the `label` and `description` aesthetics.
The textboxes are automatically placed close to the mark, but without
obscuring any of the datapoints in the layer. The placement is dynamic
so if you resize the plot you'll see that the annotation might move
around as areas become big enough or too small to fit the annotation. If
there's not enough space for the annotation without overlapping data it
will not get drawn. In these cases try resizing the plot, change the
size of the annotation, or decrease the buffer region around the marks.

## Filtering

Often marks are used to draw attention to, or annotate specific features
of the plot and it is thus not desirable to have marks around
everything. While it is possible to simply pre-filter the data used for
the mark layer, the `geom_mark_*` geoms also comes with a dedicated
`filter` aesthetic that, if set, will remove all rows where it
evalutates to `FALSE`. There are multiple benefits of using this instead
of prefiltering. First, you don't have to change your data source,
making your code more adaptable for exploration. Second, the data
removed by the filter aesthetic is remembered by the geom, and any
annotation will take care not to overlap with the removed data.

## Examples

``` r
library(ggplot2)
shape1 <- data.frame(
    x = c(0, 3, 3, 2, 2, 1, 1, 0),
    y = c(0, 0, 3, 3, 1, 1, 3, 3),
label="bracket"
)
shape2 <- data.frame(
    x = c(0, 3, 3, 0)+4,
    y = c(0, 0, 3, 3),
    label="square"
)
shape3 <- data.frame(
    x = c(0, 1.5, 3, 1.5)+8,
    y = c(1.5, 0, 1.5, 3),
    label="diamond"
)

ggplot(rbind(shape1, shape2, shape3), aes(x=x, y=y, label=label, color=label, fill=label)) +
    geom_mark_shape() +
    ylim(0, 5)


```

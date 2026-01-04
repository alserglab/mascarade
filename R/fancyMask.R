#' Generate ggplot2 layers for a labeled cluster mask
#'
#' Convenience helper that returns a list of ggplot2 components
#' that draws polygon-like outlines and
#' places cluster labels.
#' The plotting limits are expanded (via `limits.expand`) to provide
#' extra room for labels.
#'
#' @param maskTable A data.frame of mask coordinates. The first two
#'   columns are interpreted as x/y coordinates (in that order). Must contain
#'   at least the columns `cluster` (a factor) and `group` (grouping identifier
#'   passed to `geom_mark_shape()`).
#' @param ratio Optional aspect ratio passed to `ggplot2::coord_cartesian()`.
#'   Use `1` for equal scaling. Default is `NULL` (no fixed ratio).
#' @param limits.expand Numeric scalar giving the fraction of the x/y range to
#'   expand on both sides when setting plot limits. Default is `0.1` with labels and 0.05 with no labels.
#' @param linewidth Line width passed to `geom_mark_shape()` for the
#'   outline. Default is `1`.
#' @param shape.expand Expansion or contraction applied to the marked shapes,
#'   passed to `geom_mark_shape(expand = ...)`. Default is
#'   `unit(-linewidth, "pt")`.
#' @param label Boolean flag wheter the labels should be displayed.
#' @param label.fontsize Label font size passed to `geom_mark_shape()`.
#'   Default is `10`.
#' @param label.buffer Label buffer distance passed to
#'   `geom_mark_shape()`. Default is `unit(0, "cm")`.
#' @param label.fontface Label font face passed to
#'   `geom_mark_shape()`. Default is `"plain"`.
#' @param label.margin Label margin passed to
#'   `geom_mark_shape()`. Default is `margin(2, 2, 2, 2, "pt")`.
#'
#' @return A list of ggplot2 components suitable for adding to a plot with `+`,
#'   containing:
#'   1. a `ggplot2::coord_cartesian()` specification, and
#'   2. a `geom_mark_shape()` layer.
#'
#' @details
#' The first two columns of `maskTable` are used as x/y coordinates. Cluster
#' labels are taken from `maskTable$cluster`. Shapes are grouped by
#' `maskTable$group`.
#'
#' @seealso
#' * `geom_mark_shape()`
#' @examples
#' data("exampleMascarade")
#' maskTable <- generateMask(dims=exampleMascarade$dims,
#'                           clusters=exampleMascarade$clusters)
#' library(ggplot2)
#' ggplot(do.call(cbind, exampleMascarade)) +
#'     geom_point(aes(x=UMAP_1, y=UMAP_2, color=GNLY)) +
#'     fancyMask(maskTable, ratio=1) +
#'     theme_classic()
#'
#' @export
#' @importFrom rlang .data
#' @importFrom stats setNames
#' @importFrom ggplot2 aes coord_cartesian
fancyMask <- function(maskTable,
                      ratio=NULL,
                      limits.expand = ifelse(label, 0.1, 0.05),
                      linewidth=1,
                      shape.expand=linewidth*unit(-1, "pt"),
                      label=TRUE,
                      label.fontsize = 10,
                      label.buffer = unit(0, "cm"),
                      label.fontface = "plain",
                      label.margin = margin(2, 2, 2, 2, "pt")
                      ) {
    xvar <- colnames(maskTable)[1]
    yvar <- colnames(maskTable)[2]


    # expanding to give a bit more space for labels
    xyRanges <- apply(maskTable[, 1:2], 2, range)
    xyWidths <- apply(xyRanges, 2, diff)
    xyRanges <- xyRanges + c(-1, 1)  %*% t(xyWidths * limits.expand)

    # TODO: make colors controllable
    colors <- NULL
    if (is.null(colors)) {
        clusterLevels <- levels(maskTable$cluster)
        pal <- setNames(scales::hue_pal()(length(clusterLevels)), clusterLevels)
        colors <- pal[maskTable$cluster]
    }

    xlim <- range(maskTable[[xvar]])
    ylim <- range(maskTable[[yvar]])

    if (label) {
        shapes <- geom_mark_shape(data=maskTable,
                                 fill = NA,
                                 x=maskTable[[xvar]],
                                 y=maskTable[[yvar]],
                                 aes(group=group,
                                     label=cluster),
                                 colour=colors,
                                 linewidth=linewidth,
                                 expand=shape.expand,
                                 label.fontsize = label.fontsize,
                                 label.buffer = label.buffer,
                                 label.fontface = label.fontface,
                                 label.margin = label.margin,
                                 label.minwidth = 0,
                                 label.lineheight = 0,
                                 con.cap=0,
                                 con.type = "straight",
                                 con.colour = "inherit")
    } else {
        shapes <- geom_shape(data=maskTable,
                             fill = NA,
                             x=maskTable[[xvar]],
                             y=maskTable[[yvar]],
                             aes(group=group,
                                 label=cluster),
                             colour=colors,
                             linewidth=linewidth,
                             expand=shape.expand)
    }

    res <- list(
        coord_cartesian(xlim=xyRanges[,1],
                        ylim=xyRanges[,2],
                        ratio=ratio,
                        expand=FALSE), # already expanded
        shapes
    )
}

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
#' @param label Boolean flag whether the labels should be displayed.
#' @param label.largest Boolean flag. When `TRUE` (default), only the largest
#'   part of each cluster is labelled; smaller disconnected parts are drawn but
#'   not labelled. When `FALSE`, all parts are labelled. Ignored when
#'   `label = FALSE`.
#' @param label.fontsize Label font size passed to `geom_mark_shape()`.
#'   Default is `10`.
#' @param label.buffer Label buffer distance passed to
#'   `geom_mark_shape()`. Default is `unit(0, "cm")`.
#' @param label.fontface Label font face passed to
#'   `geom_mark_shape()`. Default is `"plain"`.
#' @param cols Color specification for cluster outlines (and labels). One of:
#'
#'   * `"auto"` (default) — inspects the plot at the time `fancyMask()` is
#'     added with `+`. If a layer maps `colour` to a discrete (non-numeric)
#'     variable, the mask joins that scale via `aes(colour = cluster)` so
#'     colours stay in sync regardless of `scale_color_*()` order. Otherwise
#'     (continuous colour, constant colour, or no colour aesthetic) explicit
#'     colours from `scales::hue_pal()` are baked in and the plot's scale
#'     system is left untouched.
#'   * `"inherit"` — always maps `colour` as an aesthetic (`aes(colour =
#'     cluster)`), unconditionally joining whatever colour scale is present.
#'     Useful when you want to force scale sharing; will error if the existing
#'     scale is continuous.
#'   * A palette function that accepts a single integer `n` and returns `n`
#'     colors (e.g., `scales::hue_pal()`, `rainbow`).
#'   * A single color string — applied to every cluster.
#'   * An unnamed character vector of length equal to the number of clusters —
#'     colors are assigned to clusters in factor-level order.
#'   * A named character vector — names must match cluster levels; order does
#'     not matter.
#' @param label.margin Label margin passed to
#'   `geom_mark_shape()`. Default is `margin(2, 2, 2, 2, "pt")`.
#' @param simp_ratio Fraction of the polygon bounding box area used as the
#'   label-placement simplification threshold. Cluster polygons are simplified
#'   before the label placement search by removing small concave vertices,
#'   which reduces computation while guaranteeing the simplified polygon
#'   encloses the original. Larger values simplify more aggressively;
#'   set to `0` to disable. Default is `0.001`.
#'
#' @return A list of ggplot2 components suitable for adding to a plot with `+`,
#'   containing a `ggplot2::coord_cartesian()` specification and a
#'   `geom_mark_shape()` layer.
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
#' basePlot <- ggplot(do.call(cbind, exampleMascarade)) +
#'     geom_point(aes(x=UMAP_1, y=UMAP_2, color=GNLY)) +
#'     scale_color_gradient2(low = "#404040", high="red") +
#'     theme_classic()
#'
#' basePlot + fancyMask(maskTable, ratio=1, cols=scales::hue_pal())
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
                      cols="auto",
                      label=TRUE,
                      label.largest=TRUE,
                      label.fontsize = 10,
                      label.buffer = unit(0, "cm"),
                      label.fontface = "plain",
                      label.margin = margin(2, 2, 2, 2, "pt"),
                      simp_ratio = 0.001
                      ) {

    if (identical(cols, "auto")) {
        # Defer: colour strategy is decided in ggplot_add.fancyMask once the
        # plot context (other layers and their aesthetics) is known.
        structure(
            list(maskTable    = maskTable,
                 ratio        = ratio,
                 limits.expand = limits.expand,
                 linewidth    = linewidth,
                 shape.expand = shape.expand,
                 cols         = cols,
                 label        = label,
                 label.largest = label.largest,
                 label.fontsize = label.fontsize,
                 label.buffer = label.buffer,
                 label.fontface = label.fontface,
                 label.margin = label.margin,
                 simp_ratio   = simp_ratio),
            class = "fancyMask"
        )
    } else {
        buildFancyMaskLayers(maskTable = maskTable,
                             ratio = ratio,
                             limits.expand = limits.expand,
                             linewidth = linewidth,
                             shape.expand = shape.expand,
                             cols = cols,
                             label = label,
                             label.largest = label.largest,
                             label.fontsize = label.fontsize,
                             label.buffer = label.buffer,
                             label.fontface = label.fontface,
                             label.margin = label.margin,
                             simp_ratio = simp_ratio)
    }
}

# Returns TRUE if any layer (or the global plot mapping) maps 'colour' to a
# non-numeric (discrete) variable. Used by ggplot_add.fancyMask to decide
# whether to join the existing scale or bake in explicit colours.
hasDiscreteColour <- function(plot) {
    checkMapping <- function(mapping, data) {
        col_q <- mapping[["colour"]]
        if (is.null(col_q) || is.null(data) || inherits(data, "waiver")) {
            return(FALSE)
        }
        tryCatch({
            vals <- rlang::eval_tidy(col_q, data = as.data.frame(data))
            !is.numeric(vals)
        }, error = function(e) FALSE)
    }

    if (checkMapping(plot$mapping, plot$data)) return(TRUE)

    for (layer in plot$layers) {
        layer_data <- layer$data
        if (inherits(layer_data, "waiver")) layer_data <- plot$data
        if (is.function(layer_data)) next
        if (checkMapping(layer$mapping, layer_data)) return(TRUE)
    }

    FALSE
}

#' @export
#' @importFrom ggplot2 ggplot_add
#' @importFrom rlang eval_tidy
ggplot_add.fancyMask <- function(object, plot, ...) {
    cols <- if (hasDiscreteColour(plot)) "inherit" else scales::hue_pal()
    layers <- buildFancyMaskLayers(
        maskTable     = object$maskTable,
        ratio         = object$ratio,
        limits.expand = object$limits.expand,
        linewidth     = object$linewidth,
        shape.expand  = object$shape.expand,
        cols          = cols,
        label         = object$label,
        label.largest = object$label.largest,
        label.fontsize = object$label.fontsize,
        label.buffer  = object$label.buffer,
        label.fontface = object$label.fontface,
        label.margin  = object$label.margin,
        simp_ratio    = object$simp_ratio
    )
    ggplot2::ggplot_add(layers, plot, ...)
}

getClusterLevels <- function(x) {
    if (is.factor(x)) levels(x) else unique(x)
}

resolveCols <- function(cols, clusterLevels) {
    nClusters <- length(clusterLevels)
    if (is.function(cols)) {
        setNames(cols(nClusters), clusterLevels)
    } else if (length(cols) == 1L) {
        setNames(rep(cols, nClusters), clusterLevels)
    } else if (is.null(names(cols))) {
        if (length(cols) != nClusters) {
            stop("Length of unnamed `cols` (", length(cols),
                 ") must equal the number of clusters (", nClusters, ")")
        }
        setNames(cols, clusterLevels)
    } else {
        missing <- setdiff(clusterLevels, names(cols))
        if (length(missing) > 0L) {
            stop("Named `cols` is missing entries for cluster(s): ",
                 paste(missing, collapse = ", "))
        }
        cols[clusterLevels]
    }
}

buildFancyMaskLayers <- function(maskTable, ratio, limits.expand, linewidth,
                                 shape.expand, cols, label, label.largest,
                                 label.fontsize, label.buffer, label.fontface,
                                 label.margin, simp_ratio = 0.001) {
    xvar <- colnames(maskTable)[1]
    yvar <- colnames(maskTable)[2]

    # expanding to give a bit more space for labels
    xyRanges <- apply(maskTable[, 1:2], 2, range)
    xyWidths <- apply(xyRanges, 2, diff)
    xyRanges <- xyRanges + c(-1, 1)  %*% t(xyWidths * limits.expand)

    if (label) {
        # When label.largest=TRUE, only label the first part per cluster
        # (generateMask guarantees part #1 is the largest by polygon area).
        if (label.largest) {
            isLargest <- grepl("#1$", maskTable$part)
            labelCol <- ifelse(isLargest, as.character(maskTable$cluster), NA_character_)
        } else {
            labelCol <- as.character(maskTable$cluster)
        }
        maskTable <- cbind(maskTable, .label_display = labelCol)
    }

    if (identical(cols, "inherit")) {
        # Colour is expressed as an aesthetic so ggplot2 resolves it at build
        # time — after all scale_color_*() calls have been applied, regardless
        # of the order they were added to the plot.
        if (label) {
            shapes <- geom_mark_shape(data = maskTable,
                                     fill = NA,
                                     x = maskTable[[xvar]],
                                     y = maskTable[[yvar]],
                                     aes(group = group,
                                         label = .data$.label_display,
                                         colour = cluster),
                                     linewidth = linewidth,
                                     expand = shape.expand,
                                     show.legend = FALSE,
                                     label.fontsize = label.fontsize,
                                     label.buffer = label.buffer,
                                     label.fontface = label.fontface,
                                     label.margin = label.margin,
                                     simp_ratio = simp_ratio,
                                     label.minwidth = 0,
                                     label.lineheight = 0,
                                     con.cap = 0,
                                     con.type = "straight",
                                     con.colour = "inherit")
        } else {
            shapes <- geom_shape(data = maskTable,
                                 fill = NA,
                                 x = maskTable[[xvar]],
                                 y = maskTable[[yvar]],
                                 aes(group = group,
                                     colour = cluster),
                                 linewidth = linewidth,
                                 expand = shape.expand,
                                 show.legend = FALSE)
        }

        return(list(
            coord_cartesian(xlim = xyRanges[, 1],
                            ylim = xyRanges[, 2],
                            ratio = ratio,
                            expand = FALSE),
            shapes
        ))
    }

    # processing explicit color options
    clusterLevels <- getClusterLevels(maskTable$cluster)
    pal <- resolveCols(cols, clusterLevels)
    # Use as.character() to force name-based lookup regardless of whether
    # maskTable$cluster is a factor or character vector.
    colors <- unname(pal[as.character(maskTable$cluster)])

    if (label) {
        shapes <- geom_mark_shape(data=maskTable,
                                 fill = NA,
                                 x=maskTable[[xvar]],
                                 y=maskTable[[yvar]],
                                 aes(group=group,
                                     label=.data$.label_display),
                                 colour=colors,
                                 linewidth=linewidth,
                                 expand=shape.expand,
                                 label.fontsize = label.fontsize,
                                 label.buffer = label.buffer,
                                 label.fontface = label.fontface,
                                 label.margin = label.margin,
                                 simp_ratio = simp_ratio,
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
                             aes(group=group),
                             colour=colors,
                             linewidth=linewidth,
                             expand=shape.expand)
    }

    list(
        coord_cartesian(xlim=xyRanges[,1],
                        ylim=xyRanges[,2],
                        ratio=ratio,
                        expand=FALSE), # already expanded
        shapes
    )
}

expandedRange2d <- function(x, y, fraction=0.05, fixAspectRatio=TRUE) {
    xRange <- range(x)
    xWidth <- (xRange[2] - xRange[1]) * (1 + fraction)

    yRange <- range(y)
    yWidth <- (yRange[2] - yRange[1]) * (1 + fraction)

    if (fixAspectRatio) {
        xWidth <- yWidth <- max(xWidth, yWidth)
    }

    xCenter <- mean(xRange)
    yCenter <- mean(yRange)

    return(
        c(xCenter - xWidth/2, xCenter + xWidth/2,
          yCenter - yWidth/2, yCenter + yWidth/2))
}

#' @importFrom spatstat.geom tiles tess connected as.polygonal
#' @importFrom  data.table rbindlist
borderTableFromMask <- function(curMask) {
    parts <- tiles(tess(image=connected(curMask > 0)))

    curBorderTable <- list()

    for (partIdx in seq_along(parts)) {
        part <- parts[[partIdx]]

        partBoundary <- as.polygonal(part)
        lines <- partBoundary$bdry

        curBorderTable <- c(curBorderTable, lapply(seq_along(lines), function(lineIdx) {
            curLine <- lines[[lineIdx]]
            xs <- curLine$x
            ys <- curLine$y

            # make lines closed
            xs <- c(xs, xs[1])
            ys <- c(ys, ys[1])

            # remove steps
            xs <- (head(xs, -1) + tail(xs, -1)) / 2
            ys <- (head(ys, -1) + tail(ys, -1)) / 2

            # make lines closed again
            xs <- c(xs, xs[1])
            ys <- c(ys, ys[1])

            res <- data.table(x=xs, y=ys)
            res[, part := partIdx]
            res[, group := lineIdx]
            res[]
        }))
    }
    rbindlist(curBorderTable)
}

#' Generate mask for clusters on 2D dimensional reduction plots
#'
#' Internally the function rasterizes and smoothes the density plots.
#' @param dims matrix of point coordinates.
#'      Rows are points, columns are dimensions. Only the first two columns are used.
#' @param clusters vector of cluster annotations.
#'      Should be the same length as the number of rows in `dims`.
#' @param gridSize width and height of the raster used internally
#' @param minDensity minimal required density for the grid cells to be included in the mask.
#'      Decreasing this parameter will expand masks.
#' @param smoothSigma sigma used in Gaussian smoothing represented as a fraction of plot width.
#'      Increasing this parameter can help dealing with sparse regions.
#' @param type controls the behavior of the method.
#'      When set to "partition" (default) generated masks are mutually exclusive.
#'      When set to "independent" masks can overlap.
#' @returns data.table with points representing the mask borders.
#'      Each individual border line corresponds to a single level of `group` column.
#'      Cluster assignment is in `cluster` column.
#' @importFrom data.table rbindlist data.table setnames
#' @importFrom utils head tail
#' @export
#' @examples
#' data("exampleMascarade")
#' res <- generateMask(dims=exampleMascarade$dims,
#'                     clusters=exampleMascarade$clusters)
#' \dontrun{
#' data <- data.table(exampleMascarade$dims,
#'                    cluster=exampleMascarade$clusters,
#'                    exampleMascarade$features)
#' ggplot(data, aes(x=UMAP_1, y=UMAP_2)) +
#'     geom_point(aes(color=cluster)) +
#'     geom_path(data=maskTable, aes(group=group)) +
#'     coord_fixed() +
#'     theme_classic()
#' }
generateMask <- function(dims, clusters,
                         gridSize=200,
                         minDensity=0.1,
                         smoothSigma=0.025,
                         type=c("partition", "independent")) {
    type <- match.arg(type)

    clusterLevels <- unique(clusters)

    gridRange <- expandedRange2d(dims[, 1], dims[ ,2])
    windowWidth <- gridRange[2] - gridRange[1]

    smoothSigma <- smoothSigma * windowWidth

    window <- spatstat.geom::as.mask(spatstat.geom::owin(xrange = gridRange[1:2],
                                                         yrange = gridRange[3:4]),
                                     dimyx=gridSize)

    dims <- dims[, 1:2]

    if (is.null(colnames(dims))) {
        colnames(dims) <- c("x", "y")
    }

    points <- spatstat.geom::ppp(dims[, 1], dims[, 2], window=window)


    allDensities <- lapply(clusterLevels, function(cluster) {
        res <- spatstat.geom::pixellate(points[clusters == cluster], dimyx=gridSize)
    })

    allDensitiesSmoothed <- lapply(allDensities, spatstat.explore::blur, sigma = smoothSigma)

    densityThresholds <- pmin(vapply(allDensitiesSmoothed, max, numeric(1)) / 2,
                              minDensity)

    if (type == "partition") {
        # backgroundDensity <- spatstat.geom::as.im(window) * minDensity

        allDensitiesMax <- spatstat.geom::im.apply(allDensitiesSmoothed, which.max)
        allDensitiesSmoothed <- lapply(seq_along(clusterLevels), function(i) {
            allDensitiesSmoothed[[i]] * (allDensitiesMax == i)
        })
    }

    borderTable <- rbindlist(lapply(seq_along(clusterLevels), function(i) {
        curMask <- (allDensitiesSmoothed[[i]] >= densityThresholds[i]) * allDensitiesSmoothed[[i]]
        if (sum(curMask) == 0) {
            warning(sprintf("Mask is empty for cluster %s", clusterLevels[i]))
            return(NULL)
        }
        curTable <- borderTableFromMask(curMask)
        curTable[, cluster := clusterLevels[i]]
        curTable[, part := paste0(cluster, "#", part)]
        curTable[, group := paste0(part, "#", group)]
        curTable[]
    }))

    setnames(borderTable, c("x", "y"), colnames(dims))

    return(borderTable)
}

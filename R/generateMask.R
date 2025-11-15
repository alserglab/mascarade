makeGridWindow <- function(dims, gridSize, fraction=0.05) {
    xyRanges <- apply(dims, 2, range)

    xyWidths <- (xyRanges[2,] - xyRanges[1,]) * (1 + fraction)

    xyCenters <- colMeans(xyRanges)

    gridStep <- sqrt(prod(xyWidths))/gridSize

    # switch yx and xy
    xyResolution <- ceiling(xyWidths/gridStep)

    xyWidths <- gridStep*xyResolution

    window <- spatstat.geom::as.mask(
        spatstat.geom::owin(xrange = c(xyCenters[1]-xyWidths[1]/2, xyCenters[1]+xyWidths[1]/2),
                            yrange = c(xyCenters[2]-xyWidths[2]/2, xyCenters[2]+xyWidths[2]/2)),
                                     dimyx=rev(xyResolution))
}

#' @importFrom spatstat.geom tiles tess connected as.polygonal
#' @importFrom  data.table rbindlist
#' @keywords internal
borderTableFromMask <- function(curMask, curDensity, keepMax=TRUE) {
    parts <- tiles(tess(image=connected(curMask)))

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

# robust to empty masks
splitWhichMaxLevels <- function(whichMaxDensity, nLevels) {
    lapply(seq_len(nLevels), function(i) {
        res <- (whichMaxDensity == i)
        res[res == 0] <- NA # so that as.owin works as expected
        res <- spatstat.geom::as.owin(res)
    })
}

# TODO window argument shouldn't be really needed
removeMaskIntersections <- function(curMasks, window) {
    maskWeights <- lapply(seq_along(curMasks), function(i) {
        distmap(complement.owin(curMasks[[i]]))
    })

    backgroundDensityOne <- spatstat.geom::as.im(window)

    # TODO: maybe consider cell densities here somehow
    whichMaxDensity <- spatstat.geom::im.apply(
        c(list(backgroundDensityOne*0.01), maskWeights), which.max) - 1


    curMasks <- splitWhichMaxLevels(whichMaxDensity, nLevels=length(curMasks))
    curMasks
}

getConnectedParts <- function(curMask, curDensity, minSize, absolutelyMinSize=5) {
    parts <- tiles(tess(image=connected(curMask)))
    partSizes <- vapply(parts, function(part) {
        sum(as.matrix(part) * as.matrix(curDensity))
    }, FUN.VALUE = numeric(1))

    parts <- parts[partSizes >= min(minSize, max(c(partSizes, absolutelyMinSize)))]
    unname(parts)
}


#' Generate mask for clusters on 2D dimensional reduction plots
#'
#' Internally the function rasterizes and smoothes the density plots.
#' @param dims matrix of point coordinates.
#'      Rows are points, columns are dimensions. Only the first two columns are used.
#' @param clusters vector of cluster annotations.
#'      Should be the same length as the number of rows in `dims`.
#' @param gridSize target width and height of the raster used internally
#' @param expand distance used to expand borders, represented as a fraction of sqrt(width*height). Default: 1/200.
#' @param minDensity Deprecated. Doesn't do anything.
#' @param smoothSigma Deprecated. Parameter controlling smoothing and joining close cells into groups, represented as a fraction of sqrt(width*height).
#'      Increasing this parameter can help dealing with sparse regions.
#' @param minSize Groups of less than `minSize` points are ignored, unless it is the only group for a cluster
#' @param kernel Deprecated. Doesn't do anything.
#' @param type Deprecated. Doesn't do anything.

#' @returns data.table with points representing the mask borders.
#'      Each individual border line corresponds to a single level of `group` column.
#'      Cluster assignment is in `cluster` column.
#' @importFrom data.table rbindlist data.table setnames
#' @importFrom utils head tail
#' @import spatstat.geom spatstat.explore
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
                         expand=0.005,
                         minDensity=lifecycle::deprecated(),
                         smoothSigma=NA,
                         minSize=10,
                         kernel=lifecycle::deprecated(),
                         type=lifecycle::deprecated()) {

    if (lifecycle::is_present(minDensity)) {
        lifecycle::deprecate_warn(
            when = "0.2",
            what = "generateMask(minDensity)",
            details = paste("minDensity is not used anymore.",
                            "If you need to expand the borders, use `expand` argument instead.")
        )
    }

    if (lifecycle::is_present(kernel)) {
        lifecycle::deprecate_warn(
            when = "0.2",
            what = "generateMask(kernel)"
        )
    }

    if (lifecycle::is_present(type)) {
        lifecycle::deprecate_warn(
            when = "0.2",
            what = "generateMask(type)",
            details = paste("Independent mask generation is not supported anymore",
                            "Please contact the maintainer if you need this argument to be returned.")
        )
    }


    if (!is.na(smoothSigma)) {
        lifecycle::deprecate_soft(
            when = "0.2",
            what = "generateMask(smoothSigma)",
            details = paste("Automatic calculation of smoothSigma should work in most cases.",
                            "The argument will be fully deprecated, unless an example comes up where it's useful.",
                            "Please contact the maintainer if you need this argument to be kept.")
        )
    }

    clusterLevels <- unique(clusters)

    dims <- dims[, 1:2]
    if (is.null(colnames(dims))) {
        colnames(dims) <- c("x", "y")
    }

    window <- makeGridWindow(dims, gridSize=gridSize)

    pixelSize <- window$xstep
    smoothSigma <- smoothSigma * sqrt(area(window))
    expand <- expand * sqrt(area(window))
    windowHD <- makeGridWindow(dims, gridSize=max(gridSize, 1000))


    points <- spatstat.geom::ppp(dims[, 1], dims[, 2], window=window)


    allDensities <- lapply(clusterLevels, function(cluster) {
        res <- spatstat.geom::pixellate(points[clusters == cluster], xy=window)
        res
    })

    # getting initial masks
    curMasks <- lapply(seq_along(clusterLevels), function(i) {
        partPoints <- points[clusters == clusterLevels[i]]

        partSigma <- sqrt(bw.nrd(partPoints$x) * bw.nrd(partPoints$y)) * 1.5
        if (!is.na(smoothSigma)) {
            partSigma <- sqrt(partSigma * smoothSigma)
        }

        partMask <- pixellate(partPoints, xy=window)
        partMask[partMask == 0] <- NA
        partMask <- as.owin(partMask)
        partMaskV <- dilation(partMask, r = 2*partSigma + 1.5*pixelSize, polygonal=T)
        partMaskV <- erosion(partMaskV, r = 2*partSigma, polygonal=T)
        partMask <- as.mask(partMaskV, xy=window)
        partMask
    })

    nIter <- 3

    for (iter in seq_len(nIter)) {
        allDensitiesSmoothed <- lapply(seq_along(clusterLevels), function(i) {
            # message(i)
            curMask <- curMasks[[i]]
            curDensity <- allDensities[[i]]

            smoothed <- spatstat.geom::as.im(window) * 0

            if (area(curMask) == 0) {
                # lost the cluster, don't do anything
                return(smoothed)
            }

            parts <- getConnectedParts(curMask, curDensity, minSize = minSize)

            curPoints <- points[clusters == clusterLevels[i]]


            if (iter == nIter) {
                # smoothed <- spatstat.geom::as.im(windowHD) * 0
            }

            for (part in parts) {
                partPoints <- curPoints[part][window]

                partSigma <- sqrt(bw.nrd(partPoints$x) * bw.nrd(partPoints$y)) * 1.5
                if (!is.na(smoothSigma)) {
                    partSigma <- sqrt(partSigma * smoothSigma)
                }

                partPoints <- curPoints[dilation(part, r=2*partSigma)][window]

                partMask <- pixellate(partPoints, xy=window)
                partMask[partMask == 0] <- NA
                partMask <- as.owin(partMask)
                partMaskV <- dilation(partMask, r = 2*partSigma + 1.5*pixelSize, polygonal=T)
                partMaskV <- erosion(partMaskV, r = 2*partSigma, polygonal=T)
                partMask <- as.mask(partMaskV, xy=window)

                partBorder <- setminus.owin(
                    dilation(partMask, r=pixelSize*0, tight=FALSE),
                    erosion(partMask, r=pixelSize*1.5, tight=FALSE))
                partBorder <- intersect.owin(partBorder, window)

                partDensity <- density.ppp(partPoints, sigma=partSigma, xy=window)
                t <- median(partDensity[partBorder])

                if (iter == nIter) {
                    # better but slower way of smoothing borders
                    # partDensity <- density.ppp(partPoints[windowHD], sigma=partSigma, xy=windowHD)
                }

                smoothed <- smoothed + partDensity*(partDensity > t)
            }
            smoothed
        })

        backgroundDensityOne <- spatstat.geom::as.im(as.owin(allDensitiesSmoothed[[1]]))

        whichMaxDensity <- spatstat.geom::im.apply(
            c(list(backgroundDensityOne*0.01), allDensitiesSmoothed), which.max) - 1

        curMasks <- splitWhichMaxLevels(whichMaxDensity, nLevels=length(clusterLevels))
    }

    # smooth borders and expand a little (in vector)
    # TODO: important details can be removed here
    curMasks <- lapply(curMasks, closing, r=10*pixelSize, polygonal=TRUE)

    curMasks <- lapply(curMasks, dilation, r=expand, polygonal=TRUE)

    # switch to high-res
    curMasks <- lapply(curMasks, as.mask, xy = windowHD)

    curMasks <- removeMaskIntersections(curMasks, windowHD)

    borderTable <- rbindlist(lapply(seq_along(clusterLevels), function(i) {
        curMask <- curMasks[[i]]
        if (area(curMask) == 0) {
            warning(sprintf("Mask is empty for cluster %s", clusterLevels[i]))
            return(NULL)
        }
        curTable <- borderTableFromMask(curMask, allDensities[[i]])
        curTable[, cluster := clusterLevels[i]]
        curTable[, part := paste0(cluster, "#", part)]
        curTable[, group := paste0(part, "#", group)]
        curTable[]
    }))

    setnames(borderTable, c("x", "y"), colnames(dims))

    return(borderTable)
}

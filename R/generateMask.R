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
borderTableFromMask <- function(curMask, crop=TRUE) {
    if (crop) {
        toCrop <- boundingbox(curMask)
        # to avoid one-pixel crops https://github.com/spatstat/spatstat.geom/issues/24
        toCrop <- expandRect(toCrop, curMask$xstep, curMask$ystep)
        curMask <- curMask[toCrop]
    }

    polyMask <- as.polygonal(curMask)
    parts <- xypolycomponents(polyMask)

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

# weights should share coordinates with W
splitByMaxWeight <- function(weights, W) {
    maxWeight <- as.im(0, W=W)

    for (i in seq_along(weights)) {
        ix <- match(weights[[i]]$xcol, W$xcol)
        iy <- match(weights[[i]]$yrow, W$yrow)
        stopifnot(!anyNA(ix))
        stopifnot(!anyNA(iy))

        maxWeight$v[iy, ix] <- pmax(maxWeight$v[iy, ix], weights[[i]]$v)
    }

    res <- lapply(seq_along(weights), function(i) {
        ix <- match(weights[[i]]$xcol, W$xcol)
        iy <- match(weights[[i]]$yrow, W$yrow)
        solutionset(weights[[i]] > 0 & (weights[[i]] == maxWeight$v[iy, ix]))
    })
}

# curMasks should share coordinates with W
removeMaskIntersections <- function(curMasks, W) {
    maskWeights <- lapply(seq_along(curMasks), function(i) {
        D <- distmap(complement.owin(curMasks[[i]]))
        # a bit of jitter to make numbers unique
        D$v <- D$v * (1 + runif(length(D$v))*(2**-20))
        updateGridFrom(D, curMasks[[i]])
    })

    res <- splitByMaxWeight(maskWeights, W=W)
    res
}

expandRect <- function(rw, dx = 0, dy = dx) {
    stopifnot(rw$type == "rectangle")

    rw$xrange <- rw$xrange + c(-dx, dx)
    rw$yrange <- rw$yrange + c(-dy, dy)

    rw
}

updateGridFrom <- function(to, from, eps=1e-10) {
    stopifnot(identical(dim(to), dim(from)))
    stopifnot(max(abs(to$xrange-from$xrange)) < eps)
    stopifnot(max(abs(to$yrange-from$yrange)) < eps)

    to$xcol <- from$xcol
    to$yrow <- from$yrow
    to$xstep <- from$xstep
    to$ystep <- from$ystep
    to
}


getConnectedParts <- function(curMask, curDensity, minSize, absolutelyMinSize=5) {
    toCrop <- boundingbox(curMask)
    # to avoid one-pixel crops https://github.com/spatstat/spatstat.geom/issues/24
    toCrop <- expandRect(toCrop, curMask$xstep, curMask$ystep)
    curMask <- curMask[toCrop]

    curMaskSplit <- connected(curMask)
    parts <- tiles(tess(image=curMaskSplit))

    partSizes <- vapply(parts, function(part) {
        sum(curDensity[part])
    }, FUN.VALUE = numeric(1))

    parts <- parts[partSizes >= min(minSize, max(c(partSizes, absolutelyMinSize)))]
    unname(parts)
}


getRoughMask <- function(partPoints, window, partSigma, pixelSize, crop=TRUE) {
    extD <- 2*partSigma + 1.5*pixelSize

    if (crop) {
        toCrop <- boundingbox(partPoints)
        toCrop <- expandRect(toCrop, extD, extD)
        partPoints <- partPoints[toCrop]
        window <- window[toCrop]
    }

    D <- distmap(partPoints)
    partMaskV <- solutionset(D <= 2*partSigma + 1.5*pixelSize)
    partMaskV <- erosion(partMaskV, r = 2*partSigma, polygonal=F)
    partMask <- as.mask(partMaskV, xy=window)
    partMask <- updateGridFrom(partMask, window)
    partMask
}

# values below automatic threshold are clipped out
getPartDensityClipped <- function(curPoints, part, window, smoothSigma, pixelSize) {
    partPoints <- curPoints[part]

    partSigma <- sqrt(bw.nrd(partPoints$x) * bw.nrd(partPoints$y)) * 1.5
    if (!is.na(smoothSigma)) {
        partSigma <- sqrt(partSigma * smoothSigma)
    }

    extPart <- dilation(part, r=2*partSigma)
    partPoints <- curPoints[extPart]

    toCrop <- boundingbox(partPoints)
    extD <- 2*partSigma + 1.5*pixelSize
    toCrop <- expandRect(toCrop, extD, extD)
    partPoints <- partPoints[toCrop]

    partMask <- getRoughMask(partPoints, window[toCrop], partSigma, pixelSize, crop=FALSE)


    if (min(partMask$dim) <= 3) {
        partBorder <- partMask
    } else {
        #empty mask warning
        suppressWarnings(partMaskShrinked <- erosion(partMask, r=pixelSize*2, shrink.frame = FALSE))
        partBorder <- setminus.owin(
            partMask,
            partMaskShrinked
        )
    }

    # TODO: use blur on precomputed raw density instead
    partDensity <- density.ppp(partPoints, sigma=partSigma, xy=window[toCrop])
    t <- median(partDensity[partBorder])


    partDensity <- partDensity*(partDensity > t)
    partDensity <- as.im(partDensity, W = window, na.replace = 0)
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
#' @importFrom data.table rbindlist data.table setnames :=
#' @importFrom utils head tail
#' @importFrom stats median bw.nrd
#' @import spatstat.geom spatstat.explore
#' @export
#' @examples
#' data("exampleMascarade")
#' maskTable <- generateMask(dims=exampleMascarade$dims,
#'                           clusters=exampleMascarade$clusters)
#' data <- data.frame(exampleMascarade$dims,
#'                    cluster=exampleMascarade$clusters,
#'                    exampleMascarade$features)
#' library(ggplot2)
#' ggplot(data, aes(x=UMAP_1, y=UMAP_2)) +
#'     geom_point(aes(color=cluster)) +
#'     geom_path(data=maskTable, aes(group=group)) +
#'     coord_fixed() +
#'     theme_classic()
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
    windowHD <- makeGridWindow(dims, gridSize=gridSize*5)


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

        partMask <- getRoughMask(partPoints, window, partSigma, pixelSize)
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


            for (part in parts) {
                partDensity <- getPartDensityClipped(
                    curPoints, part, window, smoothSigma, pixelSize)

                smoothed <- smoothed + partDensity
            }
            smoothed
        })

        curMasks <- splitByMaxWeight(allDensitiesSmoothed, W=window)
    }

    # smooth borders and expand a little (in vector)
    # TODO: important details can be removed here
    curMasks <- lapply(curMasks, closing, r=10*pixelSize, polygonal=TRUE)

    curMasks <- lapply(curMasks, dilation, r=expand, polygonal=TRUE)

    # switch to high-res
    curMasks <- lapply(curMasks, function(m) {
        bbox <- boundingbox(m)
        subW <- windowHD[bbox]
        res <- as.mask(m, xy = subW)

        # force the coordinate grid to be the same (for some reason they're a bit different)
        res <- updateGridFrom(res, from=subW)
    })


    curMasks <- removeMaskIntersections(curMasks, W=windowHD)

    borderTable <- rbindlist(lapply(seq_along(clusterLevels), function(i) {
        curMask <- curMasks[[i]]
        if (area(curMask) == 0) {
            warning(sprintf("Mask is empty for cluster %s", clusterLevels[i]))
            return(NULL)
        }
        curTable <- borderTableFromMask(curMask, crop=FALSE)
        curTable[, cluster := clusterLevels[i]]
        curTable[, part := paste0(cluster, "#", part)]
        curTable[, group := paste0(part, "#", group)]
        curTable[]
    }))

    setnames(borderTable, c("x", "y"), colnames(dims))

    return(borderTable)
}

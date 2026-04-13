library(ggplot2)

make_mask_table <- function() {
    data("exampleMascarade")
    generateMask(dims = exampleMascarade$dims,
                 clusters = exampleMascarade$clusters,
                 gridSize = 50)
}

get_layer_colors <- function(res) {
    # The geom layer is the second element (first is coord_cartesian).
    res[[2]]$aes_params$colour
}


test_that("fancyMask applies single color to all clusters", {
    mt <- make_mask_table()
    res <- fancyMask(mt, cols = "red")
    colors <- get_layer_colors(res)
    expect_true(all(colors == "red"))
})

test_that("fancyMask maps unnamed color vector by factor level order", {
    mt <- make_mask_table()
    clusterLevels <- levels(mt$cluster)
    col_vec <- rep("grey50", length(clusterLevels))
    col_vec[1] <- "blue"
    res <- fancyMask(mt, cols = col_vec)
    colors <- get_layer_colors(res)
    expected_pal <- setNames(col_vec, clusterLevels)
    expect_equal(unname(colors), unname(expected_pal[mt$cluster]))
})

test_that("fancyMask maps named color vector by cluster name", {
    mt <- make_mask_table()
    clusterLevels <- levels(mt$cluster)
    # Provide in reverse order to verify matching is by name, not position
    col_vec <- setNames(scales::hue_pal()(length(clusterLevels)),
                        rev(clusterLevels))
    res <- fancyMask(mt, cols = col_vec)
    colors <- get_layer_colors(res)
    expect_equal(unname(colors), unname(col_vec[as.character(mt$cluster)]))
})

test_that("fancyMask errors on wrong-length unnamed cols", {
    mt <- make_mask_table()
    expect_error(fancyMask(mt, cols = c("red", "blue")),
                 "must equal the number of clusters")
})

test_that("fancyMask errors on named cols missing a cluster", {
    mt <- make_mask_table()
    clusterLevels <- levels(mt$cluster)
    # Provide only the first two clusters, omitting the rest
    col_vec <- setNames(c("red", "blue"), clusterLevels[1:2])
    expect_error(fancyMask(mt, cols = col_vec),
                 "missing entries for cluster")
})

test_that("fancyMask works when cluster column is character, not factor", {
    data("exampleMascarade")
    # Use character clusters (not factor) to mimic Seurat metadata
    charClusters <- paste0("_", exampleMascarade$clusters)
    mt <- generateMask(dims = exampleMascarade$dims,
                       clusters = charClusters,
                       gridSize = 50)
    expect_false(is.factor(mt$cluster))

    # Direct cols should work
    res <- fancyMask(mt, cols = "red")
    colors <- get_layer_colors(res)
    expect_true(length(colors) > 0)
    expect_true(all(colors == "red"))

    # Inherit with a discrete scale should work
    plotData <- data.frame(exampleMascarade$dims,
                           cluster = factor(charClusters))
    p <- ggplot(plotData) +
        geom_point(aes(x = UMAP_1, y = UMAP_2, color = cluster), size = 0.5) +
        fancyMask(mt, cols = "inherit")
    b <- ggplot_build(p)
    mask_colours <- b$data[[length(p$layers)]]$colour
    expect_true(length(mask_colours) > 0)
    expect_true(!any(is.na(mask_colours)))
})

test_that("fancyMask renders without error with custom cols", {
    mt <- make_mask_table()
    data("exampleMascarade")
    p <- ggplot(do.call(cbind, exampleMascarade)) +
        geom_point(aes(x = UMAP_1, y = UMAP_2)) +
        fancyMask(mt, cols = "red") +
        theme_classic()
    pf <- tempfile(fileext = ".pdf")
    expect_no_error(ggsave(p, file = pf, width = 5, height = 4))
    expect_true(file.exists(pf))
})

test_that("fancyMask with cols='inherit' returns a plain list", {
    mt <- make_mask_table()
    res <- fancyMask(mt, cols = "inherit")
    expect_true(is.list(res))
    expect_false(inherits(res, "fancyMask"))
})

test_that("cols='inherit' picks up scale_color_manual values", {
    mt <- make_mask_table()
    data("exampleMascarade")
    clusterLevels <- levels(mt$cluster)
    myPal <- setNames(rainbow(length(clusterLevels)), clusterLevels)

    plotData <- as.data.frame(do.call(cbind, exampleMascarade))
    plotData$clusters <- factor(plotData$clusters, levels = clusterLevels)
    p <- ggplot(plotData) +
        geom_point(aes(x = UMAP_1, y = UMAP_2, color = clusters)) +
        scale_color_manual(values = myPal) +
        fancyMask(mt, cols = "inherit")

    # Colours resolved at build time via ggplot2's scale machinery
    b <- ggplot_build(p)
    mask_colours <- b$data[[length(p$layers)]]$colour
    expected <- unname(myPal[as.character(mt$cluster)])
    expect_equal(mask_colours, expected)
})

test_that("cols='inherit' uses ggplot2 default scale when no colour scale is set", {
    mt <- make_mask_table()
    data("exampleMascarade")

    # No colour aesthetic on the scatter — mask uses its own aes(colour=cluster)
    # and gets ggplot2's default discrete scale, not black.
    p <- ggplot(do.call(cbind, exampleMascarade)) +
        geom_point(aes(x = UMAP_1, y = UMAP_2)) +
        fancyMask(mt, cols = "inherit")

    b <- ggplot_build(p)
    mask_colours <- b$data[[length(p$layers)]]$colour
    expect_true(length(mask_colours) > 0)
    expect_false(all(mask_colours == "black"))
})

test_that("cols accepts a palette function", {
    mt <- make_mask_table()
    res <- fancyMask(mt, cols = rainbow)
    colors <- get_layer_colors(res)
    clusterLevels <- levels(mt$cluster)
    expected <- setNames(rainbow(length(clusterLevels)), clusterLevels)
    expect_equal(unname(colors), unname(expected[mt$cluster]))
})

test_that("cols='inherit' picks up colors from layer-level mapping", {
    mt <- make_mask_table()
    clusterLevels <- levels(mt$cluster)
    myPal <- setNames(rainbow(length(clusterLevels)), clusterLevels)

    # Color mapping is on the layer, not on the plot
    p <- ggplot(data.frame(x = 1, y = 1,
                           cl = factor(clusterLevels, levels = clusterLevels))) +
        geom_point(aes(x = x, y = y, color = cl)) +
        scale_color_manual(values = myPal) +
        fancyMask(mt, cols = "inherit")

    b <- ggplot_build(p)
    mask_colours <- b$data[[length(p$layers)]]$colour
    expected <- unname(myPal[as.character(mt$cluster)])
    expect_equal(mask_colours, expected)
})

make_two_part_mask_table <- function() {
    # Cluster "A" with two parts: part #1 (large square) and part #2 (small square)
    # Cluster "B" with one part
    data.frame(
        x       = c(0, 2, 2, 0, 0,   5, 5.5, 5.5, 5, 5,   10, 11, 11, 10, 10),
        y       = c(0, 0, 2, 2, 0,   0, 0,   0.5, 0.5, 0,  0,  0,  1,  1,  0),
        cluster = factor(c(rep("A", 10), rep("B", 5))),
        part    = c(rep("A#1", 5), rep("A#2", 5), rep("B#1", 5)),
        group   = c(rep("A#1#1", 5), rep("A#2#1", 5), rep("B#1#1", 5)),
        stringsAsFactors = FALSE
    )
}

test_that("label.largest=TRUE sets non-first-part labels to NA", {
    mt  <- make_two_part_mask_table()
    res <- fancyMask(mt, cols = "black", label = TRUE, label.largest = TRUE)
    layer_data <- res[[2]]$data
    expect_true(all(is.na(layer_data$.label_display[layer_data$part == "A#2"])))
    expect_false(any(is.na(layer_data$.label_display[layer_data$part == "A#1"])))
    expect_false(any(is.na(layer_data$.label_display[layer_data$part == "B#1"])))
})

test_that("label.largest=FALSE labels all parts", {
    mt  <- make_two_part_mask_table()
    res <- fancyMask(mt, cols = "black", label = TRUE, label.largest = FALSE)
    layer_data <- res[[2]]$data
    expect_false(any(is.na(layer_data$.label_display)))
})

test_that("label.largest=TRUE works with cols='inherit'", {
    mt  <- make_two_part_mask_table()
    res <- fancyMask(mt, cols = "inherit", label = TRUE, label.largest = TRUE)
    layer_data <- res[[2]]$data
    expect_true(all(is.na(layer_data$.label_display[layer_data$part == "A#2"])))
    expect_false(any(is.na(layer_data$.label_display[layer_data$part == "A#1"])))
})

test_that("cols='inherit' renders without error", {
    mt <- make_mask_table()
    data("exampleMascarade")
    clusterLevels <- levels(mt$cluster)
    myPal <- setNames(rainbow(length(clusterLevels)), clusterLevels)

    plotData <- as.data.frame(do.call(cbind, exampleMascarade))
    # cbind coerces the factor to integer codes before factor() runs, so levels
    # would be "1","2",... instead of cluster names — pass levels explicitly
    plotData$clusters <- factor(plotData$clusters, levels = clusterLevels)
    p <- ggplot(plotData) +
        geom_point(aes(x = UMAP_1, y = UMAP_2, color = clusters)) +
        scale_color_manual(values = myPal) +
        fancyMask(mt, cols = "inherit") +
        theme_classic()
    pf <- tempfile(fileext = ".pdf")
    expect_no_error(ggsave(p, file = pf, width = 5, height = 4))
    expect_true(file.exists(pf))
})

test_that("cols='inherit' border colors match point colors when factor has unused levels", {
    # Ghost levels must not inflate the palette: ggplot2's drop=TRUE means only
    # real clusters are trained, so mask and scatter colours must agree.
    data("exampleMascarade")

    clusters_with_ghosts <- factor(
        exampleMascarade$clusters,
        levels = c(levels(exampleMascarade$clusters), "ghost1", "ghost2")
    )

    mt <- generateMask(dims = exampleMascarade$dims,
                       clusters = clusters_with_ghosts,
                       gridSize = 50)

    expect_true("ghost1" %in% levels(mt$cluster))
    expect_false("ghost1" %in% as.character(mt$cluster))

    plotData <- data.frame(exampleMascarade$dims, cluster = clusters_with_ghosts)
    p <- ggplot(plotData) +
        geom_point(aes(x = UMAP_1, y = UMAP_2, color = cluster), size = 0.5) +
        fancyMask(mt, cols = "inherit")

    pb           <- ggplot_build(p)
    pt_colours   <- pb$data[[1]]$colour
    mask_colours <- pb$data[[2]]$colour

    cluster_char <- as.character(clusters_with_ghosts)
    for (cl in unique(as.character(mt$cluster))) {
        cell_idx <- which(cluster_char == cl)[1]
        mask_idx <- which(as.character(mt$cluster) == cl)[1]
        expect_equal(mask_colours[mask_idx], pt_colours[cell_idx],
                     label = paste("cluster", cl))
    }
})

test_that("cols='inherit' mask colours are identical regardless of scale_color_* order", {
    mt <- make_mask_table()
    data("exampleMascarade")
    clusterLevels <- levels(mt$cluster)
    myPal <- setNames(rainbow(length(clusterLevels)), clusterLevels)

    plotData <- as.data.frame(do.call(cbind, exampleMascarade))
    plotData$clusters <- factor(plotData$clusters, levels = clusterLevels)

    # scale before fancyMask
    p1 <- ggplot(plotData) +
        geom_point(aes(x = UMAP_1, y = UMAP_2, color = clusters)) +
        scale_color_manual(values = myPal) +
        fancyMask(mt, cols = "inherit")

    # scale after fancyMask — this was the broken case
    p2 <- ggplot(plotData) +
        geom_point(aes(x = UMAP_1, y = UMAP_2, color = clusters)) +
        fancyMask(mt, cols = "inherit") +
        scale_color_manual(values = myPal)

    mask_idx <- length(p1$layers)
    b1_colours <- ggplot_build(p1)$data[[mask_idx]]$colour
    b2_colours <- ggplot_build(p2)$data[[mask_idx]]$colour

    expect_equal(b1_colours, b2_colours)
})

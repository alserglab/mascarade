library(ggplot2)

make_mask_table <- function() {
    data("exampleMascarade")
    generateMask(dims = exampleMascarade$dims,
                 clusters = exampleMascarade$clusters,
                 gridSize = 50)
}

get_layer_colors <- function(res) {
    # The geom layer is the second element (first is coord_cartesian)
    # Note: "colour" is the ggplot2 internal field name
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
    mask_layer <- p$layers[[length(p$layers)]]
    colors <- mask_layer$aes_params$colour
    expect_true(length(colors) > 0)
    expect_true(!any(is.na(colors)))
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

test_that("fancyMask with cols='inherit' returns a fancyMask S3 object", {
    mt <- make_mask_table()
    res <- fancyMask(mt, cols = "inherit")
    expect_s3_class(res, "fancyMask")
})

test_that("cols='inherit' picks up scale_color_manual values", {
    mt <- make_mask_table()
    data("exampleMascarade")
    clusterLevels <- levels(mt$cluster)
    myPal <- setNames(rainbow(length(clusterLevels)), clusterLevels)

    plotData <- as.data.frame(do.call(cbind, exampleMascarade))
    plotData$clusters <- factor(plotData$clusters)
    p <- ggplot(plotData) +
        geom_point(aes(x = UMAP_1, y = UMAP_2, color = clusters)) +
        scale_color_manual(values = myPal) +
        fancyMask(mt, cols = "inherit")

    # The mask layer is the last layer added by ggplot_add
    mask_layer <- p$layers[[length(p$layers)]]
    colors <- mask_layer$aes_params$colour
    expected <- myPal[as.character(mt$cluster)]
    expect_equal(unname(colors), unname(expected))
})

test_that("cols='inherit' falls back to black without a color scale", {
    mt <- make_mask_table()
    data("exampleMascarade")

    p <- ggplot(do.call(cbind, exampleMascarade)) +
        geom_point(aes(x = UMAP_1, y = UMAP_2)) +
        fancyMask(mt, cols = "inherit")

    mask_layer <- p$layers[[length(p$layers)]]
    colors <- mask_layer$aes_params$colour
    expect_true(all(colors == "black"))
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

    mask_layer <- p$layers[[length(p$layers)]]
    colors <- mask_layer$aes_params$colour
    expected <- myPal[as.character(mt$cluster)]
    expect_equal(unname(colors), unname(expected))
})

test_that("cols='inherit' renders without error", {
    mt <- make_mask_table()
    data("exampleMascarade")
    clusterLevels <- levels(mt$cluster)
    myPal <- setNames(rainbow(length(clusterLevels)), clusterLevels)

    plotData <- as.data.frame(do.call(cbind, exampleMascarade))
    plotData$clusters <- factor(plotData$clusters)
    p <- ggplot(plotData) +
        geom_point(aes(x = UMAP_1, y = UMAP_2, color = clusters)) +
        scale_color_manual(values = myPal) +
        fancyMask(mt, cols = "inherit") +
        theme_classic()
    pf <- tempfile(fileext = ".pdf")
    expect_no_error(ggsave(p, file = pf, width = 5, height = 4))
    expect_true(file.exists(pf))
})

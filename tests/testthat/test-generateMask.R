test_that("generateMask works on example data", {
    data("exampleMascarade")
    res <- generateMask(dims=exampleMascarade$dims,
                        clusters=exampleMascarade$clusters,
                        gridSize=50)
    expect_true(!is.null(res))
})

test_that("generateMask errors when clusters length does not match nrow(dims)", {
    set.seed(42)
    dims <- matrix(rnorm(90 * 2), ncol = 2)
    clusters <- rep(c("A", "B"), length.out = 100)

    expect_error(
        generateMask(dims, clusters),
        "length(clusters) must equal nrow(dims)",
        fixed = TRUE
    )
})

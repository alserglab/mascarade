test_that("generateMask works on example data", {
    data("exampleMascarade")
    res <- generateMask(dims=exampleMascarade$dims,
                        clusters=exampleMascarade$clusters,
                        gridSize=50)
    expect_true(!is.null(res))
})

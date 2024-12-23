test_that("generateMask works on example data", {
    data("exampleMascarade")
    res <- generateMask(dims=exampleMascarade$dims,
                        clusters=exampleMascarade$clusters)
    expect_true(!is.null(res))
})

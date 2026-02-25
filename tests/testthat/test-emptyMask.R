test_that("generateMask handles cluster absorbed by neighbor (reproduces boundingbox error)", {
    # Cluster B has 2 identical points at the center of cluster A.
    # bw.nrd returns 0 for identical values, so the initial mask is a tiny disc.
    # During iterations, A's density dominates and B loses all pixels.
    # The resulting empty polygon causes boundingbox() to fail with:
    #   "xrange should be a vector of length 2 giving (xmin, xmax)"
    set.seed(42)
    dims <- rbind(
        matrix(rnorm(200, mean = 0, sd = 3), ncol = 2),      # cluster A: 100 points
        matrix(c(0, 0, 0, 0), ncol = 2, byrow = TRUE)         # cluster B: 2 identical points inside A
    )
    clusters <- c(rep("A", 100), "B", "B")

    # Should not throw: Error in owinInternalRect(...) : xrange should be a vector
    # of length 2 giving (xmin, xmax)
    expect_no_error(generateMask(dims, clusters))
})

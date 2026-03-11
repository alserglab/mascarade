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
    #
    # Two warnings are expected and suppressed:
    # - "data contain duplicated points": spatstat warns about the two identical
    #   (0, 0) points — intentional in this test setup
    # - "Mask is empty for cluster B": cluster B gets absorbed by A during
    #   iteration, which is exactly the scenario this test documents
    expect_no_error(suppressWarnings(generateMask(dims, clusters)))
})

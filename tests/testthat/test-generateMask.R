test_that("generateMask works on example data", {
    data("exampleMascarade")
    res <- generateMask(dims=exampleMascarade$dims,
                        clusters=exampleMascarade$clusters,
                        gridSize=50)
    expect_true(!is.null(res))
})

test_that("generateMask part #1 is the largest part by row count", {
    set.seed(42)
    # Cluster "A": large group (100 pts) far from small group (20 pts)
    large <- matrix(c(rnorm(100, mean = 0), rnorm(100, mean = 0)), ncol = 2)
    small <- matrix(c(rnorm(20,  mean = 10), rnorm(20,  mean = 10)), ncol = 2)
    # Cluster "B": single compact group so two clusters exist
    other <- matrix(c(rnorm(60,  mean = 5), rnorm(60,  mean = -5)), ncol = 2)

    dims     <- rbind(large, small, other)
    clusters <- c(rep("A", 120), rep("B", 60))

    res <- generateMask(dims, clusters, gridSize = 50, minSize = 5)

    clusterA <- res[res$cluster == "A", ]
    parts <- unique(clusterA$part)

    if (length(parts) >= 2) {
        rowsPerPart <- tapply(seq_len(nrow(clusterA)), clusterA$part, length)
        # Part ending in #1 must have at least as many rows as any other part
        part1rows  <- rowsPerPart[grep("#1$", names(rowsPerPart))]
        otherRows  <- rowsPerPart[!grepl("#1$", names(rowsPerPart))]
        expect_true(all(part1rows >= otherRows))
    }
})

test_that("generateMask warns and skips mask for single-point cluster", {
    data("exampleMascarade")
    dims2 <- rbind(exampleMascarade$dims, c(0, 0))
    clusters2 <- c(as.character(exampleMascarade$clusters), "singleton")

    expect_warning(
        res <- generateMask(dims=dims2, clusters=clusters2, gridSize=50),
        "singleton.*fewer than two points"
    )
    expect_false("singleton" %in% res$cluster)
})

test_that("generateMask handles cluster where all points share a coordinate (partSigma=0)", {
    set.seed(1)
    big <- matrix(c(rnorm(200), rnorm(200)), ncol = 2)
    # 3 close points all at same y → bw.nrd(y)=0 → partSigma=0 without fix
    small <- matrix(c(c(-0.1, 0, 0.1), c(8, 8, 8)), ncol = 2)
    dims <- rbind(big, small)
    clusters <- c(rep("A", 200), rep("B", 3))

    expect_no_error(generateMask(dims, clusters, gridSize = 50, minSize = 1))
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

library(ggplot2)

# Regression test: a polygon thin enough to be completely contracted by polyoffset(-1pt)
# is silently dropped by ggforce's shapeGrob makeContent, which leaves a gap between
# x$labeldim (n entries) and the rebuilt polygons (n-1 entries), causing
# anchors[[i]] to be out of bounds in my_place_labels.
test_that("geom_mark_shape does not crash when a thin polygon is eliminated by expansion", {
    # At ggsave(width=5, height=4), plot panel is ~112mm wide.
    # xlim range = 14 → scale ≈ 8 mm/unit.
    # thin polygon width = 0.05 units = ~0.4mm < 0.706mm (= 2 × 1pt) → eliminated.
    thin <- data.frame(x = c(0, 0.05, 0.05, 0), y = c(0, 0, 5, 5),
                       group = "thin", label = "thin")
    A    <- data.frame(x = c(3, 6, 6, 3),        y = c(0, 0, 5, 5),
                       group = "A",    label = "A")
    B    <- data.frame(x = c(8, 11, 11, 8),       y = c(0, 0, 5, 5),
                       group = "B",    label = "B")
    df <- rbind(thin, A, B)

    p <- ggplot(df, aes(x = x, y = y, group = group, label = label, color = group)) +
        geom_mark_shape(expand = unit(-1, "pt")) +
        xlim(-1, 13) + ylim(-1, 7)

    pf <- tempfile(fileext = ".png")
    expect_no_error(ggsave(p, file = pf, width = 5, height = 4))
    expect_true(file.exists(pf))
})

test_that("geom_mark_shape works", {
    shape1 <- data.frame(
        x = c(0, 3, 3, 2, 2, 1, 1, 0),
        y = c(0, 0, 3, 3, 1, 1, 3, 3),
        label="bracket"
    )
    shape2 <- data.frame(
        x = c(0, 3, 3, 0)+4,
        y = c(0, 0, 3, 3),
        label="square"
    )
    shape3 <- data.frame(
        x = c(0, 1.5, 3, 1.5)+8,
        y = c(1.5, 0, 1.5, 3),
        label="diamond"
    )

    p <- ggplot(rbind(shape1, shape2, shape3), aes(x=x, y=y, label=label, color=label, fill=label)) +
        geom_mark_shape() +
        ylim(0, 5)


    pf <- tempfile(fileext = ".pdf")
    expect_no_error(ggsave(p, file=pf, width=5, height=4))
    expect_true(file.exists(pf))
})

test_that("geom_mark_shape works with various simp_ratio values", {
    shape1 <- data.frame(
        x = c(0, 3, 3, 2, 2, 1, 1, 0),
        y = c(0, 0, 3, 3, 1, 1, 3, 3),
        label="bracket"
    )
    shape2 <- data.frame(
        x = c(0, 3, 3, 0)+4,
        y = c(0, 0, 3, 3),
        label="square"
    )
    shape3 <- data.frame(
        x = c(0, 1.5, 3, 1.5)+8,
        y = c(1.5, 0, 1.5, 3),
        label="diamond"
    )
    df <- rbind(shape1, shape2, shape3)

    for (ratio in c(0, 0.001, 0.01)) {
        p <- ggplot(df, aes(x=x, y=y, label=label, color=label, fill=label)) +
            geom_mark_shape(simp_ratio = ratio) +
            ylim(0, 5)
        pf <- tempfile(fileext = ".pdf")
        expect_no_error(ggsave(p, file=pf, width=5, height=4))
    }
})

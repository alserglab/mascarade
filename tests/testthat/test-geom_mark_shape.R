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

test_that("border colours are not shifted after a thin polygon is dropped", {
    # Regression: when ggforce drops a polygon (contracts to nothing under negative
    # expand), mark$id contains non-consecutive ids. Grid recycles mark$gp$col
    # positionally, so without remapping, each surviving polygon gets the *previous*
    # polygon's colour. The fix: mark$gp <- mark$gp[surviving] inside makeContent.shape_enc.
    #
    # Setup: three groups — "thin" is eliminated, "A" and "B" survive.
    # Before the fix: A draws in thin's colour, B draws in A's colour.
    # After  the fix: A draws in A's colour,    B draws in B's colour.
    thin_col <- "#ff0000"
    A_col    <- "#00ff00"
    B_col    <- "#0000ff"

    thin <- data.frame(x = c(0, 0.05, 0.05, 0), y = c(0, 0, 5, 5),
                       group = "thin", label = "thin")
    A    <- data.frame(x = c(3, 6, 6, 3),        y = c(0, 0, 5, 5),
                       group = "A",    label = "A")
    B    <- data.frame(x = c(8, 11, 11, 8),       y = c(0, 0, 5, 5),
                       group = "B",    label = "B")
    df <- rbind(thin, A, B)
    # Force thin first so its drop shifts the gp colour vector — the bug only
    # manifests when the dropped polygon is not the last in the factor ordering.
    df$group <- factor(df$group, levels = c("thin", "A", "B"))

    p <- ggplot(df, aes(x = x, y = y, group = group, label = label, colour = group)) +
        geom_mark_shape(expand = unit(-1, "pt"), show.legend = FALSE) +
        scale_colour_manual(values = c(thin = thin_col, A = A_col, B = B_col)) +
        xlim(-1, 13) + ylim(-1, 7) +
        theme_void()

    svg_file <- tempfile(fileext = ".svg")
    svglite::svglite(svg_file, width = 5, height = 4)
    print(p)
    dev.off()

    svg_text <- paste(readLines(svg_file), collapse = "\n")

    expect_true(grepl(A_col, svg_text, ignore.case = TRUE),
                label = "A border colour present in SVG")
    expect_true(grepl(B_col, svg_text, ignore.case = TRUE),
                label = "B border colour present in SVG")
    # Without the fix, thin's colour leaks into A's border; verify it is absent.
    expect_false(grepl(thin_col, svg_text, ignore.case = TRUE),
                 label = "dropped polygon colour absent from SVG")
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

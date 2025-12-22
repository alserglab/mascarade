library(ggplot2)

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

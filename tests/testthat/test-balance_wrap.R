# Tests for balance_wrap(): the balanced label wrapper behind label.width.

# balance_wrap measures token widths, so it needs an open graphics device.
with_device <- function(expr) {
  pdf(NULL)
  on.exit(dev.off())
  force(expr)
}

gp <- grid::gpar(fontsize = 12, fontfamily = "", lineheight = 1)
lines_of <- function(x) strsplit(x, "\n", fixed = TRUE)[[1]]

test_that("a single token is returned unchanged", {
  with_device({
    expect_identical(balance_wrap("Aerocytes", gp, 50), "Aerocytes")
  })
})

test_that("an over-long single word is never broken", {
  with_device({
    out <- balance_wrap("Supercalifragilistic", gp, 5)
    expect_false(grepl("\n", out, fixed = TRUE))
  })
})

test_that("a hyphen breaks after (never before) and stays on the upper line", {
  with_device({
    # width = 1mm forces the only break point, which is the hyphen.
    lines <- lines_of(balance_wrap("Non-Classical", gp, 1))
    expect_length(lines, 2)
    expect_match(lines[1], "-$")            # hyphen ends the upper line
    expect_identical(lines[2], "Classical")
    expect_false(any(grepl("^-", lines)))   # no line starts with a hyphen
  })
})

test_that("a tiny leading word is not orphaned on its own line", {
  with_device({
    # Greedy wrapping would strand "T"; the balanced DP merges it upward.
    lines <- lines_of(balance_wrap("T senescence-associated", gp, 30))
    expect_gt(length(lines), 1)
    expect_false(identical(lines[1], "T"))
    expect_match(lines[1], "^T senescence")
  })
})

test_that("wrapping preserves the text content", {
  with_device({
    label <- "Macrophages smoking-induced"
    out <- balance_wrap(label, gp, 30)
    # Joining lines back (hyphen breaks add no space) reproduces the input.
    rejoined <- gsub("-\n", "-", out)       # undo hyphen breaks
    rejoined <- gsub("\n", " ", rejoined)   # undo space breaks
    expect_identical(rejoined, label)
  })
})

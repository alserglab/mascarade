# https://github.com/thomasp85/ggforce/blob/main/R/mark_label.R
# `place_labels` and `make_label` function were modified and moved to `mark_label.R`

# MIT License
#
# Copyright (c) 2019 Thomas Lin Pedersen
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
#     The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.



#' @importFrom grid gpar
#' @importFrom rlang list2 inject caller_env
inherit_gp <- function(..., gp, call = caller_env()) {
  new_gp <- list2(...)
  for (par in names(new_gp)) {
    old_par <- par
    inherited_par <- new_gp[[par]]
    if (isTRUE(new_gp[[par]] == 'inherit')) {
      inherited_par <- gp[[old_par]]
    } else if (isTRUE(new_gp[[par]] == 'inherit_fill')) {
      old_par <- 'fill'
      inherited_par <- gp[[old_par]]
    } else if (isTRUE(grepl('inherit_col', new_gp[[par]]))) {
      old_par <- 'col'
      inherited_par <- gp[[old_par]]
    }
    if (is.null(inherited_par)) {
      cli::cli_abort("Can't inherit {.field {old_par}} as it is not given in the root {.cls gpar}")
    }
    new_gp[[par]] <- inherited_par
  }
  inject(gpar(!!!new_gp))
}
subset_gp <- function(gp, index, ignore = c('font')) {
  gp_names <- names(gp)
  ignore_idx <- unique0(unlist(lapply(ignore, grep, gp_names)))
  if (length(ignore_idx) > 0) {
    gp_names <- gp_names[-ignore_idx]
  }
  for (par in gp_names) {
    gp[[par]] <- rep_len(gp[[par]], max(index))[index]
  }
  gp
}



#' @importFrom grid valid.just textGrob nullGrob viewport grobWidth grobHeight
#' rectGrob gpar grid.layout unit gTree gList grobDescent
labelboxGrob <- function(label, x = unit(0.5, 'npc'), y = unit(0.5, 'npc'),
                         description = NULL, width = NULL, maxwidth = NULL,
                         min.width = 50,
                         default.units = 'mm', hjust = 0,
                         pad = margin(2, 2, 2, 2, 'mm'), gp = gpar(), desc.gp = gpar(),
                         vp = NULL) {
  width <- as_mm(width, default.units)
  maxwidth <- as_mm(maxwidth, default.units)
  min.width <- as_mm(min.width, default.units)
  pad <- as_mm(pad, default.units)
  pad[c(1, 3)] <- as_mm(pad[c(1, 3)], default.units, FALSE)
  if (is.null(width) && !is.null(maxwidth)) {
    # A soft max width caps the box, so it also lowers the effective minimum
    # (otherwise a min.width above maxwidth would force the box past the cap).
    # A strict `width`, when given, takes precedence over `maxwidth`.
    min.width <- min(min.width, maxwidth)
  }
  if (!is.null(label) && !is.na(label)) {
    if (!is.null(width)) {
      label <- wrap_text(label, gp, width - pad[2] - pad[4])
    } else if (!is.null(maxwidth)) {
      label <- balance_wrap(label, gp, maxwidth - pad[2] - pad[4])
    }
    just <- c(hjust[1], 0.5)
    lab_grob <- textGrob(label, x = just[1], y = just[2], just = just,
                         gp = gp)
  } else {
    lab_grob <- nullGrob()
  }
  if (!is.null(width)) {
    final_width <- max(width, min.width) - pad[2] - pad[4]
  } else {
    if (as_mm(grobWidth(lab_grob)) > (min.width - pad[2] - pad[4])) {
      final_width <- as_mm(grobWidth(lab_grob)) + pad[2] + pad[4]
    } else {
      final_width <- max(as_mm(grobWidth(lab_grob)), min.width) - pad[2] - pad[4]
    }
  }
  if (!is.null(description) && !is.na(description)) {
    description <- if (is.null(width) && !is.null(maxwidth)) {
      balance_wrap(description, desc.gp, final_width)
    } else {
      wrap_text(description, desc.gp, final_width)
    }
    just <- c(rep_len(hjust, 2)[2], 0.5)
    desc_grob <- textGrob(description, x = just[1], y = just[2], just = just,
                          gp = desc.gp)
    if (is.null(width)) {
      final_width_desc <- min(final_width, as_mm(grobWidth(desc_grob)))
      final_width <- as_mm(grobWidth(lab_grob))
      if (final_width < final_width_desc) {
        final_width <- final_width_desc
      }
    }
  } else {
    desc_grob <- nullGrob()
    if (is.null(width)) final_width <- as_mm(grobWidth(lab_grob))
  }
  bg_grob <- rectGrob(gp = gpar(col = NA, fill = gp$fill))
  lab_height <- as_mm(grobHeight(lab_grob), width = FALSE)
  desc_height <- as_mm(grobHeight(desc_grob), width = FALSE)
  sep_height <- if (lab_height > 0 && desc_height > 0) {
    pad[1]
  } else if (lab_height > 0) {
    # Descenders extend below grobHeight(lab_grob) (which is the ascent box). Reserve only the
    # part the bottom margin does not already cover, so a label without descenders is not padded
    # with dead space underneath.
    max(0, font_descent(gp$fontfamily, gp$fontface, gp$fontsize, gp$cex) - pad[3])
  } else {
    0
  }
  vp <- viewport(
    x = x,
    y = y,
    width = unit(final_width + pad[2] + pad[4], 'mm'),
    height = unit(pad[1] + pad[3] + lab_height + desc_height + sep_height,
                  'mm'),
    layout = grid.layout(
      5, 3,
      widths = unit(c(pad[2], final_width, pad[4]), 'mm'),
      heights = unit(c(pad[1], lab_height, sep_height, desc_height, pad[3]),
                     'mm')
    )
  )
  lab_grob$vp <- viewport(layout.pos.col = 2, layout.pos.row = 2)
  desc_grob$vp <- viewport(layout.pos.col = 2, layout.pos.row = 4)
  gTree(children = gList(bg_grob, lab_grob, desc_grob), vp = vp,
        cl = 'mark_label')
}
#' @importFrom grid widthDetails
widthDetails.mark_label <- function(x) {
  x$vp$width
}
#' @importFrom grid heightDetails
heightDetails.mark_label <- function(x) {
  x$vp$height
}
#' @importFrom grid textGrob grobWidth
wrap_text <- function(text, gp, width) {
  text <- gsub('-', '- ', text)
  text <- strsplit(text, split = ' ', fixed = TRUE)[[1]]
  text <- paste0(text, ' ')
  text <- sub('- ', '-', text)
  txt <- ''
  for (i in text) {
    oldlab <- txt
    txt <- paste0(txt, i)
    tmpGrob <- textGrob(txt, gp = gp)
    if (as_mm(grobWidth(tmpGrob)) > width) {
      txt <- paste(trimws(oldlab), i, sep = '\n')
    }
  }
  trimws(txt)
}

#' Split text into wrap tokens
#'
#' Break points are spaces and positions just *after* a hyphen (never before),
#' matching `wrap_text()`. A hyphen stays attached to the token on its left, so
#' a break there leaves the hyphen at the end of the upper line.
#'
#' @param text A length-1 character string.
#' @return A character vector of tokens.
#' @keywords internal
#' @noRd
wrap_tokens <- function(text) {
  s <- gsub('-', '- ', text, fixed = TRUE)
  toks <- strsplit(s, split = ' ', fixed = TRUE)[[1]]
  toks[toks != '']
}

#' Balance wrapped text across lines to a soft target width
#'
#' Word-wraps `text` so the lines are as even as possible and close to `width`
#' (in mm), using a minimum-raggedness dynamic program. Unlike the greedy
#' `wrap_text()`, it avoids a short dangling final line and stray orphan words,
#' and it treats `width` as a *soft* target: a line may exceed it slightly when
#' that removes an orphan, and an over-long single word is never broken (the
#' effective target grows to the widest word so it fits on its own line). Break
#' points are the same as `wrap_text()` (spaces and after-hyphens).
#'
#' The line cost is `slack^2` for a line narrower than the target and
#' `over * excess^2` for one wider than it. Squared slack penalises a nearly
#' empty line about as much as the width itself, so a tiny leading word (e.g.
#' `"T"`) is merged with its neighbour rather than stranded. Larger `over`
#' tolerates less overflow (more, shorter lines).
#'
#' @param text A length-1 character string.
#' @param gp A `grid::gpar()` describing the font; used to measure token widths.
#' @param width Soft target line width, in mm.
#' @param over Overflow penalty weight. Default `8`.
#' @return A single string with lines separated by `\n`.
#' @keywords internal
#' @noRd
#' @importFrom grid textGrob grobWidth
balance_wrap <- function(text, gp, width, over = 8) {
  toks <- wrap_tokens(text)
  n <- length(toks)
  if (n <= 1) {
    return(trimws(text))
  }
  tokenWidth <- vapply(
    toks, function(t) as_mm(grobWidth(textGrob(t, gp = gp))), numeric(1))
  spaceWidth <- as_mm(grobWidth(textGrob('x x', gp = gp))) -
    2 * as_mm(grobWidth(textGrob('x', gp = gp)))
  # A space precedes each token that shares a line with its predecessor, unless
  # that predecessor ends in '-' (a mid-word break point takes no space).
  endsHyphen <- grepl('-$', toks)
  gapBefore <- c(0, ifelse(endsHyphen[-n], 0, spaceWidth))
  cumWidth <- cumsum(tokenWidth)
  cumGap <- cumsum(gapBefore)
  lineWidth <- function(a, b) {
    before <- 0
    if (a > 1) {
      before <- cumWidth[a - 1]
    }
    tokens <- cumWidth[b] - before
    gaps <- cumGap[b] - cumGap[a]              # internal gaps only
    tokens + gaps
  }
  target <- max(width, max(tokenWidth))        # unbreakable-word rule
  lineCost <- function(a, b) {
    slack <- target - lineWidth(a, b)
    if (slack >= 0) {
      slack^2
    } else {
      over * slack^2
    }
  }

  # DP over all line counts: cost[i + 1] = min raggedness of the first i tokens.
  cost <- c(0, rep(Inf, n))
  arg <- integer(n + 1)
  for (i in seq_len(n)) {
    best <- Inf
    bestJ <- 0L
    for (j in 0:(i - 1)) {                     # last line is tokens (j + 1)..i
      cand <- cost[j + 1] + lineCost(j + 1, i)
      if (cand < best) {
        best <- cand
        bestJ <- j
      }
    }
    cost[i + 1] <- best
    arg[i + 1] <- bestJ
  }

  # Backtrack to the line-start token indices.
  starts <- integer(0)
  i <- n
  while (i > 0) {
    j <- arg[i + 1]
    starts <- c(j + 1, starts)
    i <- j
  }
  ends <- c(starts[-1] - 1, n)
  joinLine <- function(a, b) {
    out <- toks[a]
    if (b > a) {
      for (t in (a + 1):b) {
        sep <- ' '
        if (endsHyphen[t - 1]) {
          sep <- ''
        }
        out <- paste0(out, sep, toks[t])
      }
    }
    out
  }
  paste(mapply(joinLine, starts, ends), collapse = '\n')
}
#' @importFrom grid unit is.unit convertWidth convertHeight
as_mm <- function(x, def, width = TRUE) {
  if (is.null(x)) return(x)
  if (!is.unit(x)) x <- unit(x, def)
  if (width) {
    convertWidth(x, 'mm', TRUE)
  } else {
    convertHeight(x, 'mm', TRUE)
  }
}

font_descent <- function(fontfamily, fontface, fontsize, cex) {
  info <- systemfonts::font_info(
    family = fontfamily,
    italic = isTRUE(fontface %in% c(3, 4)),
    weight = if (isTRUE(fontface %in% c(2, 4))) "bold" else "normal",
    size   = fontsize * (cex %||% 1), res = 300)
  as_mm(abs(info$max_descend) * 72 / 300, "pt", FALSE)

}

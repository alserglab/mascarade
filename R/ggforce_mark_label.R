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
                         description = NULL, width = NULL, min.width = 50,
                         default.units = 'mm', hjust = 0,
                         pad = margin(2, 2, 2, 2, 'mm'), gp = gpar(), desc.gp = gpar(),
                         vp = NULL) {
  width <- as_mm(width, default.units)
  min.width <- as_mm(min.width, default.units)
  pad <- as_mm(pad, default.units)
  pad[c(1, 3)] <- as_mm(pad[c(1, 3)], default.units, FALSE)
  if (!is.null(label) && !is.na(label)) {
    if (!is.null(width)) {
      label <- wrap_text(label, gp, width - pad[2] - pad[4])
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
    description <- wrap_text(description, desc.gp, final_width)
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
  italic <- fontface >= 3
  bold <- fontface == 2 | fontface == 4
  info <- systemfonts::font_info(fontfamily, italic, bold, fontsize * (cex %||% 1), res = 300)
  as_mm(abs(info$max_descend)*72/300, 'pt', FALSE)
}

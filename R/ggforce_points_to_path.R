# ChatGPT rewrite from C++ to R
# https://github.com/thomasp85/ggforce/blob/main/src/pointPath.cpp

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

# Squared Euclidean distance between two 2D points
dist_squared <- function(p, p1) {
    x <- p1[1] - p[1]
    y <- p1[2] - p[2]
    x * x + y * y
}

# Projection of point p onto segment ab (with optional clamping to the segment)
projection <- function(a, b, p, clamp = TRUE) {
    if (identical(a, b)) return(a)

    # vector from a to b
    ab <- b - a
    # vector from a to p
    ap <- p - a

    length2 <- sum(ab * ab)
    if (length2 == 0) return(a)

    t <- sum(ab * ap) / length2

    if (clamp) {
        t <- max(0, min(1, t))
    }

    a + t * ab
}

# For a single point (x, y), find closest point on a list of paths
# path: list of matrices, each with 2 columns (x, y)
dist_to_path <- function(x, y, path, closed_poly = FALSE) {
    shortest_dist <- -1
    closest <- c(NA_real_, NA_real_)
    point <- c(x, y)

    for (i in seq_along(path)) {
        mat <- path[[i]]
        n <- nrow(mat)
        if (n < 2L) next

        # segment indices
        if (closed_poly) {
            # segments: (1,2), (2,3), ..., (n-1,n), (n,1)
            idx_from <- seq_len(n)
            idx_to   <- c(2:n, 1L)
        } else {
            # segments: (1,2), (2,3), ..., (n-1,n)
            idx_from <- seq_len(n - 1L)
            idx_to   <- seq(2L, n)
        }

        for (k in seq_along(idx_from)) {
            a <- mat[idx_from[k], ]
            b <- mat[idx_to[k], ]

            close_pt <- projection(a, b, point, clamp = TRUE)
            dist_val <- sqrt(dist_squared(point, close_pt))

            if (shortest_dist < 0 || dist_val < shortest_dist) {
                shortest_dist <- dist_val
                closest <- close_pt
            }
        }
    }

    c(closest[1], closest[2], shortest_dist)
}

# Main function: R equivalent of points_to_path()
# pos  : numeric matrix with 2 columns (x, y)
# path : list of numeric matrices (each with 2 columns)
# close: logical, whether to treat each path as closed (polygon) or open (polyline)
points_to_path <- function(pos, path, close = FALSE) {
    pos <- as.matrix(pos)
    if (ncol(pos) != 2L) {
        stop("'pos' must be a matrix with 2 columns (x, y)")
    }

    n <- nrow(pos)
    proj <- matrix(NA_real_, nrow = n, ncol = 2L)
    dist <- numeric(n)

    for (i in seq_len(n)) {
        res <- dist_to_path(pos[i, 1], pos[i, 2], path, closed_poly = close)
        proj[i, ] <- res[1:2]
        dist[i]   <- res[3]
    }

    list(
        projection = proj,
        distance  = dist
    )
}

simplify_outer <- function(poly, max_area, min_vertices = 4L) {
  x <- poly$x
  y <- poly$y
  n <- length(x)
  if (n <= min_vertices || max_area <= 0) return(poly)

  # Circular linked list: prv[i] / nxt[i] are the neighbour indices of vertex i
  prv <- c(n, seq_len(n - 1L))
  nxt <- c(seq(2L, n), 1L)

  # Signed polygon area via shoelace: positive = CCW, negative = CW
  ccw <- sum(x * y[nxt] - x[nxt] * y) / 2 > 0

  # Vectorised initial cross products: (P_i - P_prev) x (P_next - P_prev)
  # Positive = left turn (convex in CCW), Negative = right turn (concave in CCW).
  crosses <- (x - x[prv]) * (y[nxt] - y[prv]) - (y - y[prv]) * (x[nxt] - x[prv])

  # For a CCW polygon a concave vertex (right turn, cross <= 0) dips inward.
  # Removing it fills the dent so the simplified polygon encloses the original.
  # Mark non-removable vertices with Inf so which.min() skips them.
  areas <- abs(crosses) / 2
  # Mark non-removable vertices with NA so which.min() skips them;
  # dead vertices (removed) are marked Inf.
  if (ccw) areas[crosses > 0] <- NA_real_ else areas[crosses < 0] <- NA_real_

  alive <- rep(TRUE, n)
  n_alive <- n
  repeat {
    if (n_alive <= min_vertices) break
    best_i <- which.min(areas)  # NA values are ignored by which.min
    if (length(best_i) == 0L || areas[best_i] > max_area) break

    p <- prv[best_i]; nx <- nxt[best_i]
    nxt[p] <- nx; prv[nx] <- p
    alive[best_i] <- FALSE
    areas[best_i] <- NA_real_     # exclude from future which.min (dead vertex)
    n_alive <- n_alive - 1L

    # Recompute cross products only for the two affected neighbours
    for (nb in c(p, nx)) {
      a_nb <- prv[nb]; b_nb <- nxt[nb]
      cr <- (x[nb] - x[a_nb]) * (y[b_nb] - y[a_nb]) -
            (y[nb] - y[a_nb]) * (x[b_nb] - x[a_nb])
      removable <- if (ccw) cr <= 0 else cr >= 0
      areas[nb] <- if (removable) abs(cr) / 2 else NA_real_
    }
  }

  cur <- which(alive)[1L]
  rx <- numeric(n_alive); ry <- numeric(n_alive)
  for (j in seq_len(n_alive)) { rx[j] <- x[cur]; ry[j] <- y[cur]; cur <- nxt[cur] }
  list(x = rx, y = ry)
}

#' @importFrom polyclip polyoffset polyminkowski polyclip
#' @importFrom grid convertX convertY
my_place_labels <- function(rects, polygons, bounds, anchors, ghosts,
                            simp_ratio = 0.001) {
  res <- vector('list', length(rects))
  bbox <- list(
    x = c(0, bounds[1], bounds[1], 0),
    y = c(0, 0, bounds[2], bounds[2])
  )
  if (!is.null(ghosts) && length(ghosts$x) > 0) {
    ghosts$x <- convertX(ghosts$x, 'mm', TRUE)
    ghosts$y <- convertY(ghosts$y, 'mm', TRUE)
    ghosts <- Map(
      function(xmin, xmax, ymin, ymax) {
        list(x = c(xmin, xmax, xmax, xmin), y = c(ymin, ymin, ymax, ymax))
      },
      xmin = ghosts$x - 2,
      xmax = ghosts$x + 2,
      ymin = ghosts$y - 2,
      ymax = ghosts$y + 2
    )
    ghosts <- polyoffset(ghosts, 0)
    polygons <- c(polygons, ghosts)
  }
  if (simp_ratio > 0 && length(polygons) > 0) {
    all_x <- unlist(lapply(polygons, `[[`, 'x'))
    all_y <- unlist(lapply(polygons, `[[`, 'y'))
    bbox_area <- diff(range(all_x)) * diff(range(all_y))
    max_area  <- simp_ratio * bbox_area
    polygons  <- lapply(polygons, simplify_outer, max_area = max_area)
  }
  for (i in seq_along(rects)) {
    if (all(rects[[i]] == 0)) next()
    r <- rects[[i]] / 2 + 2
    rect <- list(x = c(-r[1], r[1], r[1], -r[1]),
                 y = c(-r[2], -r[2], r[2], r[2]))
    b <- polyminkowski(bbox, rect)
    for (p in polygons) {
      b <- polyclip(b, polyminkowski(p, rect)[1], 'union')
    }
    if (length(b) == 1) next()
    b <- lapply(b[-1], function(p) cbind(p$x, p$y))
    closest <- points_to_path(matrix(anchors[[i]], ncol = 2), b, TRUE)
    res[[i]] <- closest$proj
    rect$x <- rect$x + closest$proj[1]
    rect$y <- rect$y + closest$proj[2]
    polygons[[length(polygons) + 1]] <- rect
  }
  res
}
#' @importFrom polyclip polyoffset
#' @importFrom grid convertWidth convertHeight nullGrob polylineGrob
#' @importFrom stats runif
my_make_label <- function(labels, dims, polygons, ghosts, buffer, con_type,
                       con_border, con_cap, con_gp, anchor_mod, anchor_x,
                       anchor_y, arrow, simp_ratio = 0.001) {
  polygons <- lapply(polygons, function(p) {
    if (length(p$x) == 1 & length(p$y) == 1) {
      list(
        x = runif(200, p$x-0.00005, p$x+0.00005),
        y = runif(200, p$y-0.00005, p$y+0.00005)
      )
    } else {
      list(
        x = p$x,
        y = p$y
      )
    }
  })

  anchors <- lapply(seq_along(polygons), function(i) {
    x <- mean(range(polygons[[i]]$x))
    if (length(anchor_x) == length(polygons) && !is.na(anchor_x[i])) x <- anchor_x[i]
    y <- mean(range(polygons[[i]]$y))
    if (length(anchor_y) == length(polygons) && !is.na(anchor_y[i])) y <- anchor_y[i]
    c(x, y)
  })
  p_big <- polyoffset(polygons, convertWidth(buffer, 'mm', TRUE))

  area <- c(
    convertWidth(unit(1, 'npc'), 'mm', TRUE),
    convertHeight(unit(1, 'npc'), 'mm', TRUE)
  )
  labelpos <- my_place_labels(dims, p_big, area, anchors, ghosts,
                              simp_ratio = simp_ratio)
  if (all(lengths(labelpos) == 0)) {
    return(list(nullGrob()))
  }
  labels_drawn <- which(!vapply(labelpos, is.null, logical(1)))
  labels <- Map(function(lab, pos) {
    if (is.null(pos) || inherits(lab, 'null')) return(nullGrob())
    lab$vp$x <- unit(pos[1], 'mm')
    lab$vp$y <- unit(pos[2], 'mm')
    lab
  }, lab = labels, pos = labelpos)
  connect <- rlang::inject(rbind(!!!Map(function(pol, pos, dim) {
    if (is.null(pos)) return(NULL)
    dim <- dim / anchor_mod
    pos <- cbind(
      c(pos[1] - dim[1], pos[1] + dim[1], pos[1] + dim[1], pos[1] - dim[1]),
      c(pos[2] - dim[2], pos[2] - dim[2], pos[2] + dim[2], pos[2] + dim[2])
    )
    pos <- points_to_path(pos, list(cbind(pol$x, pol$y)), TRUE)
    pos$projection[which.min(pos$distance), ]
  }, pol = polygons, pos = labelpos, dim = dims)))
  labeldims <- rlang::inject(rbind(!!!dims[lengths(labelpos) != 0])) / 2
  labelpos <- rlang::inject(rbind(!!!labelpos))
  if (con_type == 'none' || !con_type %in% c('elbow', 'straight')) {
    connect <- nullGrob()
  } else {
    con_fun <- switch(con_type, elbow = elbow, straight = straight)
    connect <- con_fun(
      labelpos[, 1] - labeldims[, 1], labelpos[, 1] + labeldims[, 1],
      labelpos[, 2] - labeldims[, 2], labelpos[, 2] + labeldims[, 2],
      connect[, 1], connect[, 2]
    )
    if (con_border == 'one') {
      connect <- with_borderline(
        labelpos[, 1] - labeldims[, 1],
        labelpos[, 1] + labeldims[, 1], connect
      )
    }
    connect <- end_cap(connect, con_cap)
    connect <- zip_points(connect)
    if (!is.null(arrow)) arrow$ends <- 2L
    con_gp <- subset_gp(con_gp, labels_drawn)
    connect <- polylineGrob(connect$x, connect$y,
      id = connect$id,
      default.units = 'mm', gp = con_gp, arrow = arrow
    )
  }
  c(labels, list(connect))
}


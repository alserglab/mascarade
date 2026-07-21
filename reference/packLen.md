# Fixed-order 1-D packing minimising total leader length

For labels in a fixed top-to-bottom order, chooses stacked centre
y-positions on a fine grid that minimise the total leader length
`sqrt(dx^2 + (cy - py)^2)`, subject to a minimum centre separation of
`(h_i + h_j)/2 + gap` between neighbours. Solved exactly for the given
order by a grid dynamic program (`O(n * gridSize)`): sweeping labels
from the bottom up, `g[t]` is the least total length with the current
label at some slot at or below `t` and every lower label packed beneath
it, via a two-option transition per slot – either place the label in
slot `t` (on top of the best lower label position) or skip the slot and
place it lower down. Paired with the Hungarian order (see the seed
reorder) this reproduces the assignment optimum – exactly in the
equal-height, single-line case – because it minimises the true Euclidean
leader length rather than a squared-vertical proxy.

## Usage

``` r
packLen(dx, py, h, gap, slot, ylo, yhi)
```

## Arguments

- dx:

  Numeric horizontal pole-to-column distances, in the fixed
  top-to-bottom order.

- py:

  Numeric pole y-coordinates, in the same order.

- h:

  Numeric full box heights along the stacking axis, in the same order.

- gap:

  Numeric extra separation added between neighbouring boxes.

- slot:

  Numeric grid resolution (slot height) for the candidate centre
  positions.

- ylo, yhi:

  Numeric viewport y-bounds the grid must cover.

## Value

Numeric vector of placed centre y-positions, aligned with the inputs.

## Details

The grid spans at least the whole viewport `[ylo, yhi]` (and the pole
range), so labels may use the full vertical space; when the stacked
column is taller than that space the grid is extended symmetrically
beyond it so every box still fits.

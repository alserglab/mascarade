# Single-move conflict sweep

Multi-pass Gauss-Seidel refinement. Each pass reorders labels by their
current conflict vector (box-box, leader-leader, leader-box, then leader
length) descending, then moves each to its lexicographically best
candidate versus the current others (including staying put). Stops when
a pass makes no change.

## Usage

``` r
oneMoveSweepKernel(
  cxmin,
  cxmax,
  cymin,
  cymax,
  ex,
  ey,
  tx,
  ty,
  len,
  rows,
  init,
  maxpass
)
```

## Arguments

- cxmin, cxmax, cymin, cymax:

  Numeric padded box extents, one per candidate.

- ex, ey:

  Numeric leader start (anchor) per candidate.

- tx, ty:

  Numeric leader target (pole) per candidate.

- len:

  Numeric ranking length per candidate.

- rows:

  List of integer vectors: the 0-indexed candidate rows available to
  each label.

- init:

  Integer vector: the starting candidate index per label.

- maxpass:

  Integer cap on the number of sweeps.

## Value

Integer vector: the chosen candidate index per label.

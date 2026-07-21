# Effective length used to rank label candidates

Per candidate, the length the optimizer minimises after the conflict
counts: the base leader length, plus the arc of the leader inside any
FOREIGN cluster (routes leaders around clusters), plus `OVERFLOW_WEIGHT`
times how far the box overflows the viewport (steers labels in-bounds).

## Usage

``` r
effectiveLength(
  len,
  ex,
  ey,
  tx,
  ty,
  lab,
  polysx,
  polysy,
  cxmin,
  cxmax,
  cymin,
  cymax,
  xlo,
  xhi,
  ylo,
  yhi
)
```

## Arguments

- len:

  Numeric base leader length per candidate.

- ex, ey:

  Numeric leader start (anchor) per candidate.

- tx, ty:

  Numeric leader target (pole) per candidate.

- lab:

  Integer 1-indexed own-cluster of each candidate.

- polysx, polysy:

  Lists of parallel numeric x/y vectors, one ring per cluster.

- cxmin, cxmax, cymin, cymax:

  Numeric padded box extents per candidate.

- xlo, xhi, ylo, yhi:

  Numeric viewport bounds (already inset by label.buffer).

## Value

Numeric vector of effective lengths, one per candidate.

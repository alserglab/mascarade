# Visible leader end (first mask-boundary hit)

For each label, the leader runs from its box anchor (`ex`, `ey`) to the
cluster pole (`tx`, `ty`), which lies inside the cluster. The drawn
leader stops at the mask boundary: the first point where the
anchor-\>pole segment crosses the label's own cluster ring. When the
segment does not cross the ring the leader runs all the way to the pole.

## Usage

``` r
firstLeaderHit(ex, ey, tx, ty, polysx, polysy)
```

## Arguments

- ex, ey:

  Numeric leader-start (anchor) coordinates, one per label.

- tx, ty:

  Numeric pole (leader target) coordinates, one per label.

- polysx, polysy:

  Lists of parallel numeric x/y vectors, one ring per label (its own
  cluster).

## Value

A list with numeric `bx`, `by`: the visible leader end, one per label.

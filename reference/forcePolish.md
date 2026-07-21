# Continuous force-directed label polish

Pattern-search descent on the (squared) EFFECTIVE length under a hard
conflict guard. The effective length of a label is its centre-to-pole
leader length, plus the arc of the leader that runs inside any foreign
cluster (routing leaders around clusters), plus a soft viewport overflow
penalty – the same quantity the upstream
[`effectiveLength()`](https://alserglab.github.io/mascarade/reference/effectiveLength.md)
ranking minimises. A box-box spacing penalty is added on top. Starting
from a conflict-free layout the search only accepts conflict-free
neighbours (free-space check via the BoxFit R-tree), so feasibility is
preserved. Overflow is a SOFT term folded into the effective length, not
a hard clip, so an off-panel label can be walked back in-bounds and a
label leaves the panel only when that lowers total energy. The leader
anchor rule mirrors R's `.anchorPoint()`, so both the conflict guard and
the foreign-cluster arc match the drawn geometry.

## Usage

``` r
forcePolish(
  boxfit,
  cx0,
  cy0,
  hw,
  hh,
  tx,
  ty,
  polysx,
  polysy,
  pad,
  xlo,
  xhi,
  ylo,
  yhi,
  iters,
  step,
  MU,
  pad_tgt,
  stepmin,
  con_type,
  sq = TRUE
)
```

## Arguments

- boxfit:

  External pointer from
  [`buildBoxFit()`](https://alserglab.github.io/mascarade/reference/buildBoxFit.md)
  (the cluster keep-out).

- cx0, cy0:

  Numeric starting label-centre coordinates.

- hw, hh:

  Numeric per-label box half-sizes.

- tx, ty:

  Numeric per-label pole (leader target).

- polysx, polysy:

  Lists of parallel numeric x/y vectors, one true mask ring per cluster
  (aligned with the labels), used for the foreign-cluster arc term.

- pad:

  Numeric hard box clearance.

- xlo, xhi, ylo, yhi:

  Numeric viewport bounds.

- iters:

  Integer iteration count.

- step:

  Numeric initial pattern-search step.

- MU:

  Numeric weight of the box-box spacing penalty.

- pad_tgt:

  Numeric target inter-box spacing.

- stepmin:

  Numeric smallest step tried before abandoning a direction.

- con_type:

  Integer leader style: 0 = "ledge" (corner), otherwise
  "line"/"box"/"none".

- sq:

  Logical; if `TRUE`, the length term uses the squared distance.

## Value

A list with numeric `cx`, `cy`: the polished label centres.

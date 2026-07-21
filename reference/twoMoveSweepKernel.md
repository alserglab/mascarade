# Two-move conflict / length refinement

Per label (heaviest leader first), picks the lexicographically best move
(by change in box-box, leader-leader, leader-box conflicts, then change
in length) and applies it if it improves. A CONFLICT-FREE label uses a
length branch-and-bound: candidates are length-ascending and a shorter
c1 needs a partner only when exactly one other label currently sits in
that slot (that label is the partner – it need not itself be in
conflict), with the length bound pruning the length-increasing tail –
the exact per-step optimum. A CONFLICTED label drops the pruning and
searches ALL pairs of candidates (this label's against every other
label's), driving out conflicts present in the input; only conflicted
labels pay this cost.

## Usage

``` r
twoMoveSweepKernel(
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
  maxpass,
  sq
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

  Integer cap on the number of passes.

- sq:

  Logical; if `TRUE` the length objective uses the squared length.

## Value

Integer vector: the chosen candidate index per label.

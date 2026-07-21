# Radial free-space label candidates

From each cluster pole, marches `ndir` rays outward and emits a
candidate box centre wherever a box of the label's size (plus `pad`) is
cluster-free (a BoxFit R-tree query): at the near edge of every free
interval and every `intfill` along wide ones. Candidates may fall partly
outside the viewport – the effective-length overflow term ranks them, so
a crowded label can take a minimally-clipped edge slot instead of the
far seed.

## Usage

``` r
radialCandidates(
  boxfit,
  poi,
  hw,
  hh,
  pad,
  ndir,
  step,
  rstart,
  rmax,
  intfill,
  dedup
)
```

## Arguments

- boxfit:

  External pointer from
  [`buildBoxFit()`](https://alserglab.github.io/mascarade/reference/buildBoxFit.md).

- poi:

  K x 2 matrix of cluster poles (the ray origins).

- hw, hh:

  Numeric per-label box half-sizes.

- pad:

  Numeric hard box clearance added around each box.

- ndir:

  Integer number of rays per pole.

- step:

  Numeric radial step along each ray.

- rstart, rmax:

  Numeric first and last radius searched.

- intfill:

  Numeric spacing of extra candidates along a wide free interval.

- dedup:

  Numeric grid size for de-duplicating nearby candidates.

## Value

A data.frame with integer `label` (1-indexed) and numeric `cx`, `cy`.

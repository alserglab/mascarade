# Hungarian (Jonker-Volgenant) assignment

O(n^3) minimum-cost assignment on a square cost matrix. Used by the
boundary seed (one solve per column) to match labels to stacked slots.

## Usage

``` r
hungarian(cost)
```

## Arguments

- cost:

  Square numeric cost matrix (`cost[i, j]` = cost of assigning row i to
  column j).

## Value

Integer vector `res` where `res[i]` is the 0-indexed column assigned to
row i, minimising the total cost.

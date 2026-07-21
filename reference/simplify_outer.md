# Enclosing polygon simplification

Greedily removes small concave (inward) vertices from a polygon: a
vertex is dropped when the triangle it cuts off has area below
`max_area`. Because only concave vertices are removed the simplified
polygon ENCLOSES the original, so the box-fit keep-out built from it
stays conservative. Used to cut vertex counts before placement.

## Usage

``` r
simplify_outer(poly, max_area, min_vertices = 4L)
```

## Arguments

- poly:

  A list with numeric `x`, `y` (the polygon vertices).

- max_area:

  Numeric area threshold; vertices whose cut-off triangle is smaller are
  removed.

- min_vertices:

  Integer floor on the number of vertices kept.

## Value

A list with simplified numeric `x`, `y`.

# Build the box-fit R-tree for the cluster polygons

Constructs a Boost.Geometry R-tree over the cluster polygon envelopes
(and keeps the polygons) so the candidate and polish kernels can answer
"does this label box overlap any cluster?" quickly. Rebuilt once per
placement; the mask itself is not.

## Usage

``` r
buildBoxFit(polysx, polysy)
```

## Arguments

- polysx, polysy:

  Lists of parallel numeric x/y vectors, one ring per cluster.

## Value

An external-pointer (`XPtr<BoxFit>`) handle to the box-fit structure.

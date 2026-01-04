# Generates mask from a Seurat object. Requires `SeuratObject` package.

Generates mask from a Seurat object. Requires `SeuratObject` package.

## Usage

``` r
generateMaskSeurat(
  object,
  reduction = NULL,
  group.by = NULL,
  gridSize = 200,
  expand = 0.005,
  minSize = 10
)
```

## Arguments

- object:

  Seurat object

- reduction:

  character vector specifying which reduction to use (default:
  `DefaultDimReduc(object)`)

- group.by:

  character vector specifying which field to use for clusters (default:
  `"ident"`)

- gridSize:

  target width and height of the raster used internally

- expand:

  distance used to expand borders, represented as a fraction of
  sqrt(width\*height). Default: 1/200.

- minSize:

  Groups of less than `minSize` points are ignored, unless it is the
  only group for a cluster

## Value

data.table with points representing the mask borders. Each individual
border line corresponds to a single level of `group` column. Cluster
assignment is in `cluster` column.

## Examples

``` r
data("exampleSeurat")

maskTable <- generateMaskSeurat(exampleSeurat)

library(ggplot2)
Seurat::DimPlot(exampleSeurat) +
    geom_path(data=maskTable, aes(x=UMAP_1, y=UMAP_2, group=group)) +
    coord_fixed()
```

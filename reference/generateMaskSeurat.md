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
# only run if Seurat is installed
if (require("Seurat")) {
    data("pbmc_small")
    maskTable <- generateMaskSeurat(pbmc_small)

    library(ggplot2)
    # not the best plot, see vignettes for better examples
    DimPlot(pbmc_small) +
        geom_path(data=maskTable, aes(x=tSNE_1, y=tSNE_2, group=group))
}
#> Loading required package: Seurat
#> Loading required package: SeuratObject
#> Loading required package: sp
#> 
#> Attaching package: ‘SeuratObject’
#> The following objects are masked from ‘package:base’:
#> 
#>     intersect, t
```

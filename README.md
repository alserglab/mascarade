
[![R-CMD-check](https://github.com/alserglab/mascarade/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/alserglab/mascarade/actions/workflows/R-CMD-check.yaml)

# mascarade

`mascarade` package implements a procedure to automatically generate 2D
masks for clusters on single-cell dimensional reduction plots like t-SNE
or UMAP.

See the
[tutorial](https://alserglab.github.io/mascarade/articles/mascarade-tutorial.html)
for usage details and
[gallery](https://alserglab.github.io/mascarade/articles/mascarade-gallery.html)
for examples on different datasets.

## Installation

The package can be installed from GitHub:

``` r
remotes::install_github("alserglab/mascarade")
```

## Quick run

Loading necessary libraries:

``` r
library(mascarade)
library(ggplot2)
library(data.table)
library(Seurat)
```

Loading get the example PBMC3K dataset:

``` r
pbmc3k <- readRDS(url("https://alserglab.wustl.edu/files/mascarade/examples/pbmc3k_seurat5.rds"))
pbmc3k <- NormalizeData(pbmc3k)
pbmc3k
```

    ## An object of class Seurat 
    ## 13714 features across 2638 samples within 1 assay 
    ## Active assay: RNA (13714 features, 2000 variable features)
    ##  2 layers present: counts, data
    ##  2 dimensional reductions calculated: pca, umap

Generating masks:

``` r
maskTable <- generateMaskSeurat(pbmc3k)
```

`DimPlot` with the mask and labels:

``` r
DimPlot(pbmc3k) + NoLegend() +
    fancyMask(maskTable, ratio=1)
```

<img src="https://alserglab.github.io/mascarade/articles/mascarade-tutorial_files/figure-html/seurat-dimplot-1.png">

`DimPlot` with the just the labels

``` r
DimPlot(pbmc3k) + NoLegend() +
    fancyMask(maskTable, linewidth=0, ratio=1)
```

<img src="https://alserglab.github.io/mascarade/articles/mascarade-tutorial_files/figure-html/seurat-dimplot-noborder-1.png">

`FeaturePlot` with the mask and labels showing GNLY gene being specific
to NK cells:

``` r
FeaturePlot(pbmc3k, "GNLY", cols=c("grey90", "red")) +
    fancyMask(maskTable, ratio=1)
```

<img src="https://alserglab.github.io/mascarade/articles/mascarade-tutorial_files/figure-html/seurat-gnly-1.png">

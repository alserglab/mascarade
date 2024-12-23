
# mascarade

`mascarade` package implements a procedure to automatically generate 2D
masks for clusters on single-cell dimensional reduction plots like t-SNE
or UMAP.

## Installation

The package can be installed from GitHub:

``` r
remotes::install_github("alserglab/mascarade")
```

## Quick run

Here is quick example. See the
[vignette](https://rpubs.com/asergushichev/mascarade-tutorial) for more
details.

Loading neccessary libraries:

``` r
library(mascarade)
library(ggplot2)
library(data.table)
```

Loading example data:

``` r
data("exampleMascarade")
```

Generating masks:

``` r
maskTable <- generateMask(dims=exampleMascarade$dims, 
                          clusters=exampleMascarade$clusters)
```

Plotting with `ggplot2`:

``` r
data <- data.table(exampleMascarade$dims, 
                   cluster=exampleMascarade$clusters,
                   exampleMascarade$features)

ggplot(data, aes(x=UMAP_1, y=UMAP_2)) + 
    geom_point(aes(color=cluster)) + 
    geom_path(data=maskTable, aes(group=group)) +
    coord_fixed() + 
    theme_classic()
```

<img src="https://alserglab.wustl.edu/files/mascarade/plot.png">


# mascarade

`mascarade` package implements a procedure to automatically generate 2D
masks for clusters on single-cell dimensional reduction plots like t-SNE
or UMAP.

See the [tutorial](https://rpubs.com/asergushichev/mascarade-tutorial)
for usage details and
[gallery](https://rpubs.com/asergushichev/mascarade-gallery) for
examples on different datasets.

## Installation

The package can be installed from GitHub:

``` r
remotes::install_github("alserglab/mascarade")
```

## Quick run

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

<img src="https://alserglab.wustl.edu/files/mascarade/readme-basic-1.png">

Fancy version, showing NGLY gene being specific to NK cells:

``` r
library(ggforce)
library(ggnewscale)
fancyMask <- list(
    ggforce::geom_shape(data=maskTable, aes(group=group, color=cluster),
               linewidth=1, fill=NA, expand=unit(-1, "pt"), show.legend = FALSE),
    ggforce::geom_mark_hull(data=maskTable, fill = NA, aes(group=cluster, color=cluster, label = cluster),
                   linewidth=0,
                   radius=0, expand=0, con.cap=0, con.type = "straight",
                   label.fontsize = 10, label.buffer = unit(0, "cm"),
                   label.fontface = "plain",
                   label.minwidth = 0,
                   label.margin = margin(2, 2, 2, 2, "pt"),
                   label.lineheight = 0,
                   con.colour = "inherit",
                   show.legend = FALSE),
    # expanding to give a bit more space for labels
    scale_x_continuous(expand = expansion(mult = 0.1)),
    scale_y_continuous(expand = expansion(mult = 0.1))
)
ggplot(data, aes(x=UMAP_1, y=UMAP_2)) + 
    geom_point(aes(color=GNLY), size=0.5) +
    scale_color_gradient2(low = "#404040", high="red") + 
    new_scale_color() + 
    fancyMask +
    coord_fixed() + 
    theme_classic()
```

<img src="https://alserglab.wustl.edu/files/mascarade/readme-fancy-1.png">

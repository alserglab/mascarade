# Gallery of mascarade-generaded masks

### Loading necessary libraries

``` r

library(mascarade)
library(data.table)
library(ggplot2)
library(colorrepel)
```

### PBMC-3K UMAP

``` r

example <- readRDS(url("https://alserglab.wustl.edu/files/mascarade/examples/pbmc3k_umap.rds"))
data <- data.table(example$dims, 
                   cluster=example$clusters)

maskTable <- generateMask(dims=example$dims, 
                          clusters=example$clusters)

ggplot(data, aes(x=UMAP_1, y=UMAP_2)) + 
    geom_point(aes(color=cluster)) + 
    geom_path(data=maskTable, aes(group=group)) +
    coord_fixed() + 
    theme_classic()
```

![](mascarade-gallery_files/figure-html/unnamed-chunk-2-1.png)

With labels:

``` r

ggplot(data, aes(x=UMAP_1, y=UMAP_2)) + 
    geom_point(aes(color=cluster)) + 
    fancyMask(maskTable, ratio=1, linewidth = 0) +
    theme_classic() + theme(legend.position = "none")
```

![](mascarade-gallery_files/figure-html/unnamed-chunk-3-1.png)

### PBMC-3K t-SNE

``` r

example <- readRDS(url("https://alserglab.wustl.edu/files/mascarade/examples/pbmc3k_tsne.rds"))
data <- data.table(example$dims, 
                   cluster=example$clusters)

maskTable <- generateMask(dims=example$dims, 
                          clusters=example$clusters)

ggplot(data, aes(x=tSNE_1, y=tSNE_2)) + 
    geom_point(aes(color=cluster)) + 
    geom_path(data=maskTable, aes(group=group)) +
    coord_fixed() + 
    theme_classic()
```

![](mascarade-gallery_files/figure-html/unnamed-chunk-4-1.png)

With labels:

``` r

ggplot(data, aes(x=tSNE_1, y=tSNE_2)) + 
    geom_point(aes(color=cluster)) + 
    fancyMask(maskTable, ratio=1, linewidth = 0) +
    theme_classic() + theme(legend.position = "none")
```

![](mascarade-gallery_files/figure-html/unnamed-chunk-5-1.png)

### Aya

``` r

example <- readRDS(url("https://alserglab.wustl.edu/files/mascarade/examples/aya.rds"))
data <- data.table(example$dims, 
                   cluster=example$clusters)

maskTable <- generateMask(dims=example$dims, 
                          clusters=example$clusters)

ggplot(data, aes(x=UMAP_1, y=UMAP_2)) + 
    geom_point(aes(color=cluster), size=0.5) + 
    scale_color_repel() +
    geom_path(data=maskTable, aes(group=group)) +
    coord_fixed() + 
    theme_classic()
```

![](mascarade-gallery_files/figure-html/unnamed-chunk-6-1.png)

### Chia-Jung

``` r

example <- readRDS(url("https://alserglab.wustl.edu/files/mascarade/examples/chiajung1.rds"))
data <- data.table(example$dims, 
                   cluster=example$clusters)

maskTable <- generateMask(dims=example$dims, 
                          clusters=example$clusters)

ggplot(data, aes(x=UMAP_1, y=UMAP_2)) + 
    geom_point(aes(color=cluster), size=0.1) + 
    scale_color_repel() +
    geom_path(data=maskTable, aes(group=group)) +
    coord_fixed() + 
    theme_classic()
```

![](mascarade-gallery_files/figure-html/unnamed-chunk-7-1.png)

``` r

example <- readRDS(url("https://alserglab.wustl.edu/files/mascarade/examples/chiajung2.rds"))
data <- data.table(example$dims, 
                   cluster=example$clusters)

maskTable <- generateMask(dims=example$dims, 
                          clusters=example$clusters)

ggplot(data, aes(x=UMAP_1, y=UMAP_2)) + 
    geom_point(aes(color=cluster)) + 
    scale_color_repel() +
    geom_path(data=maskTable, aes(group=group)) +
    coord_fixed() + 
    theme_classic()
```

![](mascarade-gallery_files/figure-html/unnamed-chunk-8-1.png)

### Vladimir Shitov

Really hard case, so playing with parameters a bit to make more space.

``` r

example <- readRDS(url("https://alserglab.wustl.edu/files/mascarade/examples/vshitov.rds"))
data <- data.table(example$dims, 
                   cluster=example$clusters)

maskTable <- generateMask(dims=example$dims, 
                          clusters=example$clusters)

ggplot(data, aes(x=UMAP1, y=UMAP2)) + 
    geom_point(aes(color=cluster), size=0.1) + 
    scale_color_repel() +
    fancyMask(maskTable, 
              con.type = "line", # ledges clutter the image visually here
              ratio=1, 
              linewidth = 0, 
              limits.expand = c(0.2, 0.1), # more space to the left and right of the plot
              label.buffer = unit(1, "mm"), # smaller buffer zone around clusters
              label.width = unit(25, "mm"), # wrap long lines
              label.margin = margin(1, 1, 1, 1, "pt"), # tighter label margins
              label.fontsize = 10) +
    theme_classic() + theme(legend.position = "none")
```

![](mascarade-gallery_files/figure-html/unnamed-chunk-9-1.png)

### Session info

``` r

sessionInfo()
```

    ## R version 4.6.1 (2026-06-24)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: Ubuntu 24.04.4 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
    ## LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
    ##  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
    ##  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
    ## [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
    ## 
    ## time zone: UTC
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] colorrepel_0.5.0  ggplot2_4.0.3     data.table_1.18.4 mascarade_0.4.0  
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] dqrng_0.4.1            sass_0.4.10            generics_0.1.4        
    ##  [4] spatstat.explore_3.8-1 gtools_3.9.5           polylabelr_1.0.0      
    ##  [7] tensor_1.5.1           distances_0.1.13       spatstat.data_3.1-9   
    ## [10] lattice_0.22-9         digest_0.6.39          magrittr_2.0.5        
    ## [13] spatstat.utils_3.2-4   evaluate_1.0.5         grid_4.6.1            
    ## [16] RColorBrewer_1.1-3     fastmap_1.2.0          jsonlite_2.0.0        
    ## [19] Matrix_1.7-5           spatstat.sparse_3.2-0  purrr_1.2.2           
    ## [22] scales_1.4.0           tweenr_2.0.3           textshaping_1.0.5     
    ## [25] jquerylib_0.1.4        abind_1.4-8            cli_3.6.6             
    ## [28] rlang_1.3.0            polyclip_1.10-7        withr_3.0.3           
    ## [31] cachem_1.1.0           yaml_2.3.12            otel_0.2.0            
    ## [34] spatstat.univar_3.2-0  tools_4.6.1            deldir_2.0-4          
    ## [37] dplyr_1.2.1            spatstat.geom_3.8-1    vctrs_0.7.3           
    ## [40] R6_2.6.1               matrixStats_1.5.0      lifecycle_1.0.5       
    ## [43] fs_2.1.0               htmlwidgets_1.6.4      MASS_7.3-65           
    ## [46] ragg_1.5.2             pkgconfig_2.0.3        desc_1.4.3            
    ## [49] pkgdown_2.2.1          pillar_1.11.1          bslib_0.11.0          
    ## [52] gtable_0.3.6           Rcpp_1.1.2             glue_1.8.1            
    ## [55] ggforce_0.5.0          systemfonts_1.3.2      xfun_0.60             
    ## [58] tibble_3.3.1           tidyselect_1.2.1       knitr_1.51            
    ## [61] goftest_1.2-3          farver_2.1.2           nlme_3.1-169          
    ## [64] spatstat.random_3.5-0  htmltools_0.5.9        labeling_0.4.3        
    ## [67] rmarkdown_2.31         compiler_4.6.1         S7_0.2.2

library(data.table)
library(patchwork)
library(ggnewscale)
library(mascarade)

data <- data.table(exampleMascarade$dims,
                   cluster=as.factor(exampleMascarade$clusters),
                   exampleMascarade$features)

maskTable <- generateMask(dims=exampleMascarade$dims,
                          clusters=as.factor(exampleMascarade$clusters))


p1 <- ggplot(data, aes(x=UMAP_1, y=UMAP_2)) +
    geom_point(aes(color=cluster), size=0.25) +
    coord_fixed() +
    theme_classic()

p1


p2 <- ggplot(data, aes(x=UMAP_1, y=UMAP_2)) +
    geom_point(aes(color=GNLY), size=0.25) +
    scale_color_gradient2(low = "#404040", high="red") +
    new_scale_color() +
    geom_path(data=maskTable, aes(group=group, color=cluster), linewidth=0.5) +
    scale_color_discrete(guide="none") +
    coord_fixed() +
    theme_classic()

p <- p1 + p2
p

ggsave(p, file="combined_plot.png", width=8, height=3)

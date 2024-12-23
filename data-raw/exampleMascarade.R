## code to prepare `exampleMascarade` dataset goes here
library(SeuratData)
library(Seurat)

InstallData("pbmc3k")
LoadData("pbmc3k")

pbmc3k.final <- Seurat::UpdateSeuratObject(pbmc3k.final)

featureList <- c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP",
                 "CD8A")

exampleMascarade <- list(
    dims=Embeddings(pbmc3k.final, "umap"),
    clusters=pbmc3k.final$seurat_annotations,
    features=t(pbmc3k.final[["RNA"]]@scale.data[featureList, ])
)

usethis::use_data(exampleMascarade, overwrite = TRUE)

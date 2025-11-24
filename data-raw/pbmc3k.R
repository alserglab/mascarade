library(Seurat)
library(SeuratData)

if (!AvailableData()["pbmc3k", "Installed"]) {
    InstallData("pbmc3k")
}

LoadData("pbmc3k")

pbmc3k.final <- UpdateSeuratObject(pbmc3k.final)
pbmc3k <- DietSeurat(pbmc3k.final, layers = c("counts"), dimreducs = c("pca", "umap"))
saveRDS(pbmc3k, file="pbmc3k_seurat5.rds")


#' Generates mask from a Seurat object. Requires `SeuratObject` package.
#'
#' @param object Seurat object
#' @param reduction character vector specifying which reduction to use
#'      (default: `DefaultDimReduc(object)`)
#' @param group.by character vector specifying which field to use for clusters
#'      (default: `"ident"`)
#' @param name description
#' @returns data.table with points representing the mask borders.
#'      Each individual border line corresponds to a single level of `group` column.
#'      Cluster assignment is in `cluster` column.
#' @inheritParams generateMask
#' @export
#' @examples
#' data("exampleSeurat")
#'
#' maskTable <- generateMaskSeurat(exampleSeurat)
#'
#' Seurat::DimPlot(exampleSeurat) +
#'     geom_path(data=maskTable, aes(x=UMAP_1, y=UMAP_2, group=group)) +
#'     coord_fixed()
generateMaskSeurat <- function(object,
                               reduction = NULL,
                               group.by = NULL,
                               gridSize = 200,
                               expand = 0.005,
                               minSize = 10) {
    stopifnot(requireNamespace("SeuratObject", quietly = TRUE))

    if (is.null(reduction)) {
        reduction <- SeuratObject::DefaultDimReduc(object)
    }
    dims <- SeuratObject::Embeddings(object, reduction=reduction)

    if (is.null(group.by)) {
        group.by <- "ident"
    }
    clusters <- SeuratObject::FetchData(object, group.by)[,1]

    maskTable <- generateMask(dims = dims,
                              clusters = clusters,
                              gridSize = gridSize,
                              expand = expand,
                              minSize = minSize)

}

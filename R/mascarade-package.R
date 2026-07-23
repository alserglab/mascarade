#' Example data with UMAP points from PBMC3K dataset.
#'
#' The object is a list with three elements:
#' 1) `dims` -- matrix of UMAP coordinates of the cells,
#' 2) `clusters` -- vector of cell population annotations,
#' 3) `features` -- matrix withgene expression for several genes.
#'
#' @docType data
#' @name exampleMascarade
NULL

#' Example generated mask.
#'
#' Result of `generateMasks(dims=exampleMascarade$dims, clusters=exampleMascarade$clusters)`.
#' @docType data
#' @name exampleMascarade
NULL

#' @useDynLib mascarade, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generateMask.R
\name{generateMask}
\alias{generateMask}
\title{Generate mask for clusters on 2D dimensional reduction plots}
\usage{
generateMask(
  dims,
  clusters,
  gridSize = 200,
  minDensity = 0.1,
  smoothSigma = 0.025,
  type = c("partition", "independent")
)
}
\arguments{
\item{dims}{matrix of point coordinates.
Rows are points, columns are dimensions. Only the first two columns are used.}

\item{clusters}{vector of cluster annotations.
Should be the same length as the number of rows in `dims`.}

\item{gridSize}{width and height of the raster used internally}

\item{minDensity}{minimal required density for the grid cells to be included in the mask.
Decreasing this parameter will expand masks.}

\item{smoothSigma}{sigma used in Gaussian smoothing represented as a fraction of plot width.
Increasing this parameter can help dealing with sparse regions.}

\item{type}{controls the behavior of the method.
When set to "partition" (default) generated masks are mutually exclusive.
When set to "independent" masks can overlap.}
}
\value{
data.table with points representing the mask borders.
     Each individual border line corresponds to a single level of `group` column.
     Cluster assignment is in `cluster` column.
}
\description{
Internally the function rasterizes and smoothes the density plots.
}
\examples{
data("exampleMascarade")
res <- generateMask(dims=exampleMascarade$dims,
                    clusters=exampleMascarade$clusters)
\dontrun{
data <- data.table(exampleMascarade$dims,
                   cluster=exampleMascarade$clusters,
                   exampleMascarade$features)
ggplot(data, aes(x=UMAP_1, y=UMAP_2)) +
    geom_point(aes(color=cluster)) +
    geom_path(data=maskTable, aes(group=group)) +
    coord_fixed() +
    theme_classic()
}
}

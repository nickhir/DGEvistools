% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualization_functions.R
\name{gwena_heatmap}
\alias{gwena_heatmap}
\title{Heatmap of counts after running GWENA}
\usage{
gwena_heatmap(
  scaled_counts,
  mdata,
  info_col,
  column_cluster = NULL,
  GWENA_modules,
  heatmap_annotation = NULL,
  ...
)
}
\arguments{
\item{scaled_counts}{A matrix of z-transformed counts. Columns should correspond to samples and rows to genes.}

\item{mdata}{A dataframe with information about the samples}

\item{info_col}{A column in the \code{mdata} object, which has entries that match the column names of \code{scaled_counts}.
They have to be in the same order, otherwise an error will occure.}

\item{column_cluster}{A named vector of subjects and corresponding cluster. Use only if you already determined how to cluster patients.
If none is specified, \code{hclust} is used to perform the clustering.}

\item{GWENA_modules}{Absolute path to the identified co-expression modules. Generated using \code{GWENA::detect_modules()} followed by \code{saveRDS()}.
The whole \code{GWENA} workflow has to be performed using the same subjects and genes that are in the \code{scaled_counts}.}

\item{heatmap_annotation}{A \href{https://jokergoo.github.io/ComplexHeatmap-reference/book/heatmap-annotations.html}{HeatmapAnnotation object}.}

\item{...}{Other arguments that are passed directly to the \code{ComplexHeatmap::Heatmap} function. Might be something like \code{show_row_names=T}}
}
\value{
A \code{ComplexHeatmap} which shows the GWENA results.
}
\description{
After running GWENA (WGCNA) this function can be used to visualize the results in a heatmap.
Here, the rows get split into the different modules and columns into either predefined columns clusters or they get clustered using \code{hclust}.
}

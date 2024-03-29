% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualization_functions.R
\name{module_association}
\alias{module_association}
\title{Association between eigengene expression and phenotype}
\usage{
module_association(
  mdata_eigengene,
  phenotype,
  module_eigengene,
  xlab = ggplot2::waiver(),
  ylab = ggplot2::waiver(),
  label_x_pos = ggplot2::waiver(),
  label_y_pos = ggplot2::waiver(),
  label_size = 7,
  ymin = -0.6,
  ymax = 0.6,
  comparison = "all",
  colors = NULL
)
}
\arguments{
\item{mdata_eigengene}{A dataframe which contains both the module eigengene expression and the phenotype
for which you want to evaluate the association}

\item{phenotype}{Column of \code{mdata_eigengene} which contains the phenotype of interest.
Make sure that continuous phenotypes are of class \code{numeric} or \code{integer} and that
categorical variables are of class \code{factor} or \code{characters}.}

\item{module_eigengene}{Column of \code{mdata_eigengene} which contains the eigengene expression of the module of interest}

\item{xlab}{The x-axis title}

\item{ylab}{The y-axis title}

\item{label_x_pos}{X-axis position to write the calculated correlation coefficient. Only relevant for continuous phenotypes.}

\item{label_y_pos}{Y-axis position to write the calculated correlation coefficient. Only relevant for continuous phenotypes.}

\item{label_size}{Font size of the correlation coefficient. Only relevant for continuous phenotypes.}

\item{ymin}{Specify minimum of the Y-axis range.}

\item{ymax}{Specify maximum of the Y-axis range.}

\item{comparison}{Specify the type of statistical comparison you want to perform for categorical phenotypes.
For examples how to specify comparisons, check the \href{https://cran.r-project.org/web/packages/ggsignif/vignettes/intro.html}{homepage} of \code{ggsignif}
or the vignette. By default, all factors will be compared. If no comparisons should be performed, specify \code{NULL}.}

\item{colors}{Specify the colors of the boxplots. Number of colors must match number of factors in \code{phenotype}.}
}
\description{
The function visualizes the assocation between a modules eigengene expression and a phenotype.
For continuous phenotypes (such as age) the correlation will be plotted. For categorical phenotypes (such as gender) a boxplot will be plotted
}
\examples{
\dontrun{
library(ggpubr)
modules <- readRDS(module_save_path) # generated using GWENA::detect_modules()
mdata <- fread("path/to/file/with/metadata")
#add the eigengene expression to the metadata
mdata_eigengene <- modules$modules_eigengenes \%>\%
  rownames_to_column("sample.ID") \%>\%
  left_join(., mdata, by="sample.ID")

module_association(mdata_eigengene = mdata_eigengene,
phenotype = "sample.status",
module_eigengene = "ME2",
ymin=-0.5, ymax=0.4,
xlab = "Sample status",
ylab= "Module 2 eigengene epxression",
colors=c("firebrick","steelblue"))
}



}

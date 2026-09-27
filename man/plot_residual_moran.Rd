\name{plot_residual_moran}
\alias{plot_residual_moran}
\title{Plot held-out residual Moran's I}
\usage{plot_residual_moran(object, distances = c(50000, 100000, 150000))}
\arguments{
\item{object}{A \code{bottom_up_cv} object.}
\item{distances}{Distance thresholds in model-coordinate units.}
}
\value{Moran's I table, invisibly.}
\description{Plots residual Moran's I across distance thresholds.}

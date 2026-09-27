\name{bottom_up_residual_moran}
\alias{bottom_up_residual_moran}
\title{Moran's I for held-out residuals}
\usage{bottom_up_residual_moran(object, distances = c(50000, 100000, 150000))}
\arguments{
\item{object}{A \code{bottom_up_cv} object.}
\item{distances}{Distance thresholds in the same units as model coordinates.}
}
\value{A data frame of Moran's I statistics by distance threshold.}
\description{Computes binary-distance-band Moran's I for geographically held-out residuals.}

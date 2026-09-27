\name{fit_bottom_up}
\alias{fit_bottom_up}
\title{Fit a bottom-UpR specification}
\usage{
fit_bottom_up(spec, mesh = NULL, prior_range = c(100, 0.05),
  prior_sigma = c(1, 0.05), config = TRUE, verbose = TRUE)
}
\arguments{
\item{spec}{A \code{bottom_up_spec} object.}
\item{mesh}{Optional INLA mesh.}
\item{prior_range}{PC prior specification for spatial range.}
\item{prior_sigma}{PC prior specification for spatial marginal SD.}
\item{config}{Store INLA configurations for joint posterior sampling.}
\item{verbose}{Logical; print model-fitting progress.}
}
\value{A \code{bottom_up_fit} object.}
\description{Fits the requested or automatically selected likelihood using INLA.}

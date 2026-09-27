\name{make_bottom_up_mesh}
\alias{make_bottom_up_mesh}
\title{Build an SPDE mesh}
\usage{
make_bottom_up_mesh(data, coords = c("x", "y"), adaptive = FALSE,
  max_edge = NULL, offset = NULL, cutoff = NULL,
  adaptive_args = list())
}
\arguments{
\item{data}{Data frame containing projected coordinates.}
\item{coords}{Two projected coordinate columns.}
\item{adaptive}{Derive unspecified mesh distances from the spatial footprint of \code{data}.}
\item{max_edge}{Inner and outer maximum triangle-edge lengths; \code{NULL} permits automatic derivation in adaptive mode.}
\item{offset}{Inner and outer mesh extensions; \code{NULL} permits automatic derivation in adaptive mode.}
\item{cutoff}{Minimum point separation; \code{NULL} permits automatic derivation in adaptive mode.}
\item{adaptive_args}{Named list passed to \code{bottom_up_mesh_parameters}.}
}
\value{An INLA mesh with the parameters used stored in an attribute.}
\description{Builds a two-dimensional SPDE mesh using explicit distances or study-area-adaptive distances derived from the geographic footprint of the supplied projected locations.}

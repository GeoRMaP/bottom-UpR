\name{make_bottom_up_mesh}
\alias{make_bottom_up_mesh}
\title{Build an SPDE mesh}
\usage{
make_bottom_up_mesh(data, coords = c("x", "y"),
  max_edge = c(50, 100), offset = c(50, 100), cutoff = 30)
}
\arguments{
\item{data}{Data frame containing coordinates.}
\item{coords}{Two projected coordinate columns.}
\item{max_edge}{Inner and outer maximum triangle-edge lengths.}
\item{offset}{Inner and outer mesh extensions.}
\item{cutoff}{Minimum point separation.}
}
\value{An INLA mesh object.}
\description{Builds a two-dimensional SPDE mesh. Distance arguments use the coordinate units.}

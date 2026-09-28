\name{simulate_bottom_up_data}
\alias{simulate_bottom_up_data}
\title{Simulate synthetic bottom-up population data}
\usage{
simulate_bottom_up_data(n = 200L, seed = 123,
  include_hierarchy = TRUE, extent = 300000,
  chunkwise = FALSE, chunk_size = 50000L, verbose = TRUE)
}
\arguments{
\item{n}{Number of synthetic enumeration areas.}
\item{seed}{Random seed.}
\item{include_hierarchy}{Include source, settlement, region and EA identifiers.}
\item{extent}{Spatial side length in coordinate units.}
\item{chunkwise}{Logical; generate each variable in chunks.}
\item{chunk_size}{Rows per simulation chunk when \code{chunkwise=TRUE}.}
\item{verbose}{Logical; print simulation progress.}
}
\value{A synthetic data frame containing mapped buildings, coordinates,
covariates, PPB and count responses, and optional hierarchical grouping variables.}
\description{Generates reproducible synthetic data for examples, tutorials and package tests. It is not intended to reproduce any specific country's data-generating process.}

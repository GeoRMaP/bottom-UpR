\name{select_likelihood}
\alias{select_likelihood}
\title{Select a compatible likelihood}
\usage{
select_likelihood(spec, mesh = NULL, prior_range = c(100, 0.05),
  prior_sigma = c(1, 0.05), config = FALSE, verbose = TRUE)
}
\arguments{
\item{spec}{A \code{bottom_up_spec}.}
\item{mesh}{Optional common SPDE mesh.}
\item{prior_range}{PC range-prior specification.}
\item{prior_sigma}{PC spatial-SD-prior specification.}
\item{config}{Enable joint-posterior configuration in candidate fits so the selected candidate can be reused without an additional refit.}
\item{verbose}{Logical; print likelihood-comparison progress.}
}
\value{A list with selected likelihood, candidate scores, fits and mesh.}
\description{Compares compatible likelihoods using mean log conditional predictive ordinate on the original response scale.}

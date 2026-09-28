\name{bottom_up_inla_selector}
\alias{bottom_up_inla_selector}
\title{Training-only INLA covariate selector}
\usage{
bottom_up_inla_selector(train_data, response, buildings = NULL,
  candidate_covariates, response_type, likelihood = "auto",
  coords = c("x", "y"), hierarchical_effects = NULL,
  spatial = TRUE, mesh = NULL, level = 0.95, fallback = 0L,
  verbose = TRUE, ...)
}
\arguments{
\item{train_data}{Training-fold data only.}
\item{response}{Response column.}
\item{buildings}{Optional building/exposure column.}
\item{candidate_covariates}{Candidate fixed effects.}
\item{response_type}{Count or PPB.}
\item{likelihood}{Automatic or explicit likelihood.}
\item{coords}{Projected coordinate columns.}
\item{hierarchical_effects}{Optional hierarchy.}
\item{spatial}{Include an SPDE field in the selection model.}
\item{mesh}{Optional pre-built mesh. Reusing a common mesh avoids repeated triangulation.}
\item{level}{Posterior credible level.}
\item{fallback}{Number of strongest signals retained if none exclude zero.}
\item{verbose}{Logical; print selection progress and fallback warnings.}
\item{...}{Additional fitting arguments.}
}
\value{Character vector of selected covariate names.}
\description{Optional generic INLA selector; users may supply another selector to \code{bottom_up_cv}.}

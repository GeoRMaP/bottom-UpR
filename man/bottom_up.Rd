\name{bottom_up}
\alias{bottom_up}
\title{Fit a bottom-up population model}
\usage{
bottom_up(data, response, buildings = NULL,
  response_type = c("auto", "count", "ppb"), likelihood = "auto",
  covariates = NULL, coords = c("x", "y"), random_effects = NULL,
  hierarchical_effects = NULL, spatial = TRUE, ...)
}
\arguments{
\item{data}{Modelling data frame.}
\item{response}{Response column name.}
\item{buildings}{Optional building/exposure column name.}
\item{response_type}{Response support: automatic, count, or PPB.}
\item{likelihood}{Automatic or explicit supported likelihood.}
\item{covariates}{Character vector of numeric fixed-effect covariates.}
\item{coords}{Two projected coordinate column names.}
\item{random_effects}{Backward-compatible IID grouping-column shorthand.}
\item{hierarchical_effects}{Optional adaptable hierarchical-effect specification.}
\item{spatial}{Logical; include an SPDE spatial field.}
\item{...}{Additional arguments passed to \code{fit_bottom_up}.}
}
\value{A \code{bottom_up_fit} object.}
\description{High-level interface for Bayesian bottom-up population models.}

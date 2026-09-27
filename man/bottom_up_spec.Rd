\name{bottom_up_spec}
\alias{bottom_up_spec}
\title{Create a bottom-UpR model specification}
\usage{
bottom_up_spec(data, response, buildings = NULL,
  response_type = c("auto", "count", "ppb"), likelihood = "auto",
  covariates = NULL, coords = c("x", "y"), random_effects = NULL,
  hierarchical_effects = NULL, spatial = TRUE)
}
\arguments{
\item{data}{Modelling data frame.}
\item{response}{Response column name.}
\item{buildings}{Optional building/exposure column.}
\item{response_type}{Automatic, count, or PPB response support.}
\item{likelihood}{Automatic or explicit supported likelihood.}
\item{covariates}{Numeric fixed-effect covariate names.}
\item{coords}{Projected coordinate columns.}
\item{random_effects}{Backward-compatible IID grouping-column shorthand.}
\item{hierarchical_effects}{Optional adaptable hierarchical-effect specification.}
\item{spatial}{Logical; include an SPDE field.}
}
\value{A \code{bottom_up_spec} object.}
\description{Validates response support, likelihood compatibility, covariates and optional hierarchical effects.}

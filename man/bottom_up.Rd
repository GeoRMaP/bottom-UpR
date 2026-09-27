\name{bottom_up}
\alias{bottom_up}
\title{Fit a bottom-up population model}
\usage{
bottom_up(data, response, buildings = NULL,
  response_type = c("auto", "count", "ppb"), likelihood = "auto",
  covariates = NULL, coords = c("x", "y"), random_effects = NULL,
  hierarchical_effects = NULL, spatial = TRUE, ...)
}
\description{High-level interface for Gamma/Lognormal PPB and Poisson/Negative-Binomial COUNT models with optional SPDE and adaptable hierarchical effects.}

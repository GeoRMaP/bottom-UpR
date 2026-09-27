\name{bottom_up_inla_selector}
\alias{bottom_up_inla_selector}
\title{Training-only INLA covariate selector}
\usage{
bottom_up_inla_selector(train_data, response, buildings = NULL,
  candidate_covariates, response_type, likelihood = "auto",
  coords = c("x", "y"), hierarchical_effects = NULL,
  spatial = TRUE, level = 0.95, fallback = 0L, ...)
}
\description{Optional selector retaining fixed effects whose posterior credible intervals exclude zero, with an optional strongest-signal fallback.}

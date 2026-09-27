\name{bottom_up_spec}
\alias{bottom_up_spec}
\title{Create a bottom-UpR model specification}
\usage{
bottom_up_spec(data, response, buildings = NULL,
  response_type = c("auto", "count", "ppb"), likelihood = "auto",
  covariates = NULL, coords = c("x", "y"), random_effects = NULL,
  hierarchical_effects = NULL, spatial = TRUE)
}
\description{Creates a validated model specification. Hierarchical effects are optional and may be supplied as IID column-name shorthand or configurable effect lists.}

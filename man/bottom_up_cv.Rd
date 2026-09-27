\name{bottom_up_cv}
\alias{bottom_up_cv}
\title{Cross-validation for bottom-UpR}
\usage{
bottom_up_cv(data, response, buildings = NULL,
  response_type = c("auto", "count", "ppb"), likelihood = "auto",
  covariates = NULL, coords = c("x", "y"),
  hierarchical_effects = NULL, random_effects = NULL,
  method = c("spatial_block", "random", "loro", "loso"),
  folds = 5L, block_size = 100000, group = NULL,
  covariate_selector = NULL, standardize = TRUE,
  selector_args = list(), prediction_draws = 200L,
  mesh = NULL, mesh_args = list(), seed = 123, ...)
}
\description{Performs random, spatial-block, or grouped holdout validation with training-only scaling and optional fold-specific covariate selection.}

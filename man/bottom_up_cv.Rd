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
  mesh = NULL, mesh_args = list(), seed = 123, verbose = TRUE, ...)
}
\arguments{
\item{data}{Modelling data frame.}
\item{response}{Response column.}
\item{buildings}{Optional building/exposure column.}
\item{response_type}{Automatic, count, or PPB response type.}
\item{likelihood}{Automatic or explicit supported likelihood.}
\item{covariates}{Candidate covariates.}
\item{coords}{Projected coordinate columns.}
\item{hierarchical_effects}{Optional adaptable hierarchical effects.}
\item{random_effects}{Backward-compatible IID effect shorthand.}
\item{method}{Random, spatial-block, LORO, or LOSO validation.}
\item{folds}{Number of folds for random/spatial validation.}
\item{block_size}{Spatial block width in coordinate units.}
\item{group}{Holdout grouping column for LORO/LOSO.}
\item{covariate_selector}{Optional training-fold selector function.}
\item{standardize}{Standardize covariates using training-fold parameters.}
\item{selector_args}{Named list passed only to the selector.}
\item{prediction_draws}{Joint posterior draws used for held-out prediction.}
\item{mesh}{Optional common SPDE mesh.}
\item{mesh_args}{Arguments for mesh construction if \code{mesh} is absent.}
\item{seed}{Random seed.}
\item{verbose}{Logical; print fold construction, selection, fit, prediction, and metric progress.}
\item{...}{Additional fitting arguments.}
}
\value{A \code{bottom_up_cv} object.}
\description{Leakage-aware validation with fold-specific scaling and optional fold-specific INLA covariate selection.}

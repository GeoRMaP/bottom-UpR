\name{bottom_up_cv}
\alias{bottom_up_cv}
\title{Cross-validation for bottom-UpR}
\usage{
bottom_up_cv(data, response, buildings = NULL,
  response_type = c("auto", "count", "ppb"), likelihood = "auto",
  covariates = NULL, coords = c("x", "y"), spatial = TRUE,
  hierarchical_effects = NULL, random_effects = NULL,
  method = c("spatial_block", "random", "loro", "loso"),
  folds = 5L, block_size = 100000, group = NULL,
  covariate_selector = NULL, selection_mode = c("separate", "fold"),
  fixed_covariates = NULL, selector_spatial = FALSE, standardize = TRUE,
  selector_args = list(), prediction_draws = 200L,
  cv_profile = c("full", "fast"), fast_keep_selector = FALSE,
  mesh = NULL, mesh_args = list(), seed = 123, verbose = TRUE, ...)
}
\arguments{
\item{data}{Modelling data frame.}
\item{response}{Response column.}
\item{buildings}{Optional building/exposure column.}
\item{response_type}{Automatic, count, or PPB response type.}
\item{likelihood}{Automatic or explicit supported likelihood.}
\item{covariates}{Candidate covariates.}
\item{coords}{Coordinate columns used for spatial modelling and/or spatial-block fold assignment.}
\item{spatial}{Logical; fit an SPDE spatial field within each fold. Set to \code{FALSE} for a faster non-spatial development pass.}
\item{hierarchical_effects}{Optional adaptable hierarchical effects.}
\item{random_effects}{Backward-compatible IID effect shorthand.}
\item{method}{Random, spatial-block, LORO, or LOSO validation.}
\item{folds}{Number of folds for random/spatial validation.}
\item{block_size}{Spatial block width in coordinate units.}
\item{group}{Holdout grouping column for LORO/LOSO.}
\item{covariate_selector}{Optional covariate selector function.}
\item{selection_mode}{\code{"separate"} (default) uses a pre-selected fixed covariate set and performs no selector fits inside CV; \code{"fold"} runs selection independently in each training fold.}
\item{fixed_covariates}{Optional pre-selected covariates used when \code{selection_mode="separate"}. Defaults to \code{covariates}.}
\item{selector_spatial}{When \code{selection_mode="fold"}, fit the selector with an SPDE field. Defaults to \code{FALSE} for speed.}
\item{standardize}{Standardize covariates using training-fold parameters.}
\item{selector_args}{Named list passed only to the selector.}
\item{prediction_draws}{Joint posterior draws used for held-out prediction.}
\item{cv_profile}{\code{"full"} for the requested validation workflow or \code{"fast"} for development/debugging. Fast mode caps random/spatial CV at 3 folds, caps posterior prediction at 50 draws, and skips fold-specific covariate selection unless requested.}
\item{fast_keep_selector}{When \code{cv_profile="fast"}, retain the supplied fold-specific covariate selector.}
\item{mesh}{Optional common SPDE mesh.}
\item{mesh_args}{Arguments for mesh construction if \code{mesh} is absent.}
\item{seed}{Random seed.}
\item{verbose}{Logical; print fold construction, selection, fit, prediction, and metric progress.}
\item{...}{Additional fitting arguments.}
}
\value{A \code{bottom_up_cv} object.}
\description{Cross-validation with fold-specific scaling. Covariate selection is separate by default for speed, with explicit fold-specific nested selection available when required.}

# bottom-UpR

**bottom-UpR** is an R package for Bayesian bottom-up population modelling with
INLA-SPDE. It supports population-count and people-per-building (PPB)
formulations, optional user-defined hierarchical effects, spatial validation,
joint-posterior fine-grid prediction, uncertainty propagation and
administrative aggregation.

## Likelihoods

| Response support | Candidate likelihoods |
|---|---|
| COUNT | Poisson, Negative Binomial |
| PPB | Gamma, Lognormal |

With `likelihood = "auto"`, compatible likelihoods are compared using mean
log CPO on the original response scale.

## Installation

Install R-INLA and then bottom-UpR:

```r
install.packages(
  "INLA",
  repos = c(
    getOption("repos"),
    INLA = "https://inla.r-inla-download.org/R/stable"
  ),
  dep = TRUE
)

install.packages("remotes")
remotes::install_github("GeoRMaP/bottom-UpR")
library(bottom.UpR)
```

## Optional and adaptable hierarchical effects

Hierarchical effects are not required.

```r
fit_simple <- bottom_up(
  ea,
  response = "population",
  buildings = "buildings",
  covariates = c("roads", "water"),
  coords = c("x_m", "y_m"),
  hierarchical_effects = NULL
)
```

For a multilevel model, each effect can be any grouping column. A character
name is shorthand for an IID effect. List syntax also controls how an effect is
handled on a prediction grid.

```r
effects <- list(
  list(column = "source", model = "iid", prediction = "zero"),
  list(column = "settlement", model = "iid", prediction = "known"),
  list(column = "region_settlement", model = "iid", prediction = "known"),
  list(column = "ea_id", model = "iid", prediction = "zero")
)

fit <- bottom_up(
  ea,
  response = "ppb",
  buildings = "buildings",
  response_type = "ppb",
  likelihood = "gamma",
  covariates = covariates,
  coords = c("x_m", "y_m"),
  hierarchical_effects = effects
)
```

This design is not Cameroon-specific: the same mechanism can represent
district, facility, interviewer, survey wave, household, settlement, source or
other hierarchical groupings.

## Cross-validation

`bottom_up_cv()` supports:

- random k-fold CV;
- spatial-block CV;
- grouped leave-one-region-out style validation;
- grouped leave-one-source-out style validation.

Scaling is estimated from the training fold only. A covariate-selector callback
is also executed separately within every training fold.

```r
cv100 <- bottom_up_cv(
  data = ea,
  response = "ppb",
  buildings = "buildings",
  response_type = "ppb",
  likelihood = "gamma",
  covariates = candidate_covariates,
  coords = c("x_m", "y_m"),
  hierarchical_effects = effects,
  method = "spatial_block",
  folds = 5,
  block_size = 100000,
  covariate_selector = bottom_up_inla_selector,
  selector_args = list(fallback = 8),
  prediction_draws = 300
)

cv100$metrics
cv100$metrics_by_fold
cv100$selected_covariates
plot_observed_fitted(cv100)
```

The supplied `bottom_up_inla_selector()` is an optional package rule based on
training-only INLA fixed-effect credible intervals. It is configurable and is
not presented as the unique scientifically correct selection method. A custom
selector can be supplied instead.

## Joint-posterior 100-m prediction

`predict_bottom_up(newdata=...)` takes one set of joint INLA posterior draws
and preserves the same draw identities across all chunks. Fine-grid summaries
and administrative totals therefore propagate posterior dependence coherently.

```r
pred <- predict_bottom_up(
  fit,
  newdata = grid100m,
  population = TRUE,
  draws = 1000,
  chunk_size = 50000,
  return_draws = FALSE,
  aggregate_by = c("region", "division", "national")
)

head(pred$summary)
pred$aggregates$region
pred$aggregates$division
pred$aggregates$national
```

For PPB models, predicted PPB is multiplied by mapped buildings exactly once.
For COUNT models, buildings enter as exposure in the linear predictor and are
not multiplied again. Grid cells with zero mapped buildings are structural
zeros.

## Observation-level results and plots

```r
obs <- observation_predictions(fit, draws = 500)
head(obs)

plot_observed_fitted(fit)
plot_bottom_up_residuals(fit)
plot_posterior_fixed(fit)
plot_posterior_hyperparameters(fit)
plot_spatial_field(fit)

plot_bottom_up_exploration(
  ea,
  response = "population",
  buildings = "buildings",
  covariates = candidate_covariates
)
```

## Cameroon workflow

A complete worked script is included at:

`inst/examples/cameroon_analysis.R`

It demonstrates Gamma-PPB and negative-binomial COUNT fits, optional Cameroon
hierarchical effects, training-fold covariate selection, random and spatial CV,
LORO/LOSO, observation-level fitted results and plots, and chunked 100-m
posterior prediction with national/subnational aggregation.

## Development status

Version 0.2.0 is a development release. The repository has been updated with the
full workflow APIs, but a full `R CMD check` in an R environment with R-INLA
installed is still required before treating the package as a stable release.

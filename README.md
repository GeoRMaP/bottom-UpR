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
  list(column = "source", model = "iid", prediction = "known"),
  list(column = "settlement", model = "iid", prediction = "known"),
  list(column = "region_settlement", model = "iid", prediction = "known"),
  list(column = "ea_id", model = "iid", prediction = "known")
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
other hierarchical groupings. For a national prediction grid, effects can be
omitted without changing the fitted model, for example
`omit_effects = "source",
  marginalize_effects = "ea_id"`.

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
plot_cv_diagnostics(cv100)

# Held-out residual spatial autocorrelation
bottom_up_residual_moran(
  cv100,
  distances = c(50000, 100000, 150000)
)
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
  aggregate_by = c("region", "division", "national"),
  omit_effects = "source",
  marginalize_effects = "ea_id"
)

head(pred$summary)
pred$aggregates$region
pred$aggregates$division
pred$aggregates$national

# Raw draw-wise totals are retained too
pred$aggregate_draws$national
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


## Synthetic-data tutorial

A complete reproducible tutorial is included for users who want to learn the
package without supplying real population data.

Generate tutorial data directly in R:

```r
ea <- simulate_bottom_up_data(
  n = 250,
  seed = 2026,
  include_hierarchy = TRUE
)
```

The full executable workflow is:

`inst/examples/synthetic_tutorial.R`

The accompanying vignette is:

`vignettes/synthetic-tutorial.Rmd`

The tutorial covers PPB and count models, Gamma/Lognormal/Poisson/Negative
Binomial likelihoods, optional hierarchical effects, automatic likelihood
selection, random and spatial-block CV, LORO/LOSO, observation-level posterior
predictions, observed-versus-fitted and posterior diagnostics, residual Moran's
I, fine-grid joint-posterior prediction, structural zeros, and draw-wise
administrative aggregation.


## Adaptive mesh parameters

Mesh geometry can be scaled automatically to the spatial footprint of each data
set rather than using the same absolute distances everywhere.

```r
mesh_par <- bottom_up_mesh_parameters(
  ea,
  coords = c("x_m", "y_m")
)

mesh_par$max_edge
mesh_par$offset
mesh_par$cutoff

mesh <- make_bottom_up_mesh(
  ea,
  coords = c("x_m", "y_m"),
  adaptive = TRUE
)
```

The automatic distances are fractions of the study-area diagonal and remain in
the same units as the projected coordinates. Explicit `max_edge`, `offset`,
or `cutoff` values can still be supplied and override the corresponding
automatic value.

## Optional chunk-wise posterior prediction

Large grids can be processed chunk-by-chunk:

```r
pred_chunked <- predict_bottom_up(
  fit,
  newdata = grid100m,
  population = TRUE,
  draws = 1000,
  chunkwise = TRUE,
  chunk_size = 50000
)
```

For smaller data sets, chunking can be disabled:

```r
pred_all_at_once <- predict_bottom_up(
  fit,
  newdata = grid_small,
  draws = 500,
  chunkwise = FALSE
)
```

Both modes use one synchronized set of joint posterior draws; chunking changes
memory use, not posterior draw identity.


## Citation

If you use **bottom-UpR** in research, please cite:

> Nnanatu, Chibuzor Christopher & Tatem, Andrew J. (2026). *bottom-UpR: Bayesian Bottom-Up Population Modelling with INLA-SPDE*. R package version 0.2.0. https://github.com/GeoRMaP/bottom-UpR

From R, the package citation can be retrieved with:

```r
citation("bottom.UpR")
```


## Full tutorial

A comprehensive package tutorial is included at:

`vignettes/full-tutorial.Rmd`

It covers installation, PPB and COUNT response semantics, synthetic-data
generation, exploratory analysis, optional hierarchical effects, adaptive
SPDE meshes, all supported likelihoods, automatic likelihood selection,
observation-level posterior prediction, diagnostics, training-only covariate
selection, random/spatial/LORO/LOSO validation, residual Moran diagnostics,
fine-grid joint posterior prediction, optional chunk-wise prediction,
hierarchical-effect omission/marginalisation, structural zeros, and draw-wise
administrative aggregation.

After installation, list package vignettes with:

```r
vignette(package = "bottom.UpR")
```

and open the full tutorial with:

```r
vignette("full-tutorial", package = "bottom.UpR")
```


## Shiny gridded-result demos

Two interactive demonstrations are bundled with the package. They use synthetic
gridded posterior results so they can be launched without external data.

Install the optional app dependencies:

```r
install.packages(c("shiny", "leaflet"))
```

Launch the posterior grid explorer:

```r
run_bottom_up_shiny("grid")
```

The grid explorer provides:

- an interactive posterior population map;
- posterior mean, median, lower/upper interval, SD and CV layers;
- mapped-building display;
- structural-zero highlighting;
- region filtering;
- posterior-result distributions;
- region-level summaries.

Launch the model-comparison application:

```r
run_bottom_up_shiny("compare")
```

The comparison app demonstrates how two gridded population surfaces can be
examined interactively using:

- model A and model B population layers;
- absolute grid-cell differences;
- percentage differences;
- uncertainty (CV) comparison;
- observed spatial patterns in disagreements;
- cell-level model scatter plots;
- regional population-total comparison.

The bundled app source files are located under:

`inst/shiny/grid-results/`

and

`inst/shiny/model-comparison/`

They can be used as templates for project-specific Shiny applications that read
the outputs returned by `predict_bottom_up()`.


## Progress reports and warnings

Major modelling functions now report their current stage by default. Progress
messages are prefixed with `[bottom-UpR]`, for example:

```text
[bottom-UpR] Creating model specification
[bottom-UpR] Building SPDE mesh
[bottom-UpR] Running INLA: spatial model
[bottom-UpR] Model fit complete: likelihood=gamma
[bottom-UpR] Sampling joint posterior: draws=500
[bottom-UpR] Predicting chunk: 1/10 (rows 1-50000)
[bottom-UpR] Cross-validation fold: 3/5 [3]
```

Warnings use the same package prefix and flag conditions that should be checked,
including apparent longitude/latitude coordinates in distance-based models,
missing values, sparse hierarchical levels, CPO failures, unusual mesh
settings, large posterior-draw allocations, non-finite prediction inputs,
small validation folds, selector fallbacks, and structural-zero/exposure
issues.

Progress can be silenced without suppressing warnings:

```r
fit <- bottom_up(
  data = ea,
  response = "ppb",
  buildings = "buildings",
  covariates = candidate_covariates,
  coords = c("x_m", "y_m"),
  verbose = FALSE
)
```

The `verbose` argument is available on the major workflow functions, including
model specification/fitting, mesh construction, simulation, covariate
selection, cross-validation, posterior prediction, observation-level
prediction, aggregation, and Shiny demo launching.


## Full-data and chunk-wise synthetic simulation

The tutorial simulator can generate all rows in one block:

```r
ea_full <- simulate_bottom_up_data(
  n = 250,
  seed = 2026,
  chunkwise = FALSE
)
```

or generate the same synthetic data in chunks:

```r
ea_chunked <- simulate_bottom_up_data(
  n = 250,
  seed = 2026,
  chunkwise = TRUE,
  chunk_size = 60
)
```

With the same seed, the tutorial verifies reproducible equivalence:

```r
stopifnot(
  isTRUE(
    all.equal(
      ea_full,
      ea_chunked,
      tolerance = 0,
      check.attributes = TRUE
    )
  )
)
```

The same principle is demonstrated for posterior gridded prediction:
`chunkwise = FALSE` processes the whole grid in one block, while
`chunkwise = TRUE` processes it in memory-bounded chunks while preserving the
same posterior draw identities.

## Warning catalogue and warning demo

List all package-level warning messages/patterns, triggers, and recommended
actions:

```r
bottom_up_warning_catalogue()
```

Run the synthetic warning demonstration:

```r
warning_demo <- system.file(
  "examples",
  "synthetic_warning_demo.R",
  package = "bottom.UpR"
)

source(warning_demo)
```

The demo intentionally triggers deterministic warnings and prints the complete
warning catalogue. Data-dependent INLA diagnostics such as CPO failures or
optimizer-status warnings are listed in the catalogue rather than deliberately
forced through pathological fits.


## Faster cross-validation during development

Full `bottom_up_cv()` can be computationally expensive because each fold may
fit a fold-specific covariate-selection model, fit the final fold model, and
draw from the joint posterior for held-out prediction.

For interactive development and debugging, use:

```r
cv_fast <- bottom_up_cv(
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
  prediction_draws = 200,
  mesh = mesh,
  cv_profile = "fast",
  seed = 123
)
```

Fast mode is deliberately explicit. For random/spatial validation it uses at
most 3 folds, caps held-out posterior prediction at 50 draws, and skips the
fold-specific covariate selector by default. It emits a warning so these
results are not accidentally treated as the final validation analysis.

If you still want fold-specific selection in fast mode:

```r
cv_fast_selected <- bottom_up_cv(
  ...,
  cv_profile = "fast",
  fast_keep_selector = TRUE
)
```

For the final analysis, rerun with:

```r
cv_final <- bottom_up_cv(
  ...,
  cv_profile = "full"
)
```

Additional ways to reduce run time while developing are to use an explicit
likelihood instead of `likelihood="auto"`, temporarily reduce
`prediction_draws`, and reuse a pre-built mesh.

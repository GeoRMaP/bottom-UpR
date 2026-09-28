###############################################################################
# bottom-UpR synthetic-data tutorial
#
# This tutorial is fully reproducible and does not require Cameroon data.
# It demonstrates:
#   1. Synthetic EA data generation: full-data and chunk-wise
#   2. Complete warning catalogue and warning demonstrations
#   3. Exploratory plots
#   3. PPB and COUNT modelling
#   4. Automatic likelihood selection
#   5. Optional hierarchical effects
#   6. Random/spatial/LORO/LOSO validation
#   7. Observation-level posterior predictions
#   8. Posterior and residual plots
#   9. 100-m-like fine-grid joint-posterior prediction
#  10. Draw-wise administrative aggregation
###############################################################################

suppressPackageStartupMessages({
  library(bottom.UpR)
  library(INLA)
})

set.seed(2026)

OUT_DIR <- "synthetic_tutorial_output"
dir.create(OUT_DIR,showWarnings=FALSE,recursive=TRUE)

# =============================================================================
# 1. SYNTHETIC EA DATA
# =============================================================================

# Generate the complete synthetic data in one call.
ea_full <- simulate_bottom_up_data(
  n = 250,
  seed = 2026,
  include_hierarchy = TRUE,
  extent = 300000,
  chunkwise = FALSE,
  verbose = TRUE
)

# Generate exactly the same synthetic data in chunks. This is useful for
# tutorials that want to demonstrate chunk-oriented workflows without changing
# the data-generating process.
ea_chunked <- simulate_bottom_up_data(
  n = 250,
  seed = 2026,
  include_hierarchy = TRUE,
  extent = 300000,
  chunkwise = TRUE,
  chunk_size = 60,
  verbose = TRUE
)

# Same seed + same simulator = same synthetic data in either mode.
simulation_equal <- isTRUE(
  all.equal(
    ea_full,
    ea_chunked,
    tolerance = 0,
    check.attributes = TRUE
  )
)

print(simulation_equal)
stopifnot(simulation_equal)

# Use the full-data version for the remainder of the tutorial.
ea <- ea_full

head(ea)
summary(ea$population_count)
summary(ea$ppb)

# Save both so users can inspect them.
saveRDS(ea_full,file.path(OUT_DIR,"synthetic_ea_full.rds"))
saveRDS(ea_chunked,file.path(OUT_DIR,"synthetic_ea_chunked.rds"))

# =============================================================================
# 1B. PACKAGE WARNING CATALOGUE
# =============================================================================

warning_catalogue <- bottom_up_warning_catalogue()
print(warning_catalogue)

write.csv(
  warning_catalogue,
  file.path(OUT_DIR,"bottom_up_warning_catalogue.csv"),
  row.names=FALSE
)

# The full warning catalogue is printed above. Deterministic examples that
# intentionally trigger warnings are provided in:
#   inst/examples/synthetic_warning_demo.R

candidate_covariates <- c(
  "roads",
  "water",
  "market",
  "slope",
  "night_lights"
)

# =============================================================================
# 2. EXPLORATORY PLOTS
# =============================================================================

png(
  file.path(OUT_DIR,"01_exploratory_count.png"),
  width=1800,height=1500,res=180
)
plot_bottom_up_exploration(
  ea,
  response="population_count",
  buildings="buildings",
  covariates=candidate_covariates
)
dev.off()

png(
  file.path(OUT_DIR,"02_exploratory_ppb.png"),
  width=1800,height=1500,res=180
)
plot_bottom_up_exploration(
  ea,
  response="ppb",
  buildings="buildings",
  covariates=candidate_covariates
)
dev.off()

# =============================================================================
# 3. RESPONSE SUPPORT
# =============================================================================

detect_response_type(ea$population_count)
detect_response_type(ea$ppb)

# =============================================================================
# 4. OPTIONAL HIERARCHICAL EFFECTS
# =============================================================================

effects <- list(
  list(
    column="source",
    model="iid",
    prediction="known",
    label="survey source"
  ),
  list(
    column="settlement",
    model="iid",
    prediction="known",
    label="settlement"
  ),
  list(
    column="region_settlement",
    model="iid",
    prediction="known",
    label="region x settlement"
  ),
  list(
    column="ea_id",
    model="iid",
    prediction="known",
    label="EA residual heterogeneity"
  )
)

# To fit without hierarchy, use:
# effects <- NULL

# =============================================================================
# 5. SPDE MESH
# =============================================================================

mesh_par <- bottom_up_mesh_parameters(
  ea,
  coords=c("x_m","y_m")
)

print(mesh_par[c("max_edge","offset","cutoff","diagonal")])

mesh <- make_bottom_up_mesh(
  ea,
  coords=c("x_m","y_m"),
  adaptive=TRUE
)

saveRDS(mesh,file.path(OUT_DIR,"synthetic_mesh.rds"))

# =============================================================================
# 6. PPB MODELS: GAMMA, LOGNORMAL, AUTOMATIC
# =============================================================================

gamma_fit <- bottom_up(
  data=ea,
  response="ppb",
  buildings="buildings",
  response_type="ppb",
  likelihood="gamma",
  covariates=candidate_covariates,
  coords=c("x_m","y_m"),
  hierarchical_effects=effects,
  spatial=TRUE,
  mesh=mesh,
  prior_range=c(100000,0.05),
  prior_sigma=c(1,0.05),
  config=TRUE
)

lognormal_fit <- bottom_up(
  data=ea,
  response="ppb",
  buildings="buildings",
  response_type="ppb",
  likelihood="lognormal",
  covariates=candidate_covariates,
  coords=c("x_m","y_m"),
  hierarchical_effects=effects,
  spatial=TRUE,
  mesh=mesh,
  prior_range=c(100000,0.05),
  prior_sigma=c(1,0.05),
  config=TRUE
)

ppb_auto <- bottom_up(
  data=ea,
  response="ppb",
  buildings="buildings",
  response_type="ppb",
  likelihood="auto",
  covariates=candidate_covariates,
  coords=c("x_m","y_m"),
  hierarchical_effects=effects,
  spatial=TRUE,
  mesh=mesh,
  prior_range=c(100000,0.05),
  prior_sigma=c(1,0.05),
  config=TRUE
)

print(ppb_auto)
ppb_auto$selection$scores

# =============================================================================
# 7. COUNT MODELS: POISSON, NEGATIVE BINOMIAL, AUTOMATIC
# =============================================================================

poisson_fit <- bottom_up(
  data=ea,
  response="population_count",
  buildings="buildings",
  response_type="count",
  likelihood="poisson",
  covariates=candidate_covariates,
  coords=c("x_m","y_m"),
  hierarchical_effects=effects,
  spatial=TRUE,
  mesh=mesh,
  prior_range=c(100000,0.05),
  prior_sigma=c(1,0.05),
  config=TRUE
)

nb_fit <- bottom_up(
  data=ea,
  response="population_count",
  buildings="buildings",
  response_type="count",
  likelihood="negative binomial",
  covariates=candidate_covariates,
  coords=c("x_m","y_m"),
  hierarchical_effects=effects,
  spatial=TRUE,
  mesh=mesh,
  prior_range=c(100000,0.05),
  prior_sigma=c(1,0.05),
  config=TRUE
)

count_auto <- bottom_up(
  data=ea,
  response="population_count",
  buildings="buildings",
  response_type="count",
  likelihood="auto",
  covariates=candidate_covariates,
  coords=c("x_m","y_m"),
  hierarchical_effects=effects,
  spatial=TRUE,
  mesh=mesh,
  prior_range=c(100000,0.05),
  prior_sigma=c(1,0.05),
  config=TRUE
)

print(count_auto)
count_auto$selection$scores

# =============================================================================
# 8. OBSERVATION-LEVEL POSTERIOR RESULTS
# =============================================================================

gamma_obs <- observation_predictions(
  gamma_fit,
  draws=300,
  seed=11
)

nb_obs <- observation_predictions(
  nb_fit,
  draws=300,
  seed=12
)

write.csv(
  gamma_obs,
  file.path(OUT_DIR,"gamma_observation_predictions.csv"),
  row.names=FALSE
)

write.csv(
  nb_obs,
  file.path(OUT_DIR,"nb_observation_predictions.csv"),
  row.names=FALSE
)

# =============================================================================
# 9. OBSERVED VERSUS FITTED
# =============================================================================

png(
  file.path(OUT_DIR,"03_gamma_observed_fitted.png"),
  width=1400,height=1200,res=180
)
plot_observed_fitted(gamma_fit,draws=300)
dev.off()

png(
  file.path(OUT_DIR,"04_nb_observed_fitted.png"),
  width=1400,height=1200,res=180
)
plot_observed_fitted(nb_fit,draws=300)
dev.off()

# =============================================================================
# 10. RESIDUAL AND POSTERIOR DIAGNOSTIC PLOTS
# =============================================================================

png(
  file.path(OUT_DIR,"05_gamma_residuals.png"),
  width=1800,height=1500,res=180
)
plot_bottom_up_residuals(gamma_fit,draws=300)
dev.off()

png(
  file.path(OUT_DIR,"06_nb_residuals.png"),
  width=1800,height=1500,res=180
)
plot_bottom_up_residuals(nb_fit,draws=300)
dev.off()

png(
  file.path(OUT_DIR,"07_gamma_fixed_effects.png"),
  width=1600,height=1300,res=180
)
plot_posterior_fixed(gamma_fit)
dev.off()

png(
  file.path(OUT_DIR,"08_nb_fixed_effects.png"),
  width=1600,height=1300,res=180
)
plot_posterior_fixed(nb_fit)
dev.off()

png(
  file.path(OUT_DIR,"09_gamma_hyperparameters.png"),
  width=1600,height=1300,res=180
)
plot_posterior_hyperparameters(gamma_fit)
dev.off()

png(
  file.path(OUT_DIR,"10_gamma_source_effect.png"),
  width=1600,height=1300,res=180
)
plot_hierarchical_effects(gamma_fit,effect="source")
dev.off()

png(
  file.path(OUT_DIR,"11_gamma_spatial_field.png"),
  width=1400,height=1400,res=180
)
plot_spatial_field(gamma_fit)
dev.off()

# =============================================================================
# 11. FOLD-SPECIFIC INLA COVARIATE SELECTION
# =============================================================================

# This is a generic tutorial selector supplied by the package.
# In a real paper workflow, replace it with the exact prespecified selection rule
# if a different rule is required.

selector_args <- list(
  fallback=3,
  prior_range=c(100000,0.05),
  prior_sigma=c(1,0.05)
)

# =============================================================================
# 12. RANDOM CV
# =============================================================================

gamma_random <- bottom_up_cv(
  data=ea,
  response="ppb",
  buildings="buildings",
  response_type="ppb",
  likelihood="gamma",
  covariates=candidate_covariates,
  coords=c("x_m","y_m"),
  hierarchical_effects=effects,
  method="random",
  folds=5,
  covariate_selector=bottom_up_inla_selector,
  selector_args=selector_args,
  standardize=TRUE,
  prediction_draws=100,
  mesh=mesh,
  prior_range=c(100000,0.05),
  prior_sigma=c(1,0.05),
  seed=101
)

gamma_random$metrics
gamma_random$metrics_by_fold
gamma_random$selected_covariates

# =============================================================================
# 13. PRIMARY 100-KM SPATIAL-BLOCK CV
# =============================================================================

gamma_spatial100 <- bottom_up_cv(
  data=ea,
  response="ppb",
  buildings="buildings",
  response_type="ppb",
  likelihood="gamma",
  covariates=candidate_covariates,
  coords=c("x_m","y_m"),
  hierarchical_effects=effects,
  method="spatial_block",
  folds=5,
  block_size=100000,
  covariate_selector=bottom_up_inla_selector,
  selector_args=selector_args,
  standardize=TRUE,
  prediction_draws=100,
  mesh=mesh,
  prior_range=c(100000,0.05),
  prior_sigma=c(1,0.05),
  seed=102
)

write.csv(
  gamma_spatial100$predictions,
  file.path(OUT_DIR,"gamma_spatial100_predictions.csv"),
  row.names=FALSE
)

png(
  file.path(OUT_DIR,"12_gamma_spatial100_observed_predicted.png"),
  width=1400,height=1200,res=180
)
plot_observed_fitted(gamma_spatial100)
dev.off()

png(
  file.path(OUT_DIR,"13_gamma_spatial100_diagnostics.png"),
  width=1800,height=1500,res=180
)
plot_cv_diagnostics(gamma_spatial100)
dev.off()

moran100 <- bottom_up_residual_moran(
  gamma_spatial100,
  distances=c(50000,100000,150000)
)

write.csv(
  moran100,
  file.path(OUT_DIR,"gamma_spatial100_moran.csv"),
  row.names=FALSE
)

# =============================================================================
# 14. 50-KM AND 150-KM SPATIAL SENSITIVITY
# =============================================================================

gamma_spatial50 <- bottom_up_cv(
  data=ea,
  response="ppb",
  buildings="buildings",
  response_type="ppb",
  likelihood="gamma",
  covariates=candidate_covariates,
  coords=c("x_m","y_m"),
  hierarchical_effects=effects,
  method="spatial_block",
  folds=5,
  block_size=50000,
  covariate_selector=bottom_up_inla_selector,
  selector_args=selector_args,
  standardize=TRUE,
  prediction_draws=100,
  mesh=mesh,
  seed=103
)

gamma_spatial150 <- bottom_up_cv(
  data=ea,
  response="ppb",
  buildings="buildings",
  response_type="ppb",
  likelihood="gamma",
  covariates=candidate_covariates,
  coords=c("x_m","y_m"),
  hierarchical_effects=effects,
  method="spatial_block",
  folds=5,
  block_size=150000,
  covariate_selector=bottom_up_inla_selector,
  selector_args=selector_args,
  standardize=TRUE,
  prediction_draws=100,
  mesh=mesh,
  seed=104
)

# =============================================================================
# 15. LORO AND LOSO
# =============================================================================

gamma_loro <- bottom_up_cv(
  data=ea,
  response="ppb",
  buildings="buildings",
  response_type="ppb",
  likelihood="gamma",
  covariates=candidate_covariates,
  coords=c("x_m","y_m"),
  hierarchical_effects=effects,
  method="loro",
  group="region",
  covariate_selector=bottom_up_inla_selector,
  selector_args=selector_args,
  standardize=TRUE,
  prediction_draws=100,
  mesh=mesh,
  seed=105
)

gamma_loso <- bottom_up_cv(
  data=ea,
  response="ppb",
  buildings="buildings",
  response_type="ppb",
  likelihood="gamma",
  covariates=candidate_covariates,
  coords=c("x_m","y_m"),
  hierarchical_effects=effects,
  method="loso",
  group="source",
  covariate_selector=bottom_up_inla_selector,
  selector_args=selector_args,
  standardize=TRUE,
  prediction_draws=100,
  mesh=mesh,
  seed=106
)

# =============================================================================
# 16. COUNT-MODEL 100-KM SPATIAL CV
# =============================================================================

nb_spatial100 <- bottom_up_cv(
  data=ea,
  response="population_count",
  buildings="buildings",
  response_type="count",
  likelihood="negative binomial",
  covariates=candidate_covariates,
  coords=c("x_m","y_m"),
  hierarchical_effects=effects,
  method="spatial_block",
  folds=5,
  block_size=100000,
  covariate_selector=bottom_up_inla_selector,
  selector_args=selector_args,
  standardize=TRUE,
  prediction_draws=100,
  mesh=mesh,
  seed=107
)

# =============================================================================
# 17. VALIDATION SUMMARY
# =============================================================================

validation_summary <- rbind(
  cbind(model="Gamma PPB",design="Random",gamma_random$metrics),
  cbind(model="Gamma PPB",design="Spatial 50 km",gamma_spatial50$metrics),
  cbind(model="Gamma PPB",design="Spatial 100 km",gamma_spatial100$metrics),
  cbind(model="Gamma PPB",design="Spatial 150 km",gamma_spatial150$metrics),
  cbind(model="Gamma PPB",design="LORO",gamma_loro$metrics),
  cbind(model="Gamma PPB",design="LOSO",gamma_loso$metrics),
  cbind(model="NB COUNT",design="Spatial 100 km",nb_spatial100$metrics)
)

write.csv(
  validation_summary,
  file.path(OUT_DIR,"validation_summary.csv"),
  row.names=FALSE
)

print(validation_summary)

# =============================================================================
# 18. BUILD A SYNTHETIC FINE GRID
# =============================================================================

# A 100-m national grid would be much larger. Here we use a small regular grid
# so the tutorial runs quickly while exercising the same prediction API.

grid <- expand.grid(
  x_m=seq(5000,295000,length.out=45),
  y_m=seq(5000,295000,length.out=45)
)

ng <- nrow(grid)

set.seed(303)
grid$roads <- stats::rnorm(ng)
grid$water <- stats::rnorm(ng)
grid$market <- stats::rnorm(ng)
grid$slope <- stats::rnorm(ng)
grid$night_lights <- stats::rnorm(ng)

grid$buildings <- stats::rpois(
  ng,
  lambda=exp(2.0 + 0.15*grid$roads - 0.08*grid$slope)
)

grid$settlement <- sample(
  c("urban","periurban","rural"),
  ng,
  replace=TRUE,
  prob=c(.30,.25,.45)
)

grid$region <- sample(
  paste0("R",1:5),
  ng,
  replace=TRUE
)

grid$region_settlement <- interaction(
  grid$region,
  grid$settlement,
  drop=TRUE
)

grid$national <- "Syntheticland"

# Source and EA identifiers are deliberately absent because national-grid
# prediction below omits source and marginalizes EA-level IID heterogeneity.

# =============================================================================
# 19. GAMMA JOINT-POSTERIOR FINE-GRID PREDICTION
# =============================================================================

# All-at-once prediction.
gamma_grid_full <- predict_bottom_up(
  gamma_fit,
  newdata=grid,
  population=TRUE,
  draws=250,
  chunkwise=FALSE,
  return_draws=FALSE,
  aggregate_by=c("national","region"),
  omit_effects="source",
  marginalize_effects="ea_id",
  seed=401,
  verbose=TRUE
)

# Chunk-wise prediction of the same grid with the same posterior seed.
gamma_grid_chunked <- predict_bottom_up(
  gamma_fit,
  newdata=grid,
  population=TRUE,
  draws=250,
  chunkwise=TRUE,
  chunk_size=500,
  return_draws=FALSE,
  aggregate_by=c("national","region"),
  omit_effects="source",
  marginalize_effects="ea_id",
  seed=401,
  verbose=TRUE
)

# Chunking changes memory use, not the posterior target or draw identities.
stopifnot(
  isTRUE(all.equal(
    gamma_grid_full$summary,
    gamma_grid_chunked$summary,
    tolerance=1e-10
  )),
  isTRUE(all.equal(
    gamma_grid_full$aggregates$national,
    gamma_grid_chunked$aggregates$national,
    tolerance=1e-10
  )),
  isTRUE(all.equal(
    gamma_grid_full$aggregates$region,
    gamma_grid_chunked$aggregates$region,
    tolerance=1e-10
  ))
)

gamma_grid <- gamma_grid_chunked

head(gamma_grid$summary)
gamma_grid$aggregates$national
head(gamma_grid$aggregates$region)

# =============================================================================
# 20. NB JOINT-POSTERIOR FINE-GRID PREDICTION
# =============================================================================

nb_grid_full <- predict_bottom_up(
  nb_fit,
  newdata=grid,
  population=TRUE,
  draws=250,
  chunkwise=FALSE,
  return_draws=FALSE,
  aggregate_by=c("national","region"),
  omit_effects="source",
  marginalize_effects="ea_id",
  seed=402,
  verbose=TRUE
)

nb_grid_chunked <- predict_bottom_up(
  nb_fit,
  newdata=grid,
  population=TRUE,
  draws=250,
  chunkwise=TRUE,
  chunk_size=500,
  return_draws=FALSE,
  aggregate_by=c("national","region"),
  omit_effects="source",
  marginalize_effects="ea_id",
  seed=402,
  verbose=TRUE
)

stopifnot(
  isTRUE(all.equal(
    nb_grid_full$summary,
    nb_grid_chunked$summary,
    tolerance=1e-10
  )),
  isTRUE(all.equal(
    nb_grid_full$aggregates$national,
    nb_grid_chunked$aggregates$national,
    tolerance=1e-10
  ))
)

nb_grid <- nb_grid_chunked
nb_grid$aggregates$national

# =============================================================================
# 21. SAVE CELL SUMMARIES
# =============================================================================

gamma_cells <- cbind(
  grid[c(
    "x_m","y_m","buildings",
    "region","settlement"
  )],
  gamma_grid$summary
)

nb_cells <- cbind(
  grid[c(
    "x_m","y_m","buildings",
    "region","settlement"
  )],
  nb_grid$summary
)

saveRDS(
  gamma_cells,
  file.path(OUT_DIR,"gamma_fine_grid_population.rds")
)

saveRDS(
  nb_cells,
  file.path(OUT_DIR,"nb_fine_grid_population.rds")
)

# =============================================================================
# 22. STRUCTURAL ZERO CHECK
# =============================================================================

zero_cells <- which(grid$buildings==0)

if(length(zero_cells)) {
  stopifnot(
    all(
      gamma_grid$summary$mean[zero_cells]==0
    )
  )
  stopifnot(
    all(
      nb_grid$summary$mean[zero_cells]==0
    )
  )
}

# =============================================================================
# 23. COMPARE SYNTHETIC NATIONAL TOTALS
# =============================================================================

national_comparison <- rbind(
  cbind(
    model="Gamma PPB",
    gamma_grid$aggregates$national
  ),
  cbind(
    model="Negative Binomial COUNT",
    nb_grid$aggregates$national
  )
)

print(national_comparison)

write.csv(
  national_comparison,
  file.path(OUT_DIR,"national_population_comparison.csv"),
  row.names=FALSE
)

# =============================================================================
# 24. SAVE FITTED MODELS
# =============================================================================

saveRDS(
  gamma_fit,
  file.path(OUT_DIR,"gamma_fit.rds")
)

saveRDS(
  nb_fit,
  file.path(OUT_DIR,"nb_fit.rds")
)

saveRDS(
  ppb_auto,
  file.path(OUT_DIR,"ppb_auto_fit.rds")
)

saveRDS(
  count_auto,
  file.path(OUT_DIR,"count_auto_fit.rds")
)

message(
  "Synthetic bottom-UpR tutorial finished. Outputs: ",
  normalizePath(OUT_DIR,mustWork=FALSE)
)

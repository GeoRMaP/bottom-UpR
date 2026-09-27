###############################################################################
# Cameroon paper workflow using bottom-UpR v0.2.0
#
# Implements:
#   - Gamma PPB and Negative-Binomial COUNT formulations
#   - optional/adaptable hierarchical effects
#   - full-data INLA covariate selection for final fits
#   - fold-specific training-only INLA selection in CV
#   - random CV, 50/100/150-km spatial CV, LORO and LOSO
#   - observation-level posterior predictions and diagnostics
#   - chunked 100-m joint-posterior prediction
#   - draw-wise national/subnational aggregation
#
# IMPORTANT:
#   1. Coordinates must be projected, not longitude/latitude degrees.
#   2. Edit the column mappings below to match the analytical data.
#   3. The built-in INLA selector is a configurable package rule; replace it
#      with your own callback if the paper uses a different exact selection rule.
###############################################################################

suppressPackageStartupMessages({
  library(bottom.UpR)
  library(INLA)
})

set.seed(123)

# =============================================================================
# 0. FILES AND RUN CONTROLS
# =============================================================================

EA_FILE <- "data/cameroon_ea.csv"
GRID_FILE <- "data/cameroon_grid100m.rds"
OUT_DIR <- "output/cameroon_bottom_upR"

dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)

RUN_FULL_CV <- TRUE
RUN_GRID_PREDICTION <- TRUE

CV_DRAWS <- 300L
GRID_DRAWS <- 1000L
GRID_CHUNK_SIZE <- 50000L

# =============================================================================
# 1. COLUMN MAPPING
# =============================================================================

COL_POP        <- "population"
COL_BUILDINGS  <- "buildings"
COL_X          <- "x_m"
COL_Y          <- "y_m"
COL_SOURCE     <- "source"
COL_REGION     <- "region"
COL_DIVISION   <- "division"
COL_SETTLEMENT <- "settlement"
COL_EA_ID      <- "ea_id"

# Replace with the exact candidate names in the Cameroon analysis dataset.
COVARIATES <- c(
  "cov_conflict",
  "cov_explosions",
  "cov_water",
  "cov_herbaceous",
  "cov_local_roads",
  "cov_marketplaces",
  "cov_slope",
  "cov_night_lights"
)

# =============================================================================
# 2. SPATIAL SPECIFICATION
# =============================================================================

# Distances are metres.
MESH_ARGS <- list(
  max_edge=c(50000,100000),
  offset=c(50000,100000),
  cutoff=30000
)

# P(range < 100 km)=0.05; P(spatial SD > 1)=0.05
PRIOR_RANGE <- c(100000,0.05)
PRIOR_SIGMA <- c(1,0.05)

# =============================================================================
# 3. LOAD AND PREPARE EA DATA
# =============================================================================

ea <- read.csv(EA_FILE, stringsAsFactors=FALSE)

required <- unique(c(
  COL_POP,COL_BUILDINGS,COL_X,COL_Y,
  COL_SOURCE,COL_REGION,COL_SETTLEMENT,COL_EA_ID,
  COVARIATES
))
missing_cols <- setdiff(required,names(ea))
if(length(missing_cols))
  stop("Missing Cameroon modelling columns: ",
       paste(missing_cols,collapse=", "))

keep <- is.finite(ea[[COL_POP]]) &
        is.finite(ea[[COL_BUILDINGS]]) &
        ea[[COL_BUILDINGS]] > 0 &
        is.finite(ea[[COL_X]]) &
        is.finite(ea[[COL_Y]]) &
        complete.cases(ea[,COVARIATES,drop=FALSE])

ea <- ea[keep,,drop=FALSE]
if(!nrow(ea)) stop("No complete Cameroon EAs remain.")

ea$population_count <- as.integer(round(ea[[COL_POP]]))
ea$ppb <- ea$population_count/ea[[COL_BUILDINGS]]
ea$region_settlement <- interaction(
  ea[[COL_REGION]],ea[[COL_SETTLEMENT]],drop=TRUE
)

if(any(!is.finite(ea$ppb) | ea$ppb <= 0))
  stop("Gamma PPB modelling requires strictly positive PPB.")

# =============================================================================
# 4. OPTIONAL / ADAPTABLE HIERARCHICAL EFFECTS
# =============================================================================

# These are manuscript-style Cameroon effects, but all are optional.
#
# prediction="known": use a fitted level when it is represented in newdata.
# prediction="zero" : always omit an effect for every newdata prediction.
#
# Here all effects are retained when their levels are known. National-grid
# prediction later uses omit_effects= to suppress survey-source and EA-specific
# effects while keeping settlement and region-by-settlement effects.

CAMEROON_EFFECTS <- list(
  list(
    column=COL_SOURCE,
    model="iid",
    prediction="known",
    label="survey source"
  ),
  list(
    column=COL_SETTLEMENT,
    model="iid",
    prediction="known",
    label="settlement class"
  ),
  list(
    column="region_settlement",
    model="iid",
    prediction="known",
    label="region x settlement"
  ),
  list(
    column=COL_EA_ID,
    model="iid",
    prediction="known",
    label="EA-level residual heterogeneity"
  )
)

# To fit without hierarchical effects, simply use:
# CAMEROON_EFFECTS <- NULL
#
# To adapt the hierarchy, add/remove arbitrary grouping columns, e.g.
# list(column="district", model="iid", prediction="known")

# =============================================================================
# 5. EXPLORATORY PLOTS
# =============================================================================

grDevices::png(
  file.path(OUT_DIR,"exploratory_population_buildings.png"),
  width=1800,height=1500,res=180
)
plot_bottom_up_exploration(
  ea,
  response="population_count",
  buildings=COL_BUILDINGS,
  covariates=COVARIATES
)
grDevices::dev.off()

# =============================================================================
# 6. FINAL-FIT SCALING
# =============================================================================

# Final national scaling is estimated from all modelling EAs.
final_center <- vapply(ea[COVARIATES],mean,numeric(1),na.rm=TRUE)
final_scale <- vapply(ea[COVARIATES],stats::sd,numeric(1),na.rm=TRUE)
final_scale[!is.finite(final_scale) | final_scale==0] <- 1

ea_scaled <- ea
for(nm in COVARIATES)
  ea_scaled[[nm]] <- (ea_scaled[[nm]]-final_center[[nm]])/final_scale[[nm]]

scaler <- list(center=final_center,scale=final_scale)
saveRDS(scaler,file.path(OUT_DIR,"final_covariate_scaling.rds"))

# =============================================================================
# 7. COMMON CAMEROON SPDE MESH
# =============================================================================

mesh <- do.call(
  make_bottom_up_mesh,
  c(
    list(data=ea_scaled,coords=c(COL_X,COL_Y)),
    MESH_ARGS
  )
)
saveRDS(mesh,file.path(OUT_DIR,"cameroon_mesh.rds"))

# =============================================================================
# 8. FULL-DATA INLA COVARIATE SELECTION FOR FINAL MODELS
# =============================================================================

gamma_covariates <- bottom_up_inla_selector(
  train_data=ea_scaled,
  response="ppb",
  buildings=COL_BUILDINGS,
  candidate_covariates=COVARIATES,
  response_type="ppb",
  likelihood="gamma",
  coords=c(COL_X,COL_Y),
  hierarchical_effects=CAMEROON_EFFECTS,
  spatial=TRUE,
  fallback=min(8,length(COVARIATES)),
  mesh=mesh,
  prior_range=PRIOR_RANGE,
  prior_sigma=PRIOR_SIGMA
)

nb_covariates <- bottom_up_inla_selector(
  train_data=ea_scaled,
  response="population_count",
  buildings=COL_BUILDINGS,
  candidate_covariates=COVARIATES,
  response_type="count",
  likelihood="negative binomial",
  coords=c(COL_X,COL_Y),
  hierarchical_effects=CAMEROON_EFFECTS,
  spatial=TRUE,
  fallback=min(8,length(COVARIATES)),
  mesh=mesh,
  prior_range=PRIOR_RANGE,
  prior_sigma=PRIOR_SIGMA
)

writeLines(gamma_covariates,file.path(OUT_DIR,"gamma_selected_covariates.txt"))
writeLines(nb_covariates,file.path(OUT_DIR,"nb_selected_covariates.txt"))

# =============================================================================
# 9. FINAL GAMMA-PPB AND NB-COUNT MODELS
# =============================================================================

gamma_fit <- bottom_up(
  data=ea_scaled,
  response="ppb",
  buildings=COL_BUILDINGS,
  response_type="ppb",
  likelihood="gamma",
  covariates=gamma_covariates,
  coords=c(COL_X,COL_Y),
  hierarchical_effects=CAMEROON_EFFECTS,
  spatial=TRUE,
  mesh=mesh,
  prior_range=PRIOR_RANGE,
  prior_sigma=PRIOR_SIGMA,
  config=TRUE
)
gamma_fit$scaling <- scaler

nb_fit <- bottom_up(
  data=ea_scaled,
  response="population_count",
  buildings=COL_BUILDINGS,
  response_type="count",
  likelihood="negative binomial",
  covariates=nb_covariates,
  coords=c(COL_X,COL_Y),
  hierarchical_effects=CAMEROON_EFFECTS,
  spatial=TRUE,
  mesh=mesh,
  prior_range=PRIOR_RANGE,
  prior_sigma=PRIOR_SIGMA,
  config=TRUE
)
nb_fit$scaling <- scaler

saveRDS(gamma_fit,file.path(OUT_DIR,"cameroon_gamma_ppb_fit.rds"))
saveRDS(nb_fit,file.path(OUT_DIR,"cameroon_nb_count_fit.rds"))

# =============================================================================
# 10. OBSERVATION-LEVEL POSTERIOR RESULTS
# =============================================================================

gamma_obs <- observation_predictions(gamma_fit,draws=500)
nb_obs <- observation_predictions(nb_fit,draws=500)

gamma_obs[[COL_EA_ID]] <- ea[[COL_EA_ID]]
nb_obs[[COL_EA_ID]] <- ea[[COL_EA_ID]]

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

# Observed versus fitted scatter plots.
grDevices::png(
  file.path(OUT_DIR,"gamma_observed_vs_fitted.png"),
  width=1400,height=1200,res=180
)
plot_observed_fitted(gamma_fit,draws=500)
grDevices::dev.off()

grDevices::png(
  file.path(OUT_DIR,"nb_observed_vs_fitted.png"),
  width=1400,height=1200,res=180
)
plot_observed_fitted(nb_fit,draws=500)
grDevices::dev.off()

# Residual diagnostics.
grDevices::png(
  file.path(OUT_DIR,"gamma_residual_diagnostics.png"),
  width=1800,height=1500,res=180
)
plot_bottom_up_residuals(gamma_fit,draws=500)
grDevices::dev.off()

grDevices::png(
  file.path(OUT_DIR,"nb_residual_diagnostics.png"),
  width=1800,height=1500,res=180
)
plot_bottom_up_residuals(nb_fit,draws=500)
grDevices::dev.off()

# Posterior fixed effects.
grDevices::png(
  file.path(OUT_DIR,"gamma_posterior_fixed_effects.png"),
  width=1600,height=1400,res=180
)
plot_posterior_fixed(gamma_fit)
grDevices::dev.off()

grDevices::png(
  file.path(OUT_DIR,"nb_posterior_fixed_effects.png"),
  width=1600,height=1400,res=180
)
plot_posterior_fixed(nb_fit)
grDevices::dev.off()

# Posterior hyperparameters.
grDevices::png(
  file.path(OUT_DIR,"gamma_posterior_hyperparameters.png"),
  width=1700,height=1400,res=180
)
plot_posterior_hyperparameters(gamma_fit)
grDevices::dev.off()

grDevices::png(
  file.path(OUT_DIR,"nb_posterior_hyperparameters.png"),
  width=1700,height=1400,res=180
)
plot_posterior_hyperparameters(nb_fit)
grDevices::dev.off()

# Spatial posterior means at mesh nodes.
grDevices::png(
  file.path(OUT_DIR,"gamma_spatial_field.png"),
  width=1400,height=1400,res=180
)
plot_spatial_field(gamma_fit)
grDevices::dev.off()

grDevices::png(
  file.path(OUT_DIR,"nb_spatial_field.png"),
  width=1400,height=1400,res=180
)
plot_spatial_field(nb_fit)
grDevices::dev.off()

# Posterior hierarchical effects (first effect shown here; pass effect= to select others).
grDevices::png(
  file.path(OUT_DIR,"gamma_hierarchical_source_effect.png"),
  width=1600,height=1400,res=180
)
plot_hierarchical_effects(gamma_fit,effect=COL_SOURCE)
grDevices::dev.off()

grDevices::png(
  file.path(OUT_DIR,"nb_hierarchical_source_effect.png"),
  width=1600,height=1400,res=180
)
plot_hierarchical_effects(nb_fit,effect=COL_SOURCE)
grDevices::dev.off()

# Save posterior summaries.
write.csv(
  gamma_fit$inla$summary.fixed,
  file.path(OUT_DIR,"gamma_fixed_effects.csv")
)
write.csv(
  nb_fit$inla$summary.fixed,
  file.path(OUT_DIR,"nb_fixed_effects.csv")
)
write.csv(
  gamma_fit$inla$summary.hyperpar,
  file.path(OUT_DIR,"gamma_hyperparameters.csv")
)
write.csv(
  nb_fit$inla$summary.hyperpar,
  file.path(OUT_DIR,"nb_hyperparameters.csv")
)

# =============================================================================
# 11. CROSS-VALIDATION HELPER
# =============================================================================

run_cv <- function(response,response_type,likelihood,method,
                   block_size=NULL,group=NULL,tag) {

  selector_args <- list(
    fallback=min(8,length(COVARIATES)),
    prior_range=PRIOR_RANGE,
    prior_sigma=PRIOR_SIGMA
  )

  args <- list(
    data=ea,
    response=response,
    buildings=COL_BUILDINGS,
    response_type=response_type,
    likelihood=likelihood,
    covariates=COVARIATES,
    coords=c(COL_X,COL_Y),
    hierarchical_effects=CAMEROON_EFFECTS,
    method=method,
    folds=5L,
    group=group,
    covariate_selector=bottom_up_inla_selector,
    selector_args=selector_args,
    standardize=TRUE,
    prediction_draws=CV_DRAWS,
    mesh=mesh,
    seed=123,
    prior_range=PRIOR_RANGE,
    prior_sigma=PRIOR_SIGMA
  )

  if(!is.null(block_size)) args$block_size <- block_size

  z <- do.call(bottom_up_cv,args)

  saveRDS(z,file.path(OUT_DIR,paste0(tag,".rds")))
  write.csv(
    z$predictions,
    file.path(OUT_DIR,paste0(tag,"_predictions.csv")),
    row.names=FALSE
  )
  write.csv(
    z$metrics_by_fold,
    file.path(OUT_DIR,paste0(tag,"_metrics_by_fold.csv")),
    row.names=FALSE
  )

  grDevices::png(
    file.path(OUT_DIR,paste0(tag,"_observed_vs_predicted.png")),
    width=1400,height=1200,res=180
  )
  plot_observed_fitted(z)
  grDevices::dev.off()

  grDevices::png(
    file.path(OUT_DIR,paste0(tag,"_diagnostics.png")),
    width=1800,height=1500,res=180
  )
  plot_cv_diagnostics(z)
  grDevices::dev.off()

  moran <- bottom_up_residual_moran(
    z,distances=c(50000,100000,150000)
  )
  write.csv(
    moran,
    file.path(OUT_DIR,paste0(tag,"_residual_moran.csv")),
    row.names=FALSE
  )

  grDevices::png(
    file.path(OUT_DIR,paste0(tag,"_residual_moran.png")),
    width=1400,height=1100,res=180
  )
  plot_residual_moran(z,distances=c(50000,100000,150000))
  grDevices::dev.off()

  z
}

# =============================================================================
# 12. RANDOM / SPATIAL / LORO / LOSO VALIDATION
# =============================================================================

if(RUN_FULL_CV) {

  # Primary 100-km spatial-block CV.
  gamma_cv100 <- run_cv(
    "ppb","ppb","gamma","spatial_block",
    block_size=100000,tag="gamma_cv_spatial_100km"
  )
  nb_cv100 <- run_cv(
    "population_count","count","negative binomial","spatial_block",
    block_size=100000,tag="nb_cv_spatial_100km"
  )

  # Spatial sensitivity.
  gamma_cv50 <- run_cv(
    "ppb","ppb","gamma","spatial_block",
    block_size=50000,tag="gamma_cv_spatial_50km"
  )
  nb_cv50 <- run_cv(
    "population_count","count","negative binomial","spatial_block",
    block_size=50000,tag="nb_cv_spatial_50km"
  )

  gamma_cv150 <- run_cv(
    "ppb","ppb","gamma","spatial_block",
    block_size=150000,tag="gamma_cv_spatial_150km"
  )
  nb_cv150 <- run_cv(
    "population_count","count","negative binomial","spatial_block",
    block_size=150000,tag="nb_cv_spatial_150km"
  )

  # Random five-fold comparator.
  gamma_cv_random <- run_cv(
    "ppb","ppb","gamma","random",
    tag="gamma_cv_random"
  )
  nb_cv_random <- run_cv(
    "population_count","count","negative binomial","random",
    tag="nb_cv_random"
  )

  # Leave-one-region-out.
  gamma_loro <- run_cv(
    "ppb","ppb","gamma","loro",
    group=COL_REGION,tag="gamma_loro"
  )
  nb_loro <- run_cv(
    "population_count","count","negative binomial","loro",
    group=COL_REGION,tag="nb_loro"
  )

  # Leave-one-source-out.
  gamma_loso <- run_cv(
    "ppb","ppb","gamma","loso",
    group=COL_SOURCE,tag="gamma_loso"
  )
  nb_loso <- run_cv(
    "population_count","count","negative binomial","loso",
    group=COL_SOURCE,tag="nb_loso"
  )

  cv_summary <- rbind(
    cbind(model="Gamma",design="Spatial 50 km",gamma_cv50$metrics),
    cbind(model="NB",design="Spatial 50 km",nb_cv50$metrics),
    cbind(model="Gamma",design="Spatial 100 km",gamma_cv100$metrics),
    cbind(model="NB",design="Spatial 100 km",nb_cv100$metrics),
    cbind(model="Gamma",design="Spatial 150 km",gamma_cv150$metrics),
    cbind(model="NB",design="Spatial 150 km",nb_cv150$metrics),
    cbind(model="Gamma",design="Random 5-fold",gamma_cv_random$metrics),
    cbind(model="NB",design="Random 5-fold",nb_cv_random$metrics),
    cbind(model="Gamma",design="LORO",gamma_loro$metrics),
    cbind(model="NB",design="LORO",nb_loro$metrics),
    cbind(model="Gamma",design="LOSO",gamma_loso$metrics),
    cbind(model="NB",design="LOSO",nb_loso$metrics)
  )

  write.csv(
    cv_summary,
    file.path(OUT_DIR,"cameroon_cv_summary.csv"),
    row.names=FALSE
  )
}

# =============================================================================
# 13. 100-M JOINT-POSTERIOR PREDICTION AND ADMINISTRATIVE AGGREGATION
# =============================================================================

if(RUN_GRID_PREDICTION && file.exists(GRID_FILE)) {

  grid100 <- if(grepl("\\.rds$",GRID_FILE,ignore.case=TRUE))
    readRDS(GRID_FILE) else read.csv(GRID_FILE,stringsAsFactors=FALSE)

  grid_required <- unique(c(
    COL_BUILDINGS,COL_X,COL_Y,COL_REGION,COL_SETTLEMENT,
    COVARIATES
  ))
  miss <- setdiff(grid_required,names(grid100))
  if(length(miss))
    stop("100-m grid is missing: ",paste(miss,collapse=", "))

  grid100$region_settlement <- interaction(
    grid100[[COL_REGION]],grid100[[COL_SETTLEMENT]],drop=TRUE
  )
  grid100$national <- "Cameroon"

  aggregate_cols <- c("national",COL_REGION)
  if(COL_DIVISION %in% names(grid100))
    aggregate_cols <- c(aggregate_cols,COL_DIVISION)

  gamma_grid <- predict_bottom_up(
    gamma_fit,
    newdata=grid100,
    population=TRUE,
    draws=GRID_DRAWS,
    chunk_size=GRID_CHUNK_SIZE,
    return_draws=FALSE,
    aggregate_by=aggregate_cols,
    omit_effects=c(COL_SOURCE,COL_EA_ID),
    seed=1001
  )

  nb_grid <- predict_bottom_up(
    nb_fit,
    newdata=grid100,
    population=TRUE,
    draws=GRID_DRAWS,
    chunk_size=GRID_CHUNK_SIZE,
    return_draws=FALSE,
    aggregate_by=aggregate_cols,
    omit_effects=c(COL_SOURCE,COL_EA_ID),
    seed=1001
  )

  # 100-m cell summaries.
  gamma_cell <- cbind(
    grid100[,intersect(
      c(COL_X,COL_Y,COL_REGION,COL_DIVISION,COL_SETTLEMENT,COL_BUILDINGS),
      names(grid100)
    ),drop=FALSE],
    gamma_grid$summary
  )
  nb_cell <- cbind(
    grid100[,intersect(
      c(COL_X,COL_Y,COL_REGION,COL_DIVISION,COL_SETTLEMENT,COL_BUILDINGS),
      names(grid100)
    ),drop=FALSE],
    nb_grid$summary
  )

  saveRDS(gamma_cell,file.path(OUT_DIR,"gamma_100m_population.rds"))
  saveRDS(nb_cell,file.path(OUT_DIR,"nb_100m_population.rds"))

  # Posterior totals are calculated from draw-wise sums, not by summing
  # marginal interval endpoints.
  for(g in names(gamma_grid$aggregates)) {
    write.csv(
      gamma_grid$aggregates[[g]],
      file.path(OUT_DIR,paste0("gamma_population_",g,".csv")),
      row.names=FALSE
    )
  }

  for(g in names(nb_grid$aggregates)) {
    write.csv(
      nb_grid$aggregates[[g]],
      file.path(OUT_DIR,paste0("nb_population_",g,".csv")),
      row.names=FALSE
    )
  }

} else if(RUN_GRID_PREDICTION) {
  message(
    "GRID_FILE was not found. Final EA models/CV were produced, ",
    "but 100-m prediction was skipped: ",GRID_FILE
  )
}

message(
  "\nCameroon bottom-UpR workflow completed. Outputs written to: ",
  normalizePath(OUT_DIR,mustWork=FALSE)
)

###############################################################################
# Cameroon application using bottom-UpR
#
# Fits the two paper formulations:
#   1) Gamma model for people per building (PPB)
#   2) Negative-binomial model for population COUNT with log(buildings) exposure
#
# The script assumes that EA centroid coordinates have already been projected
# into metres. Do NOT use longitude/latitude degrees for the SPDE mesh.
###############################################################################

suppressPackageStartupMessages({
  library(bottom.UpR)
  library(INLA)
})

set.seed(123)

# ---------------------------------------------------------------------------
# 1. USER SETTINGS
# ---------------------------------------------------------------------------

EA_FILE <- "data/cameroon_ea.csv"
OUT_DIR <- "output/cameroon"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Map these names to your analytical dataset.
COL_POP        <- "population"
COL_BUILDINGS  <- "buildings"
COL_X          <- "x_m"             # projected easting, metres
COL_Y          <- "y_m"             # projected northing, metres
COL_SOURCE     <- "source"
COL_REGION     <- "region"
COL_SETTLEMENT <- "settlement"
COL_EA_ID      <- "ea_id"

# Supply the final manuscript covariates here.
# IMPORTANT: for cross-validation, selection must be repeated inside each
# training fold; do not preselect on the full dataset and reuse the set in CV.
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

# Spatial specification. Units are metres.
MESH_MAX_EDGE <- c(50000, 100000)
MESH_OFFSET   <- c(50000, 100000)
MESH_CUTOFF   <- 30000

# PC priors:
# P(range < 100 km) = 0.05
# P(spatial SD > 1) = 0.05
PRIOR_RANGE <- c(100000, 0.05)
PRIOR_SIGMA <- c(1, 0.05)

# ---------------------------------------------------------------------------
# 2. LOAD AND CHECK EA DATA
# ---------------------------------------------------------------------------

ea <- read.csv(EA_FILE, stringsAsFactors = FALSE)

required <- unique(c(
  COL_POP, COL_BUILDINGS, COL_X, COL_Y,
  COL_SOURCE, COL_REGION, COL_SETTLEMENT, COL_EA_ID,
  COVARIATES
))
missing_cols <- setdiff(required, names(ea))
if (length(missing_cols)) {
  stop("Missing required Cameroon columns: ",
       paste(missing_cols, collapse = ", "))
}

# Positive building exposure is required in model-training EAs.
keep <- is.finite(ea[[COL_POP]]) &
        is.finite(ea[[COL_BUILDINGS]]) &
        ea[[COL_BUILDINGS]] > 0 &
        is.finite(ea[[COL_X]]) &
        is.finite(ea[[COL_Y]])

keep <- keep & complete.cases(ea[, COVARIATES, drop = FALSE])
ea <- ea[keep, , drop = FALSE]

if (!nrow(ea)) stop("No complete Cameroon EAs remain after data checks.")

# ---------------------------------------------------------------------------
# 3. DERIVE BOTH RESPONSE FORMULATIONS
# ---------------------------------------------------------------------------

ea$population_count <- as.integer(round(ea[[COL_POP]]))
ea$ppb <- ea$population_count / ea[[COL_BUILDINGS]]

if (any(ea$ppb <= 0 | !is.finite(ea$ppb))) {
  stop("Gamma PPB training response must be strictly positive.")
}

# ---------------------------------------------------------------------------
# 4. HIERARCHICAL EFFECT INDICES
# ---------------------------------------------------------------------------

# bottom-UpR currently expects integer-safe grouping indices for robust INLA use.
ea$source_id <- as.integer(factor(ea[[COL_SOURCE]]))
ea$settlement_id <- as.integer(factor(ea[[COL_SETTLEMENT]]))
ea$region_settlement_id <- as.integer(
  interaction(ea[[COL_REGION]], ea[[COL_SETTLEMENT]], drop = TRUE)
)
ea$ea_iid <- seq_len(nrow(ea))

random_effects <- c(
  "source_id",
  "settlement_id",
  "region_settlement_id",
  "ea_iid"
)

# ---------------------------------------------------------------------------
# 5. STANDARDISE CONTINUOUS COVARIATES
# ---------------------------------------------------------------------------

# For the final national fit, scaling is estimated from all modelling EAs.
# During CV this must instead be estimated from training EAs only.
scaling <- data.frame(
  covariate = COVARIATES,
  mean = vapply(ea[COVARIATES], mean, numeric(1), na.rm = TRUE),
  sd = vapply(ea[COVARIATES], sd, numeric(1), na.rm = TRUE),
  row.names = NULL
)
scaling$sd[!is.finite(scaling$sd) | scaling$sd == 0] <- 1

for (j in seq_len(nrow(scaling))) {
  nm <- scaling$covariate[j]
  ea[[nm]] <- (ea[[nm]] - scaling$mean[j]) / scaling$sd[j]
}

write.csv(scaling, file.path(OUT_DIR, "covariate_scaling.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# 6. BUILD ONE COMMON CAMEROON SPDE MESH
# ---------------------------------------------------------------------------

mesh <- make_bottom_up_mesh(
  data = ea,
  coords = c(COL_X, COL_Y),
  max_edge = MESH_MAX_EDGE,
  offset = MESH_OFFSET,
  cutoff = MESH_CUTOFF
)

saveRDS(mesh, file.path(OUT_DIR, "cameroon_mesh.rds"))

# ---------------------------------------------------------------------------
# 7. GAMMA PPB MODEL
# ---------------------------------------------------------------------------

gamma_fit <- bottom_up(
  data = ea,
  response = "ppb",
  buildings = COL_BUILDINGS,
  response_type = "ppb",
  likelihood = "gamma",
  covariates = COVARIATES,
  coords = c(COL_X, COL_Y),
  random_effects = random_effects,
  spatial = TRUE,
  mesh = mesh,
  prior_range = PRIOR_RANGE,
  prior_sigma = PRIOR_SIGMA,
  config = TRUE
)

saveRDS(gamma_fit, file.path(OUT_DIR, "cameroon_gamma_ppb_fit.rds"))

# ---------------------------------------------------------------------------
# 8. NEGATIVE-BINOMIAL COUNT MODEL
# ---------------------------------------------------------------------------

nb_fit <- bottom_up(
  data = ea,
  response = "population_count",
  buildings = COL_BUILDINGS,
  response_type = "count",
  likelihood = "negative binomial",
  covariates = COVARIATES,
  coords = c(COL_X, COL_Y),
  random_effects = random_effects,
  spatial = TRUE,
  mesh = mesh,
  prior_range = PRIOR_RANGE,
  prior_sigma = PRIOR_SIGMA,
  config = TRUE
)

saveRDS(nb_fit, file.path(OUT_DIR, "cameroon_nb_count_fit.rds"))

# ---------------------------------------------------------------------------
# 9. OPTIONAL WITHIN-SUPPORT AUTOMATIC LIKELIHOOD COMPARISONS
# ---------------------------------------------------------------------------

# PPB: Gamma versus Lognormal
ppb_auto <- bottom_up(
  data = ea,
  response = "ppb",
  buildings = COL_BUILDINGS,
  response_type = "ppb",
  likelihood = "auto",
  covariates = COVARIATES,
  coords = c(COL_X, COL_Y),
  random_effects = random_effects,
  spatial = TRUE,
  mesh = mesh,
  prior_range = PRIOR_RANGE,
  prior_sigma = PRIOR_SIGMA,
  config = TRUE
)

# COUNT: Poisson versus Negative Binomial
count_auto <- bottom_up(
  data = ea,
  response = "population_count",
  buildings = COL_BUILDINGS,
  response_type = "count",
  likelihood = "auto",
  covariates = COVARIATES,
  coords = c(COL_X, COL_Y),
  random_effects = random_effects,
  spatial = TRUE,
  mesh = mesh,
  prior_range = PRIOR_RANGE,
  prior_sigma = PRIOR_SIGMA,
  config = TRUE
)

saveRDS(ppb_auto, file.path(OUT_DIR, "cameroon_ppb_auto_fit.rds"))
saveRDS(count_auto, file.path(OUT_DIR, "cameroon_count_auto_fit.rds"))

# ---------------------------------------------------------------------------
# 10. MODEL SUMMARIES
# ---------------------------------------------------------------------------

model_table <- data.frame(
  formulation = c(
    "Gamma PPB",
    "Negative Binomial COUNT",
    "Automatic PPB",
    "Automatic COUNT"
  ),
  selected_likelihood = c(
    gamma_fit$likelihood,
    nb_fit$likelihood,
    ppb_auto$likelihood,
    count_auto$likelihood
  ),
  WAIC = c(
    gamma_fit$inla$waic$waic,
    nb_fit$inla$waic$waic,
    ppb_auto$inla$waic$waic,
    count_auto$inla$waic$waic
  )
)

write.csv(model_table,
          file.path(OUT_DIR, "cameroon_model_summary.csv"),
          row.names = FALSE)

write.csv(
  gamma_fit$inla$summary.fixed,
  file.path(OUT_DIR, "gamma_fixed_effects.csv")
)
write.csv(
  gamma_fit$inla$summary.hyperpar,
  file.path(OUT_DIR, "gamma_hyperparameters.csv")
)
write.csv(
  nb_fit$inla$summary.fixed,
  file.path(OUT_DIR, "nb_fixed_effects.csv")
)
write.csv(
  nb_fit$inla$summary.hyperpar,
  file.path(OUT_DIR, "nb_hyperparameters.csv")
)

if (!is.null(ppb_auto$selection)) {
  write.csv(
    ppb_auto$selection$scores,
    file.path(OUT_DIR, "ppb_likelihood_selection.csv"),
    row.names = FALSE
  )
}
if (!is.null(count_auto$selection)) {
  write.csv(
    count_auto$selection$scores,
    file.path(OUT_DIR, "count_likelihood_selection.csv"),
    row.names = FALSE
  )
}

# ---------------------------------------------------------------------------
# 11. FITTED-SUPPORT CHECKS
# ---------------------------------------------------------------------------

gamma_ea <- predict_bottom_up(gamma_fit, population = TRUE)
nb_ea <- predict_bottom_up(nb_fit, population = TRUE)

gamma_ea$observed <- ea$population_count
nb_ea$observed <- ea$population_count

metrics <- function(obs, pred) {
  data.frame(
    RMSE = sqrt(mean((pred - obs)^2, na.rm = TRUE)),
    MAE = mean(abs(pred - obs), na.rm = TRUE),
    total_ratio = sum(pred, na.rm = TRUE) / sum(obs, na.rm = TRUE)
  )
}

fit_metrics <- rbind(
  cbind(model = "Gamma PPB", metrics(gamma_ea$observed, gamma_ea$mean)),
  cbind(model = "Negative Binomial COUNT",
        metrics(nb_ea$observed, nb_ea$mean))
)

write.csv(gamma_ea, file.path(OUT_DIR, "gamma_fitted_EA.csv"), row.names = FALSE)
write.csv(nb_ea, file.path(OUT_DIR, "nb_fitted_EA.csv"), row.names = FALSE)
write.csv(fit_metrics, file.path(OUT_DIR, "fitted_support_metrics.csv"),
          row.names = FALSE)

# ---------------------------------------------------------------------------
# 12. IMPORTANT PAPER-WORKFLOW NOTE
# ---------------------------------------------------------------------------

message(
  paste0(
    "\nCameroon final models fitted.\n",
    "For the paper's validation, covariate selection and scaling must be ",
    "repeated within each training fold for random CV, 50/100/150-km spatial ",
    "block CV, LORO and LOSO.\n",
    "The current bottom-UpR v0.1 fitted-support predictor is NOT yet the full ",
    "100-m newdata posterior prediction engine. Do not use fitted EA predictions ",
    "as a substitute for the paper's 100-m joint-posterior workflow.\n"
  )
)

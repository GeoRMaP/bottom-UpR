###############################################################################
# bottom-UpR synthetic warning demonstration
#
# This script intentionally triggers deterministic package warnings so users can
# learn what they mean. It also prints the complete warning catalogue, including
# warnings that depend on INLA diagnostics and therefore cannot be forced
# reliably in a deterministic tutorial.
###############################################################################

suppressPackageStartupMessages({
  library(bottom.UpR)
})

options(warn=1)

cat("\n=== COMPLETE WARNING CATALOGUE ===\n")
warning_catalogue <- bottom_up_warning_catalogue()
print(warning_catalogue,row.names=FALSE)

cat("\n=== 1. SMALL SYNTHETIC SAMPLE ===\n")
small <- simulate_bottom_up_data(
  n=30,
  seed=1,
  verbose=TRUE
)

cat("\n=== 2. CHUNK SIZE TOO LARGE FOR CHUNKED SIMULATION ===\n")
small_chunk <- simulate_bottom_up_data(
  n=30,
  seed=1,
  chunkwise=TRUE,
  chunk_size=100,
  verbose=TRUE
)

cat("\n=== 3. RESPONSE-TYPE OVERRIDE ===\n")
override_data <- data.frame(
  y=c(1,2,3,4,5),
  x_m=seq(0,4000,length.out=5),
  y_m=seq(0,4000,length.out=5)
)

invisible(
  bottom_up_spec(
    override_data,
    response="y",
    response_type="ppb",
    coords=c("x_m","y_m"),
    spatial=FALSE,
    verbose=TRUE
  )
)

cat("\n=== 4. COUNT MODEL WITHOUT BUILDING EXPOSURE ===\n")
invisible(
  bottom_up_spec(
    override_data,
    response="y",
    response_type="count",
    spatial=FALSE,
    verbose=TRUE
  )
)

cat("\n=== 5. PPB MODEL WITHOUT BUILDINGS FOR POPULATION RECONSTRUCTION ===\n")
ppb_data <- override_data
ppb_data$ppb <- c(1.1,1.4,2.2,1.8,2.7)

invisible(
  bottom_up_spec(
    ppb_data,
    response="ppb",
    response_type="ppb",
    spatial=FALSE,
    verbose=TRUE
  )
)

cat("\n=== 6. MISSING MODEL INPUT ===\n")
missing_data <- ppb_data
missing_data$cov1 <- c(0.1,NA,0.3,0.4,0.5)

invisible(
  bottom_up_spec(
    missing_data,
    response="ppb",
    response_type="ppb",
    covariates="cov1",
    spatial=FALSE,
    verbose=TRUE
  )
)

cat("\n=== 7. SPARSE HIERARCHICAL LEVEL ===\n")
sparse_data <- ppb_data
sparse_data$source <- c("A","A","A","A","rare")

invisible(
  bottom_up_spec(
    sparse_data,
    response="ppb",
    response_type="ppb",
    hierarchical_effects=list("source"),
    spatial=FALSE,
    verbose=TRUE
  )
)

cat("\n=== 8. COORDINATES THAT LOOK LIKE LONGITUDE/LATITUDE ===\n")
geo <- data.frame(
  ppb=seq(1.1,3,length.out=25),
  lon=seq(8,15,length.out=25),
  lat=seq(2,12,length.out=25)
)

invisible(
  bottom_up_spec(
    geo,
    response="ppb",
    response_type="ppb",
    coords=c("lon","lat"),
    spatial=TRUE,
    verbose=TRUE
  )
)

cat("\n=== 9. UNUSUAL ADAPTIVE MESH EDGE ORDER ===\n")
invisible(
  bottom_up_mesh_parameters(
    data.frame(
      x_m=c(0,1000,2000,3000),
      y_m=c(0,1500,500,2500)
    ),
    coords=c("x_m","y_m"),
    inner_edge_fraction=0.20,
    outer_edge_fraction=0.10,
    verbose=TRUE
  )
)

cat("\n=== 10. MISSING AGGREGATION GROUP ===\n")
draws <- matrix(
  c(10,11,12,20,21,22,30,31,32),
  nrow=3,
  byrow=TRUE
)

print(
  aggregate_bottom_up(
    draws,
    group=c("A",NA,"B"),
    verbose=TRUE
  )
)

cat("\n=== 11. NON-FINITE POSTERIOR DRAW IN AGGREGATION ===\n")
draws_bad <- draws
draws_bad[2,2] <- NA_real_

print(
  aggregate_bottom_up(
    draws_bad,
    group=c("A","A","B"),
    verbose=TRUE
  )
)

cat("\n=== 12. MORAN DISTANCE BAND WITH NO LINKS ===\n")
fake_cv <- structure(
  list(
    predictions=data.frame(
      coord_x=c(0,100000,200000,300000),
      coord_y=c(0,100000,200000,300000),
      observed=c(10,12,11,13),
      predicted=c(9,11,12,13)
    )
  ),
  class="bottom_up_cv"
)

print(
  bottom_up_residual_moran(
    fake_cv,
    distances=c(1,50000),
    verbose=TRUE
  )
)

if(requireNamespace("INLA",quietly=TRUE)) {
  cat("\n=== 13. PREDICTION-SPECIFIC WARNINGS ===\n")

  fit_data <- simulate_bottom_up_data(
    n=60,
    seed=55,
    verbose=FALSE
  )

  fit <- bottom_up(
    data=fit_data,
    response="population_count",
    buildings="buildings",
    response_type="count",
    likelihood="poisson",
    covariates=c("roads","water"),
    coords=c("x_m","y_m"),
    spatial=FALSE,
    config=TRUE,
    verbose=FALSE
  )

  # One-chunk warning.
  invisible(
    predict_bottom_up(
      fit,
      newdata=fit_data[1:10,],
      draws=10,
      chunkwise=TRUE,
      chunk_size=100,
      return_draws=FALSE,
      verbose=TRUE
    )
  )

  # Invalid exposure warnings. These are treated as structural-zero safeguards
  # for prediction, but the underlying exposure data should normally be fixed.
  bad_grid <- fit_data[1:10,]
  bad_grid$buildings[1] <- NA_real_
  bad_grid$buildings[2] <- -1

  invisible(
    predict_bottom_up(
      fit,
      newdata=bad_grid,
      draws=10,
      chunkwise=TRUE,
      chunk_size=5,
      return_draws=FALSE,
      verbose=TRUE
    )
  )
} else {
  cat(
    "\nINLA is not installed, so prediction-specific warning examples were skipped.\n"
  )
}

cat("\n=== WARNINGS THAT ARE DIAGNOSTIC/DATA-DEPENDENT ===\n")
cat(
  paste(
    "- INLA CPO failures\n",
    "- non-zero INLA optimization mode status\n",
    "- non-finite automatic likelihood-selection scores\n",
    "- selector fallback when no credible interval excludes zero\n",
    "- zero covariates selected inside a validation fold\n",
    "- very small cross-validation folds\n",
    "- very large return_draws allocations\n",
    sep=""
  )
)
cat(
  "These are included in bottom_up_warning_catalogue() but are not forced here ",
  "because doing so would require deliberately pathological fits or very large ",
  "memory allocations.\n"
)

cat("\nWarning demonstration complete.\n")

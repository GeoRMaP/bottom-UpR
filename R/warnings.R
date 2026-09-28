#' List bottom-UpR warning messages and recommended actions
#'
#' Returns a catalogue of package-level warning messages that can be emitted by
#' the main modelling workflow. The catalogue is intended for tutorials,
#' reproducibility reports, and troubleshooting. Exact INLA warnings originating
#' inside the INLA package itself are not reproduced here.
#'
#' @return A data frame with warning category, message, trigger, and suggested
#'   action.
#' @export
bottom_up_warning_catalogue <- function() {
  data.frame(
    category=c(
      "coordinates","missing_data","hierarchy","response","response",
      "mesh","mesh","fit","fit","likelihood","selection","selection",
      "cross_validation","cross_validation","prediction","prediction",
      "prediction","prediction","aggregation","aggregation","moran",
      "simulation"
    ),
    message=c(
      "Coordinates appear to be longitude/latitude degrees. SPDE mesh distances, block sizes, and range priors use coordinate units; project coordinates before distance-based spatial modelling.",
      "Missing values detected in model or prediction inputs. Preprocess explicitly before fitting or prediction.",
      "A hierarchical effect contains levels with fewer than 3 observations; posterior estimates for sparse levels may be weakly identified.",
      "An explicit response_type overrides automatic response detection.",
      "COUNT models without buildings omit the log(buildings) exposure offset; PPB models without buildings cannot reconstruct population automatically.",
      "Inner max_edge is not smaller than outer max_edge; check mesh settings.",
      "cutoff exceeds inner max_edge; this may over-thin input locations.",
      "INLA reported CPO failures for some observations; inspect fit$cpo$failure.",
      "INLA optimization returned a non-zero mode status; inspect fit$mode.",
      "Likelihood comparison produced non-finite mean log CPO scores.",
      "No candidate credible intervals excluded zero; fallback covariates are retained when fallback > 0.",
      "No candidate covariates were supplied, or none were found in INLA fixed-effect summaries.",
      "Some validation folds contain fewer than 3 held-out observations; metrics may be unstable.",
      "A validation fold selected zero covariates; an intercept/spatial/hierarchical model is fitted.",
      "chunkwise=TRUE but chunk_size is at least nrow(newdata); prediction uses one chunk.",
      "return_draws=TRUE may allocate a very large cell-by-draw matrix.",
      "Non-finite building exposure in prediction data is treated as a structural zero.",
      "Negative building exposure in prediction data is treated as a structural zero.",
      "Missing aggregation group labels are excluded from posterior aggregation.",
      "Non-finite posterior draw values can produce non-finite aggregate summaries.",
      "A Moran distance band contains no neighbour links.",
      "Small synthetic sample size may make spatial and cross-validation examples unstable; chunkwise simulation with chunk_size >= n uses one chunk."
    ),
    trigger=c(
      "Coordinate names/values look geographic while distance-based modelling is requested.",
      "Required model or prediction columns contain NA values.",
      "Any hierarchical grouping level has fewer than 3 records.",
      "response_type is explicitly set to a value different from automatic detection.",
      "Buildings are absent for COUNT exposure or PPB population reconstruction.",
      "max_edge[1] >= max_edge[2].",
      "cutoff > max_edge[1].",
      "INLA returns positive CPO failure indicators.",
      "INLA returns a non-zero optimization mode status.",
      "All compared likelihood mean-log-CPO scores are non-finite.",
      "Selector finds no interval excluding zero and fallback is requested.",
      "Candidate covariate input/summary matching is empty.",
      "Fold construction creates a holdout with fewer than 3 observations.",
      "Fold-specific selector returns no covariates.",
      "chunkwise is TRUE and chunk_size >= number of prediction rows.",
      "nrow(newdata) * draws exceeds the package large-allocation threshold.",
      "Prediction buildings contain NA/Inf.",
      "Prediction buildings contain negative values.",
      "Aggregation group contains NA.",
      "Posterior draw matrix contains NA/Inf.",
      "No pair of held-out observations falls within a requested distance threshold.",
      "Synthetic n < 50, or chunkwise simulation uses chunk_size >= n."
    ),
    suggested_action=c(
      "Project coordinates to an appropriate CRS and use matching distance units.",
      "Impute, filter, or otherwise preprocess missing values deliberately.",
      "Combine sparse levels, simplify the hierarchy, or interpret them cautiously.",
      "Confirm the override is scientifically intended.",
      "Supply mapped buildings when exposure/population reconstruction is required.",
      "Use a smaller inner edge than outer edge.",
      "Reduce cutoff or increase the inner mesh edge.",
      "Inspect influential records and CPO diagnostics.",
      "Inspect optimization diagnostics and refit if necessary.",
      "Inspect failed candidate fits/CPO values rather than relying on automatic selection.",
      "Review the selector rule and fallback size.",
      "Check covariate names and selector inputs.",
      "Change folds/blocking or interpret fold-level metrics cautiously.",
      "Confirm an intercept/spatial-only fold model is acceptable.",
      "Use a smaller chunk_size if multiple chunks are desired.",
      "Use return_draws=FALSE and draw-wise aggregate_by for large grids.",
      "Correct exposure values; structural-zero treatment is a fallback safeguard.",
      "Correct invalid negative exposure values.",
      "Supply complete grouping labels or accept exclusion of missing groups.",
      "Inspect the source prediction draws before aggregation.",
      "Increase the distance threshold or inspect the spatial sampling geometry.",
      "Increase n for realistic tutorials or reduce chunk_size to demonstrate chunking."
    ),
    stringsAsFactors=FALSE
  )
}

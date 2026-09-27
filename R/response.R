#' Detect population response type
#'
#' Classifies a response as COUNT or PPB. Non-negative integer-valued responses
#' are classified as counts; strictly positive non-integer responses are
#' classified as people-per-building (PPB). Because rounded PPB can be
#' integer-valued, users can override detection with `response_type`.
#'
#' @param y Numeric response vector.
#' @return `"count"` or `"ppb"`.
#' @export
detect_response_type <- function(y) {
  .assert_numeric(y,"y")
  z <- y[!is.na(y)]
  if(!length(z)) stop("Response contains no observed values.",call.=FALSE)
  if(any(!is.finite(z))) stop("Response contains non-finite values.",call.=FALSE)
  if(any(z < 0)) stop("Population responses cannot be negative.",call.=FALSE)

  if(.is_integerish(z)) return("count")
  if(any(z <= 0))
    stop("PPB likelihoods require strictly positive responses.",call.=FALSE)
  "ppb"
}

#' Create a bottom-UpR model specification
#'
#' Hierarchical effects are optional. Each effect can be supplied as a column
#' name (shorthand for an IID effect) or as a list, for example
#' `list(column="source", model="iid", prediction="zero")`. This makes the
#' hierarchy adaptable to source, settlement, region-by-settlement, household,
#' facility, interviewer, district, or other user-defined effects.
#'
#' @param data Data frame.
#' @param response Response column name.
#' @param buildings Optional mapped-building exposure column.
#' @param response_type `"auto"`, `"count"`, or `"ppb"`.
#' @param likelihood `"auto"` or a supported likelihood.
#' @param covariates Character vector of fixed-effect columns.
#' @param coords Two coordinate column names for spatial models.
#' @param random_effects Backward-compatible character-vector shorthand for IID
#'   hierarchical effects.
#' @param hierarchical_effects Optional list of arbitrary hierarchical effects.
#' @param spatial Include an SPDE spatial field?
#' @param verbose Logical; print progress messages and validation warnings.
#' @return A `bottom_up_spec` object.
#' @export
bottom_up_spec <- function(data, response, buildings=NULL,
                           response_type=c("auto","count","ppb"),
                           likelihood="auto", covariates=NULL,
                           coords=c("x","y"), random_effects=NULL,
                           hierarchical_effects=NULL,
                           spatial=TRUE,
                           verbose=TRUE) {
  .bottom_up_progress(verbose,"Creating model specification")
  response_type <- match.arg(response_type)
  response <- .as_name(response)
  if(!is.null(buildings)) buildings <- .as_name(buildings)

  effects <- .normalize_hierarchical_effects(
    hierarchical_effects=hierarchical_effects,
    random_effects=random_effects
  )
  effect_cols <- if(length(effects)) vapply(effects,`[[`,character(1),"column") else character()

  needed <- unique(c(response,buildings,covariates,
                     if(spatial) coords else NULL,effect_cols))
  .validate_columns(data,needed)
  .warn_missingness(data,needed)
  if(spatial) .warn_coordinate_scale(data,coords)
  if(length(covariates)) {
    bad_cov <- !vapply(data[covariates],is.numeric,logical(1))
    if(any(bad_cov))
      stop("Covariates must currently be numeric. Non-numeric: ",
           paste(covariates[bad_cov],collapse=", "),call.=FALSE)
  }

  y <- data[[response]]
  detected <- detect_response_type(y)
  type <- if(response_type=="auto") detected else response_type
  if(response_type!="auto" && response_type!=detected)
    .bottom_up_warn(
      "response_type='",response_type,
      "' overrides automatic detection ('",detected,"')."
    )
  z <- y[!is.na(y)]

  if(type=="count" && (any(z < 0) || !.is_integerish(z)))
    stop("COUNT response must be non-negative and integer-valued.",call.=FALSE)
  if(type=="ppb" && any(z <= 0))
    stop("PPB response must be strictly positive.",call.=FALSE)
  if(type=="count" && is.null(buildings))
    .bottom_up_warn("COUNT model has no building exposure column; the model will not include a log(buildings) offset.")
  if(type=="ppb" && is.null(buildings))
    .bottom_up_warn("PPB model has no buildings column; population reconstruction from PPB will not be available automatically.")

  if(!identical(likelihood,"auto")) {
    likelihood <- .match_likelihood(likelihood)
    valid <- if(type=="count") c("poisson","nbinomial") else c("gamma","lognormal")
    if(!likelihood %in% valid)
      stop("Likelihood ",likelihood," is incompatible with response type ",type,
           ".",call.=FALSE)
  }

  .warn_sparse_groups(data,effects)
  .bottom_up_progress(
    verbose,"Model specification ready",
    paste0("response_type=",type,", likelihood=",likelihood,
           ", covariates=",length(covariates %||% character()),
           ", hierarchical_effects=",length(effects),
           ", spatial=",spatial)
  )

  structure(list(
    data=data,response=response,buildings=buildings,
    response_type=type,detected_type=detected,
    likelihood=likelihood,covariates=covariates %||% character(),
    coords=coords,
    random_effects=random_effects,
    hierarchical_effects=effects,
    spatial=spatial
  ),class="bottom_up_spec")
}

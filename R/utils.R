.require_inla <- function() {
  if (!requireNamespace("INLA", quietly = TRUE)) {
    stop(
      "bottom-UpR requires the INLA package for model fitting. ",
      "Install it from the official R-INLA repository; see README.md.",
      call. = FALSE
    )
  }
}

.assert_numeric <- function(x, name) {
  if (!is.numeric(x)) stop(name, " must be numeric.", call.=FALSE)
}

.is_integerish <- function(x, tol=sqrt(.Machine$double.eps)) {
  x <- x[is.finite(x)]
  length(x) > 0L && all(abs(x-round(x)) <= tol)
}

.safe_log <- function(x, name) {
  if (any(!is.finite(x) | x <= 0, na.rm=TRUE))
    stop(name, " must be strictly positive.", call.=FALSE)
  log(x)
}

.as_name <- function(x) {
  if (is.character(x) && length(x)==1L) return(x)
  stop("Column arguments must be supplied as single character names.",call.=FALSE)
}

.validate_columns <- function(data, cols) {
  miss <- setdiff(cols,names(data))
  if(length(miss)) stop("Missing columns: ",paste(miss,collapse=", "),call.=FALSE)
}

.match_likelihood <- function(x) {
  z <- tolower(gsub("[ _-]","",x))
  map <- c(gamma="gamma", poisson="poisson",
           negativebinomial="nbinomial", nb="nbinomial",
           nbinomial="nbinomial", lognormal="lognormal",
           lognorm="lognormal")
  if(!z %in% names(map))
    stop("Unsupported likelihood: ",x,
         ". Use Gamma, Poisson, Negative Binomial, or Lognormal.",call.=FALSE)
  unname(map[[z]])
}

.mean_log_cpo <- function(fit, y, likelihood) {
  cpo <- fit$cpo$cpo
  ok <- is.finite(cpo) & cpo > 0 & is.finite(y)
  if(!any(ok)) return(-Inf)
  ans <- log(cpo[ok])
  # Gaussian likelihood is fitted to log(y). Transform density back to y:
  # p_Y(y) = p_logY(log y) / y.
  if(likelihood=="lognormal") ans <- ans - log(y[ok])
  mean(ans)
}

.extract_fixed_names <- function(covariates) {
  if(is.null(covariates)) character() else as.character(covariates)
}


.bottom_up_progress <- function(enabled=TRUE, stage, detail=NULL) {
  if (!isTRUE(enabled)) return(invisible(NULL))
  msg <- paste0("[bottom-UpR] ",stage)
  if (!is.null(detail) && nzchar(as.character(detail)))
    msg <- paste0(msg,": ",detail)
  message(msg)
  invisible(NULL)
}

.bottom_up_warn <- function(..., call.=FALSE, immediate.=FALSE) {
  warning("[bottom-UpR] ",...,call.=call.,immediate.=immediate.)
}

.warn_coordinate_scale <- function(data,coords) {
  if(is.null(coords) || length(coords)!=2L) return(invisible(NULL))
  if(!all(coords %in% names(data))) return(invisible(NULL))
  xy <- as.matrix(data[,coords,drop=FALSE])
  z <- xy[is.finite(xy)]
  if(!length(z)) return(invisible(NULL))
  xr <- range(xy[,1],na.rm=TRUE)
  yr <- range(xy[,2],na.rm=TRUE)
  # Heuristic only: values bounded like longitude/latitude may indicate degrees.
  if(all(xr >= -180 & xr <= 180) && all(yr >= -90 & yr <= 90)) {
    .bottom_up_warn(
      "Coordinates appear to be longitude/latitude degrees. ",
      "SPDE mesh distances, block sizes, and range priors use coordinate units; ",
      "project coordinates before distance-based spatial modelling."
    )
  }
  invisible(NULL)
}

.warn_sparse_groups <- function(data,effects,min_n=3L) {
  if(!length(effects)) return(invisible(NULL))
  for(e in effects) {
    if(!e$column %in% names(data)) next
    tab <- table(data[[e$column]],useNA="no")
    if(length(tab) && any(tab < min_n)) {
      .bottom_up_warn(
        "Hierarchical effect '",e$column,
        "' contains levels with fewer than ",min_n,
        " observations; posterior estimates for sparse levels may be weakly identified."
      )
    }
  }
  invisible(NULL)
}

.warn_missingness <- function(data,cols) {
  cols <- intersect(cols,names(data))
  if(!length(cols)) return(invisible(NULL))
  nmiss <- vapply(data[cols],function(x) sum(is.na(x)),integer(1))
  bad <- names(nmiss)[nmiss>0]
  if(length(bad)) {
    .bottom_up_warn(
      "Missing values detected in: ",paste(bad,collapse=", "),
      ". INLA/model functions may drop or fail on incomplete rows; preprocess explicitly."
    )
  }
  invisible(NULL)
}

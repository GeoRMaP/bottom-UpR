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

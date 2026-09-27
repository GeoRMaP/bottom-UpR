#' Predict from a bottom-UpR model
#'
#' This lightweight prediction method returns fitted-support posterior summaries.
#' For PPB models with a building column it can convert PPB to population. For
#' national gridded prediction with a new SPDE A matrix, use the chunked workflow
#' described in the vignette; a fuller newdata method is planned for v0.2.
#'
#' @param object A `bottom_up_fit`.
#' @param population Convert PPB to population using mapped buildings?
#' @return Data frame of posterior summaries.
#' @export
predict_bottom_up <- function(object,population=TRUE) {
  if(!inherits(object,"bottom_up_fit")) stop("object must be bottom_up_fit.",call.=FALSE)
  sm <- object$inla$summary.fitted.values
  out <- data.frame(mean=sm$mean,sd=sm$sd,
                    lower=sm$`0.025quant`,upper=sm$`0.975quant`)
  if(object$likelihood=="lognormal") {
    # These are summaries on log(PPB) support; transform quantiles/median-like
    # summaries explicitly. Exact expected PPB should be obtained from posterior
    # draws when needed.
    out$mean <- exp(out$mean)
    out$lower <- exp(out$lower); out$upper <- exp(out$upper)
    out$sd <- NA_real_
  }
  if(object$spec$response_type=="ppb" && population &&
     !is.null(object$spec$buildings)) {
    b <- object$spec$data[[object$spec$buildings]]
    out[c("mean","lower","upper")] <- out[c("mean","lower","upper")] * b
  }
  out
}

#' Aggregate posterior population draws
#'
#' @param draws Cell-by-draw numeric matrix.
#' @param group Group identifier for each cell.
#' @return Administrative posterior summaries.
#' @export
aggregate_bottom_up <- function(draws,group) {
  if(!is.matrix(draws)) draws <- as.matrix(draws)
  if(nrow(draws)!=length(group)) stop("group must have one value per cell.",call.=FALSE)
  lev <- unique(group[!is.na(group)])
  ans <- lapply(lev,function(g) {
    z <- colSums(draws[group==g,,drop=FALSE])
    data.frame(group=as.character(g),mean=mean(z),sd=stats::sd(z),
               lower=stats::quantile(z,.025,names=FALSE),
               upper=stats::quantile(z,.975,names=FALSE),
               cv=stats::sd(z)/mean(z))
  })
  do.call(rbind,ans)
}

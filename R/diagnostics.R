#' Observation-level posterior predictions
#'
#' Generates joint-posterior predictions at the modelled observation locations
#' and returns observed population, posterior fitted population and residuals.
#'
#' @param object A `bottom_up_fit`.
#' @param draws Number of joint posterior samples.
#' @param seed Random seed.
#' @return Data frame of observation-level posterior results.
#' @export
observation_predictions <- function(object,draws=500,seed=123) {
  if(!inherits(object,"bottom_up_fit"))
    stop("object must be bottom_up_fit.",call.=FALSE)

  # Training data are already on the scale used for model fitting, so suppress
  # any external CV scaling when predicting back at observations.
  old_scaling <- object$scaling
  object$scaling <- NULL

  pr <- predict_bottom_up(
    object,
    newdata=object$spec$data,
    population=TRUE,
    draws=draws,
    return_draws=FALSE,
    seed=seed
  )$summary

  object$scaling <- old_scaling

  obs <- if(object$spec$response_type=="ppb" &&
            !is.null(object$spec$buildings)) {
    object$spec$data[[object$spec$response]] *
      object$spec$data[[object$spec$buildings]]
  } else {
    object$spec$data[[object$spec$response]]
  }

  data.frame(
    row_id=seq_along(obs),
    observed=obs,
    fitted=pr$mean,
    fitted_median=pr$median,
    lower=pr$lower,
    upper=pr$upper,
    residual=obs-pr$mean,
    covered95=obs>=pr$lower & obs<=pr$upper
  )
}

#' Plot observed versus fitted values
#'
#' Accepts either a fitted model or a cross-validation object.
#'
#' @param object A `bottom_up_fit` or `bottom_up_cv`.
#' @param draws Posterior samples when plotting a fitted model.
#' @param log_scale Use log1p axes.
#' @param main Plot title.
#' @param ... Additional arguments passed to code{plot}.
#' @return Invisibly returns plotted data.
#' @export
plot_observed_fitted <- function(object,draws=500,log_scale=FALSE,
                                 main=NULL,...) {
  if(inherits(object,"bottom_up_fit")) {
    d <- observation_predictions(object,draws=draws)
    x <- d$observed; y <- d$fitted
    if(is.null(main)) main <- paste("Observed vs fitted:",object$likelihood)
  } else if(inherits(object,"bottom_up_cv")) {
    d <- object$predictions
    x <- d$observed; y <- d$predicted
    if(is.null(main)) main <- paste("Observed vs CV prediction:",object$method)
  } else {
    stop("object must be bottom_up_fit or bottom_up_cv.",call.=FALSE)
  }

  if(log_scale) {
    x <- log1p(pmax(x,0)); y <- log1p(pmax(y,0))
    xlab <- "log1p(observed population)"
    ylab <- "log1p(predicted population)"
  } else {
    xlab <- "Observed population"; ylab <- "Predicted population"
  }

  graphics::plot(x,y,xlab=xlab,ylab=ylab,main=main,...)
  graphics::abline(0,1,lty=2)
  ok <- is.finite(x) & is.finite(y)
  if(sum(ok)>2) graphics::abline(stats::lm(y[ok]~x[ok]),lty=3)
  invisible(d)
}

#' Plot residual diagnostics
#'
#' @param object Fitted model.
#' @param draws Joint posterior samples for observation predictions.
#' @return Invisibly returns observation-level diagnostics.
#' @export
plot_bottom_up_residuals <- function(object,draws=500) {
  d <- observation_predictions(object,draws=draws)
  old <- graphics::par(no.readonly=TRUE)
  on.exit(graphics::par(old),add=TRUE)
  graphics::par(mfrow=c(2,2))

  graphics::plot(
    d$fitted,d$residual,
    xlab="Posterior fitted population",ylab="Observed - fitted",
    main="Residuals vs fitted"
  )
  graphics::abline(h=0,lty=2)

  graphics::hist(
    d$residual,breaks="FD",
    xlab="Residual",main="Residual distribution"
  )

  stats::qqnorm(d$residual,main="Residual Q-Q plot")
  stats::qqline(d$residual,lty=2)

  graphics::plot(
    d$observed,d$residual,
    xlab="Observed population",ylab="Observed - fitted",
    main="Residuals vs observed"
  )
  graphics::abline(h=0,lty=2)

  invisible(d)
}

#' Plot posterior fixed effects
#'
#' Forest plot of posterior means and 95 percent credible intervals.
#'
#' @param object Fitted model.
#' @param include_intercept Include intercept.
#' @param main Plot title.
#' @return Invisibly returns the fixed-effect summary.
#' @export
plot_posterior_fixed <- function(object,include_intercept=FALSE,
                                 main="Posterior fixed effects") {
  if(!inherits(object,"bottom_up_fit"))
    stop("object must be bottom_up_fit.",call.=FALSE)

  sm <- object$inla$summary.fixed
  if(!include_intercept)
    sm <- sm[!rownames(sm) %in% c("(Intercept)","Intercept"),,drop=FALSE]
  if(!nrow(sm)) {
    warning("No fixed effects to plot.")
    return(invisible(sm))
  }

  y <- seq_len(nrow(sm))
  graphics::plot(
    sm$mean,y,
    xlim=range(c(sm$`0.025quant`,sm$`0.975quant`),finite=TRUE),
    yaxt="n",ylab="",xlab="Posterior coefficient",main=main
  )
  graphics::segments(sm$`0.025quant`,y,sm$`0.975quant`,y)
  graphics::points(sm$mean,y,pch=19)
  graphics::axis(2,at=y,labels=rownames(sm),las=1)
  graphics::abline(v=0,lty=2)
  invisible(sm)
}

#' Plot posterior hyperparameters
#'
#' @param object Fitted model.
#' @param main Plot title.
#' @return Invisibly returns hyperparameter summaries.
#' @export
plot_posterior_hyperparameters <- function(object,
                                           main="Posterior hyperparameters") {
  if(!inherits(object,"bottom_up_fit"))
    stop("object must be bottom_up_fit.",call.=FALSE)
  sm <- object$inla$summary.hyperpar
  if(is.null(sm) || !nrow(sm)) return(invisible(sm))

  y <- seq_len(nrow(sm))
  graphics::plot(
    sm$mean,y,
    xlim=range(c(sm$`0.025quant`,sm$`0.975quant`),finite=TRUE),
    yaxt="n",ylab="",xlab="Posterior value",main=main
  )
  graphics::segments(sm$`0.025quant`,y,sm$`0.975quant`,y)
  graphics::points(sm$mean,y,pch=19)
  graphics::axis(2,at=y,labels=rownames(sm),las=1)
  invisible(sm)
}

#' Plot posterior spatial-field means at mesh nodes
#'
#' @param object Spatial fitted model.
#' @param main Plot title.
#' @return Invisibly returns node coordinates and posterior mean.
#' @export
plot_spatial_field <- function(object,main="Posterior spatial field") {
  if(!inherits(object,"bottom_up_fit") || is.null(object$mesh))
    stop("A spatial bottom_up_fit is required.",call.=FALSE)

  sm <- object$inla$summary.random$.bottom_spatial
  if(is.null(sm)) stop("Spatial-field summary was not found.",call.=FALSE)

  loc <- object$mesh$loc[,1:2,drop=FALSE]
  n <- min(nrow(loc),nrow(sm))
  z <- sm$mean[seq_len(n)]

  graphics::plot(
    loc[seq_len(n),1],loc[seq_len(n),2],
    pch=21,
    cex=.6 + 1.4*(z-min(z,na.rm=TRUE))/
      pmax(diff(range(z,na.rm=TRUE)),.Machine$double.eps),
    xlab="x",ylab="y",main=main
  )
  invisible(data.frame(x=loc[seq_len(n),1],y=loc[seq_len(n),2],mean=z))
}

#' Exploratory plots for bottom-up population data
#'
#' @param data Data frame.
#' @param response Population count or PPB response column.
#' @param buildings Optional building-count column.
#' @param covariates Optional covariates for a correlation image.
#' @return Invisibly returns summary information.
#' @export
plot_bottom_up_exploration <- function(data,response,buildings=NULL,
                                       covariates=NULL) {
  .validate_columns(data,c(response,buildings,covariates))
  y <- data[[response]]

  old <- graphics::par(no.readonly=TRUE)
  on.exit(graphics::par(old),add=TRUE)
  graphics::par(mfrow=c(2,2))

  graphics::hist(y,breaks="FD",main=paste("Response:",response),xlab=response)

  if(!is.null(buildings)) {
    b <- data[[buildings]]
    graphics::hist(b,breaks="FD",main=paste("Buildings:",buildings),xlab=buildings)
    graphics::plot(b,y,xlab="Buildings",ylab=response,
                   main="Response vs buildings")
    ppb <- y/b
    ppb[!is.finite(ppb)] <- NA
    graphics::hist(ppb,breaks="FD",main="People per building",xlab="PPB")
  } else if(length(covariates)) {
    k <- covariates[1]
    graphics::plot(data[[k]],y,xlab=k,ylab=response,
                   main=paste(response,"vs",k))
    graphics::boxplot(y,main="Response boxplot",ylab=response)
    graphics::plot(stats::ecdf(y),main="Response ECDF",xlab=response,ylab="F(x)")
  } else {
    graphics::boxplot(y,main="Response boxplot",ylab=response)
    graphics::plot(stats::ecdf(y),main="Response ECDF",xlab=response,ylab="F(x)")
    graphics::plot(seq_along(y),y,xlab="Observation",ylab=response,
                   main="Response by observation")
  }

  invisible(list(
    n=nrow(data),
    response_summary=summary(y),
    buildings_summary=if(is.null(buildings)) NULL else summary(data[[buildings]])
  ))
}


#' Plot posterior hierarchical effects
#'
#' @param object Fitted model.
#' @param effect Effect column, label, or internal name. If omitted, the first
#'   hierarchical effect is plotted.
#' @param main Optional title.
#' @return Invisibly returns plotted posterior summaries with original levels.
#' @export
plot_hierarchical_effects <- function(object,effect=NULL,main=NULL) {
  if(!inherits(object,"bottom_up_fit"))
    stop("object must be bottom_up_fit.",call.=FALSE)
  eff <- object$spec$hierarchical_effects
  if(!length(eff)) stop("Model has no hierarchical effects.",call.=FALSE)

  if(is.null(effect)) {
    i <- 1L
  } else {
    hit <- vapply(eff,function(e)
      any(c(e$column,e$label,e$internal)==effect),logical(1))
    if(!any(hit)) stop("Requested hierarchical effect was not found.",call.=FALSE)
    i <- which(hit)[1]
  }

  e <- eff[[i]]
  sm <- object$inla$summary.random[[e$internal]]
  if(is.null(sm)) stop("Posterior summary for effect was not found.",call.=FALSE)

  mp <- object$effect_maps[[i]]$map
  inv <- names(mp)[match(sm$ID,unname(mp))]
  inv[is.na(inv)] <- as.character(sm$ID[is.na(inv)])

  y <- seq_len(nrow(sm))
  if(is.null(main)) main <- paste("Hierarchical effect:",e$label)

  graphics::plot(
    sm$mean,y,
    xlim=range(c(sm$`0.025quant`,sm$`0.975quant`),finite=TRUE),
    yaxt="n",ylab="",xlab="Posterior effect",main=main
  )
  graphics::segments(sm$`0.025quant`,y,sm$`0.975quant`,y)
  graphics::points(sm$mean,y,pch=19)
  graphics::axis(2,at=y,labels=inv,las=1,cex.axis=.7)
  graphics::abline(v=0,lty=2)

  invisible(data.frame(level=inv,sm,row.names=NULL))
}

#' Cross-validation diagnostic plots
#'
#' @param object A `bottom_up_cv` object.
#' @return Invisibly returns the cross-validation predictions.
#' @export
plot_cv_diagnostics <- function(object) {
  if(!inherits(object,"bottom_up_cv"))
    stop("object must be bottom_up_cv.",call.=FALSE)

  d <- object$predictions
  d$residual <- d$observed-d$predicted

  old <- graphics::par(no.readonly=TRUE)
  on.exit(graphics::par(old),add=TRUE)
  graphics::par(mfrow=c(2,2))

  graphics::plot(
    d$observed,d$predicted,
    xlab="Observed population",ylab="CV predicted population",
    main=paste("Observed vs predicted:",object$method)
  )
  graphics::abline(0,1,lty=2)

  graphics::plot(
    d$predicted,d$residual,
    xlab="CV predicted population",ylab="Observed - predicted",
    main="Held-out residuals"
  )
  graphics::abline(h=0,lty=2)

  graphics::hist(
    d$residual,breaks="FD",
    xlab="Held-out residual",main="Residual distribution"
  )

  cover <- stats::aggregate(
    d$observed>=d$lower & d$observed<=d$upper,
    list(fold=d$fold),
    mean,na.rm=TRUE
  )
  graphics::barplot(
    cover$x,names.arg=cover$fold,ylim=c(0,1),
    xlab="Fold",ylab="95% interval coverage",
    main="Coverage by fold"
  )
  graphics::abline(h=.95,lty=2)

  invisible(d)
}

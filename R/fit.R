.prepare_fit_data <- function(spec, likelihood) {
  d <- spec$data
  y <- d[[spec$response]]

  if(likelihood=="lognormal") {
    if(any(y <= 0,na.rm=TRUE)) stop("Lognormal response must be positive.",call.=FALSE)
    d$.bottom_y <- log(y)
  } else d$.bottom_y <- y

  # Count models use mapped buildings as an exposure when supplied.
  if(spec$response_type=="count" && !is.null(spec$buildings)) {
    b <- d[[spec$buildings]]
    if(any(!is.finite(b) | b <= 0,na.rm=TRUE))
      stop("Count exposure must be strictly positive for observed records.",call.=FALSE)
    d$.bottom_log_exposure <- log(b)
  } else d$.bottom_log_exposure <- 0
  d
}

.make_formula <- function(spec, spatial=spec$spatial) {
  rhs <- c("1",spec$covariates)
  if(spec$response_type=="count" && !is.null(spec$buildings))
    rhs <- c(rhs,"offset(.bottom_log_exposure)")
  if(length(spec$random_effects)) {
    rhs <- c(rhs,paste0("f(",spec$random_effects,", model='iid')"))
  }
  if(spatial) rhs <- c(rhs,"f(.bottom_spatial, model=.bottom_spde)")
  stats::as.formula(paste(".bottom_y ~",paste(rhs,collapse=" + ")),
                    env=parent.frame())
}

.fit_candidate <- function(spec, likelihood, mesh=NULL,
                           prior_range=c(100,.05),prior_sigma=c(1,.05),
                           config=FALSE, compute_waic=TRUE) {
  .require_inla()
  d <- .prepare_fit_data(spec,likelihood)
  family <- switch(likelihood,
                   gamma="gamma",poisson="poisson",
                   nbinomial="nbinomial",lognormal="gaussian")

  if(spec$spatial) {
    if(is.null(mesh)) mesh <- make_bottom_up_mesh(d,spec$coords)
    spde <- INLA::inla.spde2.pcmatern(mesh,alpha=2,
                                      prior.range=prior_range,
                                      prior.sigma=prior_sigma)
    A <- INLA::inla.spde.make.A(mesh,loc=as.matrix(d[,spec$coords,drop=FALSE]))
    idx <- INLA::inla.spde.make.index(".bottom_spatial",mesh$n)
    fixed <- d
    fixed$.bottom_spatial <- NULL
    stk <- INLA::inla.stack(data=list(.bottom_y=d$.bottom_y),
                            A=list(A,1),
                            effects=list(idx,fixed),tag="est")
    dat <- INLA::inla.stack.data(stk)
    dat$.bottom_spde <- spde
    form <- .make_formula(spec,TRUE)
    fit <- INLA::inla(form,family=family,data=dat,
                      control.predictor=list(A=INLA::inla.stack.A(stk),compute=TRUE),
                      control.compute=list(cpo=TRUE,waic=compute_waic,
                                           config=config))
  } else {
    form <- .make_formula(spec,FALSE)
    fit <- INLA::inla(form,family=family,data=d,
                      control.predictor=list(compute=TRUE),
                      control.compute=list(cpo=TRUE,waic=compute_waic,
                                           config=config))
    spde <- NULL
  }
  list(inla=fit,mesh=mesh,spde=spde,likelihood=likelihood,
       family=family,formula=form)
}

#' Automatically select a supported likelihood
#'
#' `bottom-UpR` first determines the response support. COUNT responses compare
#' Poisson and negative binomial models; PPB responses compare Gamma and
#' lognormal models. Candidates are ranked by mean leave-one-out log predictive
#' density from INLA's conditional predictive ordinates (CPO). The lognormal
#' candidate is fitted as Gaussian on log(PPB), with the Jacobian adjustment
#' applied when its CPO density is returned to the original PPB scale.
#'
#' @param spec A `bottom_up_spec`.
#' @param mesh Optional common mesh.
#' @param prior_range PC prior range specification.
#' @param prior_sigma PC prior spatial-SD specification.
#' @return List containing selected likelihood and candidate scores.
#' @export
select_likelihood <- function(spec,mesh=NULL,
                              prior_range=c(100,.05),
                              prior_sigma=c(1,.05)) {
  if(!inherits(spec,"bottom_up_spec")) stop("spec must be bottom_up_spec.",call.=FALSE)
  candidates <- if(spec$response_type=="count")
    c("poisson","nbinomial") else c("gamma","lognormal")

  # All candidates share the same mesh so likelihood comparison does not also
  # change the spatial discretisation.
  if(spec$spatial && is.null(mesh))
    mesh <- make_bottom_up_mesh(spec$data,spec$coords)

  fits <- lapply(candidates,function(z)
    .fit_candidate(spec,z,mesh,prior_range,prior_sigma,config=FALSE))
  y <- spec$data[[spec$response]]
  scores <- vapply(fits,function(z).mean_log_cpo(z$inla,y,z$likelihood),numeric(1))
  tab <- data.frame(likelihood=candidates,mean_log_cpo=scores)
  tab <- tab[order(tab$mean_log_cpo,decreasing=TRUE),,drop=FALSE]
  list(selected=tab$likelihood[1],scores=tab,fits=fits,mesh=mesh)
}

#' Fit a bottom-UpR population model
#'
#' @param spec A `bottom_up_spec`.
#' @param mesh Optional mesh.
#' @param prior_range PC range prior.
#' @param prior_sigma PC spatial-SD prior.
#' @param config Enable INLA joint-posterior configuration.
#' @return A `bottom_up_fit`.
#' @export
fit_bottom_up <- function(spec,mesh=NULL,
                          prior_range=c(100,.05),
                          prior_sigma=c(1,.05),config=TRUE) {
  if(spec$likelihood=="auto") {
    sel <- select_likelihood(spec,mesh,prior_range,prior_sigma)
    likelihood <- sel$selected; mesh <- sel$mesh
  } else {
    likelihood <- spec$likelihood; sel <- NULL
  }
  z <- .fit_candidate(spec,likelihood,mesh,prior_range,prior_sigma,config=config)
  structure(c(z,list(spec=spec,selection=sel)),class="bottom_up_fit")
}

#' High-level bottom-up population model
#'
#' The user supplies a response. Integer-valued non-negative responses are
#' treated as COUNT; positive continuous responses as PPB. `bottom-UpR`
#' automatically compares the two compatible likelihoods and selects the one
#' with the highest mean CPO log predictive density. Use `response_type` or
#' `likelihood` to override automatic behaviour.
#'
#' @inheritParams bottom_up_spec
#' @param ... Additional arguments passed to `fit_bottom_up`.
#' @return A fitted `bottom_up_fit`.
#' @export
bottom_up <- function(data,response,buildings=NULL,
                      response_type=c("auto","count","ppb"),
                      likelihood="auto",covariates=NULL,
                      coords=c("x","y"),random_effects=NULL,
                      spatial=TRUE,...) {
  spec <- bottom_up_spec(data,response,buildings,response_type,
                         likelihood,covariates,coords,random_effects,spatial)
  fit_bottom_up(spec,...)
}

#' @export
print.bottom_up_fit <- function(x,...) {
  cat("bottom-UpR model\n")
  cat("  response type :",toupper(x$spec$response_type),"\n")
  cat("  likelihood    :",x$likelihood,"\n")
  cat("  spatial       :",x$spec$spatial,"\n")
  if(!is.null(x$selection)) {
    cat("\nAutomatic likelihood comparison (higher is better):\n")
    print(x$selection$scores,row.names=FALSE)
  }
  invisible(x)
}

#' @export
summary.bottom_up_fit <- function(object,...) {
  list(response_type=object$spec$response_type,
       likelihood=object$likelihood,
       likelihood_selection=if(is.null(object$selection)) NULL else object$selection$scores,
       fixed=object$inla$summary.fixed,
       hyperparameters=object$inla$summary.hyperpar,
       waic=object$inla$waic)
}

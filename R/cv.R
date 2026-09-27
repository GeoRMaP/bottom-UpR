.fit_scaler <- function(data,covariates) {
  if (!length(covariates)) return(list(center=numeric(),scale=numeric()))
  center <- vapply(data[covariates],mean,numeric(1),na.rm=TRUE)
  scalev <- vapply(data[covariates],stats::sd,numeric(1),na.rm=TRUE)
  scalev[!is.finite(scalev) | scalev==0] <- 1
  list(center=center,scale=scalev)
}

.apply_scaler <- function(data,scaler) {
  for (nm in names(scaler$center)) {
    data[[nm]] <- (data[[nm]]-scaler$center[[nm]])/scaler$scale[[nm]]
  }
  data
}

.make_random_folds <- function(n,k,seed=123) {
  set.seed(seed)
  sample(rep(seq_len(k),length.out=n))
}

.make_spatial_folds <- function(data,coords,block_size,k,seed=123) {
  xy <- as.matrix(data[,coords,drop=FALSE])
  if(any(!is.finite(xy))) stop("Spatial CV coordinates must be finite.",call.=FALSE)
  if(length(block_size)==1L) block_size <- rep(block_size,2)
  bx <- floor((xy[,1]-min(xy[,1]))/block_size[1])
  by <- floor((xy[,2]-min(xy[,2]))/block_size[2])
  block <- interaction(bx,by,drop=TRUE)
  tab <- sort(table(block),decreasing=TRUE)
  load <- rep(0,k)
  map <- stats::setNames(integer(length(tab)),names(tab))
  set.seed(seed)
  for (b in names(tab)) {
    choices <- which(load==min(load))
    j <- sample(choices,1)
    map[[b]] <- j
    load[j] <- load[j]+tab[[b]]
  }
  unname(map[as.character(block)])
}

.cv_metrics <- function(obs,pred,lower=NULL,upper=NULL) {
  ok <- is.finite(obs) & is.finite(pred)
  cal_intercept <- NA_real_
  cal_slope <- NA_real_
  if(sum(ok) >= 3L && stats::sd(pred[ok]) > 0) {
    cal <- stats::lm(obs[ok] ~ pred[ok])
    cal_intercept <- unname(stats::coef(cal)[1])
    cal_slope <- unname(stats::coef(cal)[2])
  }

  ans <- data.frame(
    n=sum(ok),
    RMSE=sqrt(mean((pred[ok]-obs[ok])^2)),
    MAE=mean(abs(pred[ok]-obs[ok])),
    RMSLE=sqrt(mean((log1p(pmax(pred[ok],0))-log1p(pmax(obs[ok],0)))^2)),
    bias=mean(pred[ok]-obs[ok]),
    total_ratio=sum(pred[ok])/sum(obs[ok]),
    calibration_intercept=cal_intercept,
    calibration_slope=cal_slope
  )
  if(!is.null(lower) && !is.null(upper)) {
    kk <- ok & is.finite(lower) & is.finite(upper)
    ans$coverage95 <- mean(obs[kk]>=lower[kk] & obs[kk]<=upper[kk])
  }
  ans
}

#' INLA credible-interval covariate selector
#'
#' This is an optional package selection rule, not a claim about any particular
#' manuscript. It fits one training-only bottom-UpR model containing all
#' candidates and retains fixed effects whose posterior credible interval
#' excludes zero.
#'
#' @param train_data Training data only.
#' @param response Response column.
#' @param buildings Optional building column.
#' @param candidate_covariates Candidate fixed effects.
#' @param response_type Response type.
#' @param likelihood Likelihood or `"auto"`.
#' @param coords Coordinate columns.
#' @param hierarchical_effects Optional hierarchy.
#' @param spatial Use SPDE in the selector model?
#' @param level Credible level.
#' @param fallback Number of strongest effects to retain if none exclude zero.
#' @param ... Passed to `bottom_up`.
#' @return Character vector of selected covariates.
#' @export
bottom_up_inla_selector <- function(train_data,response,buildings=NULL,
                                    candidate_covariates,
                                    response_type,
                                    likelihood="auto",
                                    coords=c("x","y"),
                                    hierarchical_effects=NULL,
                                    spatial=TRUE,
                                    level=.95,
                                    fallback=0L,
                                    verbose=TRUE,...) {
  .bottom_up_progress(verbose,"Starting INLA covariate selection",
                      paste0("candidates=",length(candidate_covariates)))
  if(!length(candidate_covariates)) {
    .bottom_up_warn("No candidate covariates supplied; selector returns an empty set.")
    return(character())
  }
  if(level <= 0 || level >= 1)
    stop("level must lie strictly between 0 and 1.",call.=FALSE)

  fit <- bottom_up(
    data=train_data,
    response=response,
    buildings=buildings,
    response_type=response_type,
    likelihood=likelihood,
    covariates=candidate_covariates,
    coords=coords,
    hierarchical_effects=hierarchical_effects,
    spatial=spatial,
    config=FALSE,
    verbose=verbose,
    ...
  )

  sm <- fit$inla$summary.fixed
  lo_name <- grep("quant",names(sm),value=TRUE)[1]
  hi_name <- grep("quant",names(sm),value=TRUE)[length(grep("quant",names(sm),value=TRUE))]
  if (is.na(lo_name) || is.na(hi_name))
    stop("Could not identify INLA fixed-effect credible interval columns.",call.=FALSE)

  rn <- rownames(sm)
  keep <- intersect(candidate_covariates,rn)
  if(!length(keep)) {
    .bottom_up_warn("None of the candidate covariates were found in INLA fixed-effect summaries; returning all candidates.")
    return(candidate_covariates)
  }

  selected <- keep[
    (sm[keep,lo_name] > 0 & sm[keep,hi_name] > 0) |
    (sm[keep,lo_name] < 0 & sm[keep,hi_name] < 0)
  ]

  if(!length(selected) && fallback > 0L) {
    .bottom_up_warn(
      "No candidate credible intervals excluded zero; applying fallback and retaining the strongest ",
      fallback," effect(s)."
    )
    score <- abs(sm[keep,"mean"]/sm[keep,"sd"])
    selected <- names(sort(score,decreasing=TRUE))[seq_len(min(fallback,length(score)))]
  }
  .bottom_up_progress(verbose,"INLA covariate selection complete",
                      paste0("selected=",length(selected),
                             if(length(selected)) paste0(" [",paste(selected,collapse=", "),"]") else ""))
  selected
}

#' Cross-validation for bottom-UpR
#'
#' Implements random k-fold, spatial-block, leave-one-region-out (LORO), and
#' leave-one-source-out (LOSO) style validation. LORO and LOSO are generic
#' grouped holdouts: supply the relevant grouping column through `group`.
#'
#' Covariate scaling is estimated from training data only. A user-supplied
#' `covariate_selector` is also called separately inside every training fold,
#' so fold-specific INLA selection can be performed without leakage.
#'
#' @param data Modelling data.
#' @param response Response column.
#' @param buildings Optional building exposure.
#' @param response_type `"auto"`, `"count"`, or `"ppb"`.
#' @param likelihood Likelihood or `"auto"`.
#' @param covariates Candidate covariates.
#' @param coords Projected coordinate columns.
#' @param hierarchical_effects Optional adaptable hierarchical effects.
#' @param random_effects Backward-compatible IID shorthand.
#' @param method Validation method.
#' @param folds Number of folds for random/spatial CV.
#' @param block_size Spatial block width in coordinate units.
#' @param group Grouping column for LORO/LOSO.
#' @param covariate_selector Optional training-only selector function. It must
#'   return a character vector of selected covariate names.
#' @param standardize Standardise covariates using training-fold means/SDs.
#' @param selector_args Named list passed only to `covariate_selector`.
#' @param prediction_draws Joint posterior draws per fold.
#' @param mesh Optional common mesh. If omitted, random/spatial CV constructs a
#'   full-location mesh once; grouped holdouts also reuse that geometry by
#'   default. The mesh uses locations only, not outcomes.
#' @param mesh_args List passed to `make_bottom_up_mesh`.
#' @param seed Random seed.
#' @param verbose Logical; print progress messages for fold construction,
#'   selection, fitting, prediction, and metrics.
#' @param ... Additional arguments passed to `bottom_up` and selector.
#' @return A `bottom_up_cv` object.
#' @export
bottom_up_cv <- function(data,response,buildings=NULL,
                         response_type=c("auto","count","ppb"),
                         likelihood="auto",covariates=NULL,
                         coords=c("x","y"),
                         hierarchical_effects=NULL,
                         random_effects=NULL,
                         method=c("spatial_block","random","loro","loso"),
                         folds=5L,block_size=100000,
                         group=NULL,
                         covariate_selector=NULL,
                         standardize=TRUE,
                         selector_args=list(),
                         prediction_draws=200L,
                         mesh=NULL,
                         mesh_args=list(),
                         seed=123,
                         verbose=TRUE,...) {
  response_type <- match.arg(response_type)
  method <- match.arg(method)
  .bottom_up_progress(verbose,"Starting cross-validation",
                      paste0("method=",method,", rows=",nrow(data)))
  .warn_coordinate_scale(data,coords)
  .warn_missingness(data,unique(c(response,buildings,covariates,coords,group)))

  effects <- .normalize_hierarchical_effects(
    hierarchical_effects=hierarchical_effects,
    random_effects=random_effects
  )

  fit_dots <- list(...)

  spec0 <- bottom_up_spec(
    data=data,response=response,buildings=buildings,
    response_type=response_type,likelihood=likelihood,
    covariates=covariates,coords=coords,
    hierarchical_effects=effects,spatial=TRUE,
    verbose=verbose
  )
  rt <- spec0$response_type

  if(method %in% c("random","spatial_block") && as.integer(folds) < 2L)
    stop("folds must be at least 2.",call.=FALSE)
  if(method=="spatial_block" && any(block_size <= 0))
    stop("block_size must be positive.",call.=FALSE)

  if(method=="random") {
    fold_id <- .make_random_folds(nrow(data),as.integer(folds),seed)
    fold_labels <- seq_len(as.integer(folds))
  } else if(method=="spatial_block") {
    fold_id <- .make_spatial_folds(data,coords,block_size,as.integer(folds),seed)
    fold_labels <- sort(unique(fold_id))
  } else {
    if(is.null(group)) {
      stop("For LORO/LOSO, supply the holdout column through group=.",call.=FALSE)
    }
    .validate_columns(data,group)
    fold_id <- as.character(data[[group]])
    fold_labels <- unique(fold_id[!is.na(fold_id)])
  }

  fold_sizes <- table(fold_id)
  if(any(fold_sizes < 3L))
    .bottom_up_warn("Some validation folds contain fewer than 3 held-out observations; metrics may be unstable.")
  .bottom_up_progress(verbose,"Validation folds ready",
                      paste0("folds=",length(fold_labels),
                             ", sizes=",paste(as.integer(fold_sizes),collapse="/")))

  if(is.null(mesh)) {
    mesh <- do.call(
      make_bottom_up_mesh,
      c(list(data=data,coords=coords,verbose=verbose),mesh_args)
    )
  }

  all_pred <- vector("list",length(fold_labels))
  selected_sets <- vector("list",length(fold_labels))
  likelihoods <- character(length(fold_labels))
  fold_metrics <- vector("list",length(fold_labels))

  for(i in seq_along(fold_labels)) {
    lab <- fold_labels[[i]]
    .bottom_up_progress(verbose,"Cross-validation fold",
                        paste0(i,"/",length(fold_labels)," [",lab,"]"))
    test_idx <- fold_id==lab
    train_idx <- !test_idx & !is.na(test_idx)
    train_raw <- data[train_idx,,drop=FALSE]
    test_raw <- data[test_idx,,drop=FALSE]
    if(!nrow(train_raw) || !nrow(test_raw))
      stop("A validation fold has an empty training or test partition.",call.=FALSE)
    .bottom_up_progress(verbose,"Fold partition",
                        paste0("train=",nrow(train_raw),", test=",nrow(test_raw)))

    scaler <- if(standardize) .fit_scaler(train_raw,covariates) else
      list(center=numeric(),scale=numeric())
    train <- if(standardize) .apply_scaler(train_raw,scaler) else train_raw

    selected <- covariates
    if(is.function(covariate_selector)) {
      sel_call <- c(
        list(
          train_data=train,
          response=response,
          buildings=buildings,
          candidate_covariates=covariates,
          response_type=rt,
          likelihood=likelihood,
          coords=coords,
          hierarchical_effects=effects,
          verbose=verbose
        ),
        selector_args
      )
      selected <- do.call(covariate_selector,sel_call)
      selected <- intersect(as.character(selected),covariates)
      if(!length(selected))
        .bottom_up_warn("Fold ",lab," selected zero covariates; fitting an intercept/spatial/hierarchical model.")
    }
    .bottom_up_progress(verbose,"Fold covariates",
                        if(length(selected)) paste(selected,collapse=", ") else "none")

    fit_call <- c(
      list(
        data=train,
        response=response,
        buildings=buildings,
        response_type=rt,
        likelihood=likelihood,
        covariates=selected,
        coords=coords,
        hierarchical_effects=effects,
        spatial=TRUE,
        mesh=mesh,
        config=TRUE,
        verbose=verbose
      ),
      fit_dots
    )
    fit <- do.call(bottom_up,fit_call)
    fit$scaling <- scaler

    pr <- predict_bottom_up(
      fit,
      newdata=test_raw,
      population=TRUE,
      draws=prediction_draws,
      return_draws=FALSE,
      seed=seed+i,
      verbose=verbose
    )$summary

    obs <- if(rt=="ppb" && !is.null(buildings)) {
      test_raw[[response]]*test_raw[[buildings]]
    } else {
      test_raw[[response]]
    }

    pp <- data.frame(
      row_id=which(test_idx),
      fold=as.character(lab),
      coord_x=test_raw[[coords[1]]],
      coord_y=test_raw[[coords[2]]],
      observed=obs,
      predicted=pr$mean,
      lower=pr$lower,
      median=pr$median,
      upper=pr$upper,
      selected_likelihood=fit$likelihood,
      stringsAsFactors=FALSE
    )

    all_pred[[i]] <- pp
    selected_sets[[i]] <- selected
    likelihoods[[i]] <- fit$likelihood
    fold_metrics[[i]] <- cbind(
      fold=as.character(lab),
      .cv_metrics(pp$observed,pp$predicted,pp$lower,pp$upper)
    )
    .bottom_up_progress(
      verbose,"Fold complete",
      paste0("RMSE=",signif(fold_metrics[[i]]$RMSE,4),
             ", MAE=",signif(fold_metrics[[i]]$MAE,4))
    )
  }

  predictions <- do.call(rbind,all_pred)
  metrics_by_fold <- do.call(rbind,fold_metrics)
  overall <- .cv_metrics(
    predictions$observed,predictions$predicted,
    predictions$lower,predictions$upper
  )

  if(any(!is.finite(predictions$predicted)))
    .bottom_up_warn("Non-finite held-out predictions were produced; inspect fold-level results.")
  .bottom_up_progress(
    verbose,"Cross-validation complete",
    paste0("RMSE=",signif(overall$RMSE,4),
           ", MAE=",signif(overall$MAE,4),
           ", coverage95=",if("coverage95" %in% names(overall)) signif(overall$coverage95,4) else "NA")
  )

  structure(
    list(
      method=method,
      predictions=predictions,
      metrics=overall,
      metrics_by_fold=metrics_by_fold,
      selected_covariates=stats::setNames(selected_sets,as.character(fold_labels)),
      selected_likelihoods=stats::setNames(likelihoods,as.character(fold_labels)),
      fold_id=fold_id,
      block_size=if(method=="spatial_block") block_size else NULL,
      group=group,
      mesh=mesh,
      response_type=rt
    ),
    class="bottom_up_cv"
  )
}

#' @export
print.bottom_up_cv <- function(x,...) {
  cat("bottom-UpR cross-validation\n")
  cat("  method :",x$method,"\n")
  cat("  folds  :",length(unique(x$predictions$fold)),"\n\n")
  print(x$metrics,row.names=FALSE)
  invisible(x)
}


#' Moran's I for cross-validation residuals
#'
#' Computes binary-distance-band Moran's I for geographically held-out residuals.
#' Coordinate units must match the supplied distance thresholds.
#'
#' @param object A `bottom_up_cv` object.
#' @param distances Distance thresholds in coordinate units.
#' @param verbose Logical; print progress messages and warnings.
#' @return Data frame of Moran's I statistics.
#' @export
bottom_up_residual_moran <- function(object,
                                     distances=c(50000,100000,150000),
                                     verbose=TRUE) {
  if(!inherits(object,"bottom_up_cv"))
    stop("object must be bottom_up_cv.",call.=FALSE)
  .bottom_up_progress(verbose,"Computing residual Moran's I",
                      paste0("distance bands=",paste(distances,collapse=", ")))
  if(any(!is.finite(distances)) || any(distances <= 0))
    stop("distances must contain positive finite thresholds.",call.=FALSE)

  d <- object$predictions
  ok <- is.finite(d$coord_x) & is.finite(d$coord_y) &
        is.finite(d$observed) & is.finite(d$predicted)
  if(any(!ok))
    .bottom_up_warn("Non-finite held-out records are excluded from Moran's I calculation.")
  d <- d[ok,,drop=FALSE]
  if(nrow(d) < 3L) stop("Too few valid held-out predictions.",call.=FALSE)

  xy <- as.matrix(d[,c("coord_x","coord_y")])
  dd <- as.matrix(stats::dist(xy))
  z <- d$observed-d$predicted
  z <- z-mean(z)
  denom <- sum(z^2)

  out <- lapply(distances,function(th) {
    w <- (dd > 0) & (dd <= th)
    s0 <- sum(w)
    I <- if(s0 > 0 && denom > 0)
      nrow(d)/s0 * sum(w * tcrossprod(z))/denom else NA_real_
    data.frame(distance=th,moran_I=I,n=nrow(d),links=s0)
  })
  ans <- do.call(rbind,out)
  if(any(ans$links == 0))
    .bottom_up_warn("One or more Moran distance bands contain no neighbour links.")
  .bottom_up_progress(verbose,"Residual Moran's I complete",
                      paste0("n=",nrow(d)))
  ans
}

#' Plot residual Moran's I across distance bands
#'
#' @param object A `bottom_up_cv`.
#' @param distances Distance thresholds.
#' @return Invisibly returns Moran's I table.
#' @export
plot_residual_moran <- function(object,
                                distances=c(50000,100000,150000)) {
  x <- bottom_up_residual_moran(object,distances)
  graphics::plot(
    x$distance,x$moran_I,type="b",
    xlab="Distance threshold",ylab="Held-out residual Moran's I",
    main="Residual spatial autocorrelation"
  )
  graphics::abline(h=0,lty=2)
  invisible(x)
}

.latent_vector <- function(sample) {
  z <- sample$latent
  if (is.matrix(z) || is.data.frame(z)) z <- z[,1]
  as.numeric(z)
}

.latent_names <- function(sample) {
  z <- sample$latent
  rn <- rownames(z)
  if (is.null(rn) && !is.null(names(z))) rn <- names(z)
  rn
}

.sample_named <- function(sample, exact=NULL, pattern=NULL, default=0) {
  nm <- .latent_names(sample)
  z <- .latent_vector(sample)
  if (is.null(nm)) return(default)

  if (!is.null(exact)) {
    k <- which(nm==exact | startsWith(nm,paste0(exact,":")))
    if (length(k)) return(z[k[1]])
  }

  if (!is.null(pattern)) {
    k <- grep(pattern,nm)
    if (length(k)) return(z[k])
  }

  default
}

.sample_field <- function(sample,prefix,n) {
  nm <- .latent_names(sample)
  z <- .latent_vector(sample)
  if (is.null(nm)) return(rep(0,n))
  k <- which(startsWith(nm,paste0(prefix,":")))
  if (!length(k)) return(rep(0,n))
  ids <- suppressWarnings(as.integer(sub(".*:","",nm[k])))
  out <- rep(0,n)
  ok <- is.finite(ids) & ids >= 1 & ids <= n
  out[ids[ok]] <- z[k[ok]]
  out
}

.gaussian_variance <- function(sample) {
  hp <- sample$hyperpar
  if (is.null(hp)) return(0)
  if (is.matrix(hp)) hp <- hp[,1]
  nm <- names(hp)
  if (is.null(nm)) nm <- rownames(sample$hyperpar)
  if (is.null(nm)) return(0)
  k <- grep("precision.*gaussian|gaussian.*precision|precision for.*observ",nm,
            ignore.case=TRUE)
  if (!length(k)) return(0)
  p <- as.numeric(hp[k[1]])
  if (!is.finite(p) || p <= 0) return(0)
  1/p
}

.effect_omitted <- function(e,omit_effects=NULL) {
  identical(e$prediction,"zero") ||
    (!is.null(omit_effects) &&
       any(c(e$column,e$label,e$internal) %in% omit_effects))
}

.prepare_newdata_effects <- function(object,newdata,omit_effects=NULL) {
  if (!length(object$spec$hierarchical_effects)) return(newdata)

  for (i in seq_along(object$spec$hierarchical_effects)) {
    e <- object$spec$hierarchical_effects[[i]]
    if (.effect_omitted(e,omit_effects)) {
      newdata[[e$internal]] <- NA_integer_
      next
    }
    m <- object$effect_maps[[i]]$map
    id <- unname(m[as.character(newdata[[e$column]])])
    newdata[[e$internal]] <- as.integer(id)
  }
  newdata
}

.apply_fit_scaling <- function(object,newdata) {
  sc <- object$scaling
  if (is.null(sc) || !length(sc$center)) return(newdata)
  for (nm in names(sc$center)) {
    if (!nm %in% names(newdata))
      stop("newdata is missing scaled covariate: ",nm,call.=FALSE)
    newdata[[nm]] <- (newdata[[nm]] - sc$center[[nm]]) / sc$scale[[nm]]
  }
  newdata
}

.init_aggregate_accumulators <- function(newdata,aggregate_by,n_draws) {
  if (is.null(aggregate_by) || !length(aggregate_by)) return(NULL)
  .validate_columns(newdata,aggregate_by)
  out <- vector("list",length(aggregate_by))
  names(out) <- aggregate_by
  for (g in aggregate_by) {
    lev <- unique(as.character(newdata[[g]][!is.na(newdata[[g]])]))
    out[[g]] <- list(
      levels=lev,
      draws=matrix(0,nrow=length(lev),ncol=n_draws,
                   dimnames=list(lev,NULL))
    )
  }
  out
}

.update_aggregate_accumulators <- function(acc,newdata,drawmat,rows) {
  if (is.null(acc)) return(acc)
  for (g in names(acc)) {
    grp <- as.character(newdata[[g]][rows])
    for (lev in acc[[g]]$levels) {
      k <- which(grp==lev)
      if (length(k))
        acc[[g]]$draws[lev,] <- acc[[g]]$draws[lev,] +
          colSums(drawmat[k,,drop=FALSE])
    }
  }
  acc
}

.summarise_draw_matrix <- function(x) {
  data.frame(
    mean=rowMeans(x),
    sd=apply(x,1,stats::sd),
    lower=apply(x,1,stats::quantile,probs=.025,names=FALSE),
    median=apply(x,1,stats::quantile,probs=.5,names=FALSE),
    upper=apply(x,1,stats::quantile,probs=.975,names=FALSE),
    cv=apply(x,1,stats::sd)/pmax(rowMeans(x),.Machine$double.eps)
  )
}

#' Predict from a bottom-UpR model
#'
#' With `newdata=NULL`, returns INLA fitted-support summaries. With
#' `newdata`, draws once from the joint INLA posterior and reuses the same
#' posterior draw identities across all chunks. This allows coherent 100-m
#' prediction and draw-wise administrative aggregation.
#'
#' For PPB models, expected PPB is multiplied by mapped buildings exactly once
#' when `population=TRUE`. For count models, building exposure enters through
#' the log offset and is not multiplied again. Cells with non-positive building
#' exposure are structural zeros.
#'
#' @param object A `bottom_up_fit`.
#' @param newdata Optional prediction data frame.
#' @param population Return population rather than PPB for PPB models.
#' @param draws Number of joint posterior samples.
#' @param chunk_size Number of prediction rows processed per chunk.
#' @param return_draws Retain cell-by-draw matrix. For large national grids this
#'   may require substantial memory.
#' @param aggregate_by Optional character vector of grouping columns in
#'   `newdata`. Administrative totals are accumulated draw-wise.
#' @param omit_effects Optional character vector of hierarchical effect columns,
#'   labels, or internal names to omit for this prediction target. This is useful
#'   for national grids where survey-source or EA-specific effects should not be
#'   assigned.
#' @param seed Random seed for posterior sampling.
#' @return If `newdata=NULL`, a data frame. Otherwise a
#'   `bottom_up_prediction` object with `summary`, optional `draws`, and
#'   `aggregates`.
#' @export
predict_bottom_up <- function(object,newdata=NULL,population=TRUE,
                              draws=500,chunk_size=50000,
                              return_draws=FALSE,
                              aggregate_by=NULL,
                              omit_effects=NULL,
                              seed=123) {
  if(!inherits(object,"bottom_up_fit"))
    stop("object must be bottom_up_fit.",call.=FALSE)

  if (is.null(newdata)) {
    sm <- object$inla$summary.fitted.values
    out <- data.frame(
      mean=sm$mean,
      sd=sm$sd,
      lower=sm$`0.025quant`,
      median=sm$`0.5quant`,
      upper=sm$`0.975quant`
    )
    if(object$likelihood=="lognormal") {
      out$mean <- exp(out$mean)
      out$lower <- exp(out$lower)
      out$median <- exp(out$median)
      out$upper <- exp(out$upper)
      out$sd <- NA_real_
    }
    if(object$spec$response_type=="ppb" && population &&
       !is.null(object$spec$buildings)) {
      b <- object$spec$data[[object$spec$buildings]]
      out[c("mean","lower","median","upper")] <-
        out[c("mean","lower","median","upper")] * b
    }
    return(out)
  }

  .require_inla()
  if (!is.data.frame(newdata)) newdata <- as.data.frame(newdata)

  needed <- unique(c(
    object$spec$covariates,
    if(object$spec$spatial) object$spec$coords else NULL,
    if(!is.null(object$spec$buildings)) object$spec$buildings else NULL,
    if(length(object$spec$hierarchical_effects)) {
      eff_keep <- Filter(function(e) !.effect_omitted(e,omit_effects),
                         object$spec$hierarchical_effects)
      if(length(eff_keep)) vapply(eff_keep,`[[`,character(1),"column") else NULL
    } else NULL,
    aggregate_by
  ))
  .validate_columns(newdata,needed)

  newdata <- .apply_fit_scaling(object,newdata)
  newdata <- .prepare_newdata_effects(object,newdata,omit_effects)

  n <- nrow(newdata)
  if (!n) stop("newdata has zero rows.",call.=FALSE)
  draws <- as.integer(draws)
  if (draws < 2) stop("draws must be at least 2.",call.=FALSE)
  chunk_size <- max(1L,as.integer(chunk_size))

  set.seed(seed)
  ps <- INLA::inla.posterior.sample(draws,object$inla)

  keep_draws <- if(return_draws)
    matrix(NA_real_,nrow=n,ncol=draws) else NULL

  summary_out <- matrix(
    NA_real_,nrow=n,ncol=6,
    dimnames=list(NULL,c("mean","sd","lower","median","upper","cv"))
  )

  acc <- .init_aggregate_accumulators(newdata,aggregate_by,draws)

  starts <- seq.int(1L,n,by=chunk_size)

  for (st in starts) {
    en <- min(n,st+chunk_size-1L)
    rows <- st:en
    nd <- newdata[rows,,drop=FALSE]
    nr <- nrow(nd)

    # Spatial projector for this chunk.
    A <- NULL
    if (object$spec$spatial) {
      A <- INLA::inla.spde.make.A(
        object$mesh,
        loc=as.matrix(nd[,object$spec$coords,drop=FALSE])
      )
    }

    dm <- matrix(0,nrow=nr,ncol=draws)

    structural_zero <- rep(FALSE,nr)
    if (!is.null(object$spec$buildings)) {
      b <- nd[[object$spec$buildings]]
      structural_zero <- !is.finite(b) | b <= 0
    } else {
      b <- rep(1,nr)
    }

    for (s in seq_len(draws)) {
      smp <- ps[[s]]

      intercept <- .sample_named(smp,exact=".bottom_intercept",default=NA_real_)
      if (!is.finite(intercept))
        intercept <- .sample_named(smp,exact="(Intercept)",default=NA_real_)
      if (!is.finite(intercept))
        intercept <- .sample_named(smp,exact="Intercept",default=0)
      eta <- rep(intercept,nr)

      for (nm in object$spec$covariates) {
        beta <- .sample_named(smp,exact=nm,default=0)
        eta <- eta + beta*nd[[nm]]
      }

      if (length(object$spec$hierarchical_effects)) {
        for (i in seq_along(object$spec$hierarchical_effects)) {
          e <- object$spec$hierarchical_effects[[i]]
          if (.effect_omitted(e,omit_effects)) next
          id <- nd[[e$internal]]
          field <- .sample_field(
            smp,
            e$internal,
            length(object$effect_maps[[i]]$map)
          )
          ok <- !is.na(id) & id >= 1 & id <= length(field)
          if (any(ok)) eta[ok] <- eta[ok] + field[id[ok]]
        }
      }

      if (object$spec$spatial) {
        w <- .sample_field(smp,".bottom_spatial",object$mesh$n)
        eta <- eta + as.numeric(A %*% w)
      }

      if (object$spec$response_type=="count" &&
          !is.null(object$spec$buildings)) {
        eta[!structural_zero] <-
          eta[!structural_zero] + log(b[!structural_zero])
      }

      mu <- if (object$likelihood=="lognormal") {
        exp(eta + .5*.gaussian_variance(smp))
      } else {
        exp(eta)
      }

      if (object$spec$response_type=="ppb" && population &&
          !is.null(object$spec$buildings)) {
        mu <- mu*b
      }

      mu[structural_zero] <- 0
      dm[,s] <- mu
    }

    summary_out[rows,] <- as.matrix(.summarise_draw_matrix(dm))
    if(return_draws) keep_draws[rows,] <- dm
    acc <- .update_aggregate_accumulators(acc,newdata,dm,rows)
  }

  summ <- data.frame(
    row_id=seq_len(n),
    summary_out,
    check.names=FALSE
  )

  aggregates <- NULL
  if (!is.null(acc)) {
    aggregates <- lapply(acc,function(z) {
      ans <- .summarise_draw_matrix(z$draws)
      cbind(group=rownames(z$draws),ans,row.names=NULL)
    })
  }

  structure(
    list(
      summary=summ,
      draws=keep_draws,
      aggregates=aggregates,
      n_draws=draws,
      population=population,
      response_type=object$spec$response_type,
      likelihood=object$likelihood
    ),
    class="bottom_up_prediction"
  )
}

#' Aggregate posterior population draws
#'
#' Sums draw-wise values before computing uncertainty summaries.
#'
#' @param draws Cell-by-draw numeric matrix.
#' @param group Group identifier for each cell.
#' @return Administrative posterior summaries.
#' @export
aggregate_bottom_up <- function(draws,group) {
  if(!is.matrix(draws)) draws <- as.matrix(draws)
  if(nrow(draws)!=length(group))
    stop("group must have one value per cell.",call.=FALSE)

  lev <- unique(group[!is.na(group)])
  ans <- lapply(lev,function(g) {
    z <- colSums(draws[group==g,,drop=FALSE])
    data.frame(
      group=as.character(g),
      mean=mean(z),
      sd=stats::sd(z),
      lower=stats::quantile(z,.025,names=FALSE),
      median=stats::quantile(z,.5,names=FALSE),
      upper=stats::quantile(z,.975,names=FALSE),
      cv=stats::sd(z)/pmax(mean(z),.Machine$double.eps)
    )
  })
  do.call(rbind,ans)
}

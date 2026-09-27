#' Simulate synthetic data for bottom-up population tutorials
#'
#' Generates a reproducible synthetic enumeration-area data set with mapped
#' buildings, spatial coordinates, continuous covariates, optional hierarchical
#' grouping variables, a positive PPB response, and an integer population-count
#' response. It is designed for examples, tutorials, and tests rather than as a
#' scientific data-generating model for any particular country.
#'
#' @param n Number of enumeration areas.
#' @param seed Random seed.
#' @param include_hierarchy Include source, settlement, region and EA identifiers.
#' @param extent Spatial side length in coordinate units.
#' @param verbose Logical; print progress messages.
#' @return A data frame with both code{ppb} and code{population} responses.
#' @export
simulate_bottom_up_data <- function(n=200L,seed=123,
                                    include_hierarchy=TRUE,
                                    extent=300000,
                                    verbose=TRUE) {
  .bottom_up_progress(verbose,"Simulating tutorial data",
                      paste0("n=",n,", extent=",extent,", seed=",seed))
  n <- as.integer(n)
  if(n < 20L) stop("n must be at least 20 for the tutorial simulator.",call.=FALSE)
  if(!is.finite(extent) || extent <= 0)
    stop("extent must be a positive finite number.",call.=FALSE)
  if(n < 50L)
    .bottom_up_warn("Small synthetic sample size; spatial and cross-validation examples may be unstable.")
  set.seed(seed)

  x <- stats::runif(n,0,extent)
  y <- stats::runif(n,0,extent)

  roads <- stats::rnorm(n)
  water <- stats::rnorm(n)
  market <- stats::rnorm(n)
  slope <- stats::rnorm(n)
  night_lights <- stats::rnorm(n)

  buildings <- pmax(
    1L,
    stats::rpois(
      n,
      lambda=exp(3.2 + 0.18*roads - 0.12*slope + 0.10*night_lights)
    )
  )

  source <- sample(c("A","B","C"),n,replace=TRUE,prob=c(.45,.35,.20))
  settlement <- sample(c("urban","periurban","rural"),n,replace=TRUE,
                       prob=c(.30,.25,.45))
  region <- sample(paste0("R",1:5),n,replace=TRUE)

  source_re <- c(A=0.00,B=0.12,C=-0.10)
  settlement_re <- c(urban=0.25,periurban=0.08,rural=-0.18)
  region_re <- setNames(c(-.12,.06,.10,-.04,.00),paste0("R",1:5))

  spatial_signal <- 0.25*sin(x/70000) + 0.20*cos(y/85000)

  eta_ppb <- 1.15 +
    0.22*roads -
    0.18*water +
    0.15*market -
    0.08*slope +
    0.12*night_lights +
    spatial_signal +
    source_re[source] +
    settlement_re[settlement] +
    region_re[region]

  ppb_mean <- exp(eta_ppb)
  gamma_shape <- 10
  ppb <- stats::rgamma(
    n,
    shape=gamma_shape,
    scale=ppb_mean/gamma_shape
  )

  count_eta <- 0.20 +
    0.20*roads -
    0.14*water +
    0.11*market -
    0.06*slope +
    0.10*night_lights +
    spatial_signal +
    source_re[source] +
    settlement_re[settlement] +
    region_re[region]

  mu_count <- buildings*exp(count_eta)
  population <- stats::rnbinom(n,size=12,mu=mu_count)

  out <- data.frame(
    ea_id=sprintf("EA%04d",seq_len(n)),
    x_m=x,
    y_m=y,
    buildings=as.integer(buildings),
    roads=roads,
    water=water,
    market=market,
    slope=slope,
    night_lights=night_lights,
    ppb=ppb,
    population_from_ppb=as.integer(round(ppb*buildings)),
    population_count=as.integer(population),
    population=as.integer(population),
    stringsAsFactors=FALSE
  )

  if(include_hierarchy) {
    out$source <- source
    out$settlement <- settlement
    out$region <- region
    out$region_settlement <- interaction(
      out$region,out$settlement,drop=TRUE
    )
  }

  .bottom_up_progress(
    verbose,"Synthetic tutorial data ready",
    paste0("rows=",nrow(out),", hierarchy=",include_hierarchy,
           ", total_buildings=",sum(out$buildings))
  )
  out
}

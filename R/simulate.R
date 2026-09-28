#' Simulate synthetic data for bottom-up population tutorials
#'
#' Generates a reproducible synthetic enumeration-area data set with mapped
#' buildings, spatial coordinates, continuous covariates, optional hierarchical
#' grouping variables, a positive PPB response, and an integer population-count
#' response. Simulation can be performed in one block or in chunks. With the
#' same seed, both modes use the same variable-wise random-number stream and are
#' intended to return identical synthetic data.
#'
#' @param n Number of enumeration areas.
#' @param seed Random seed.
#' @param include_hierarchy Include source, settlement, region and EA identifiers.
#' @param extent Spatial side length in coordinate units.
#' @param chunkwise Logical; generate each variable in chunks.
#' @param chunk_size Number of rows per simulation chunk when
#'   code{chunkwise=TRUE}.
#' @param verbose Logical; print progress messages.
#' @return A data frame with both code{ppb} and code{population} responses.
#' @export
simulate_bottom_up_data <- function(n=200L,seed=123,
                                    include_hierarchy=TRUE,
                                    extent=300000,
                                    chunkwise=FALSE,
                                    chunk_size=50000L,
                                    verbose=TRUE) {
  n <- as.integer(n)
  chunk_size <- max(1L,as.integer(chunk_size))

  .bottom_up_progress(
    verbose,
    "Simulating tutorial data",
    paste0(
      "n=",n,", extent=",extent,", seed=",seed,
      ", chunkwise=",isTRUE(chunkwise),
      if(isTRUE(chunkwise)) paste0(", chunk_size=",chunk_size) else ""
    )
  )

  if(n < 20L)
    stop("n must be at least 20 for the tutorial simulator.",call.=FALSE)
  if(!is.finite(extent) || extent <= 0)
    stop("extent must be a positive finite number.",call.=FALSE)
  if(n < 50L)
    .bottom_up_warn(
      "Small synthetic sample size; spatial and cross-validation examples may be unstable."
    )
  if(isTRUE(chunkwise) && chunk_size >= n)
    .bottom_up_warn(
      "chunkwise=TRUE but chunk_size is at least n; simulation will use one chunk."
    )

  draw_chunks <- function(FUN) {
    if(!isTRUE(chunkwise)) return(FUN(n))
    starts <- seq.int(1L,n,by=chunk_size)
    ans <- vector("list",length(starts))
    for(i in seq_along(starts)) {
      st <- starts[[i]]
      k <- min(chunk_size,n-st+1L)
      .bottom_up_progress(
        verbose,
        "Simulating chunk",
        paste0(i,"/",length(starts)," (rows ",st,"-",st+k-1L,")")
      )
      ans[[i]] <- FUN(k)
    }
    unlist(ans,use.names=FALSE)
  }

  categorical <- function(k,values,prob=NULL) {
    if(is.null(prob)) prob <- rep(1/length(values),length(values))
    prob <- prob/sum(prob)
    u <- stats::runif(k)
    values[findInterval(u,c(0,cumsum(prob)),rightmost.closed=TRUE)]
  }

  set.seed(seed)

  x <- draw_chunks(function(k) stats::runif(k,0,extent))
  y <- draw_chunks(function(k) stats::runif(k,0,extent))

  roads <- draw_chunks(function(k) stats::rnorm(k))
  water <- draw_chunks(function(k) stats::rnorm(k))
  market <- draw_chunks(function(k) stats::rnorm(k))
  slope <- draw_chunks(function(k) stats::rnorm(k))
  night_lights <- draw_chunks(function(k) stats::rnorm(k))

  building_lambda <- exp(
    3.2 + 0.18*roads - 0.12*slope + 0.10*night_lights
  )
  pos <- 1L
  buildings <- integer(n)
  starts <- if(isTRUE(chunkwise)) seq.int(1L,n,by=chunk_size) else 1L
  for(i in seq_along(starts)) {
    st <- starts[[i]]
    idx <- st:min(n,st+chunk_size-1L)
    buildings[idx] <- pmax(
      1L,
      stats::rpois(length(idx),lambda=building_lambda[idx])
    )
    pos <- pos + length(idx)
  }

  source <- draw_chunks(function(k)
    categorical(k,c("A","B","C"),c(.45,.35,.20)))
  settlement <- draw_chunks(function(k)
    categorical(k,c("urban","periurban","rural"),c(.30,.25,.45)))
  region <- draw_chunks(function(k)
    categorical(k,paste0("R",1:5)))

  source_re <- c(A=0.00,B=0.12,C=-0.10)
  settlement_re <- c(urban=0.25,periurban=0.08,rural=-0.18)
  region_re <- stats::setNames(c(-.12,.06,.10,-.04,.00),paste0("R",1:5))

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

  ppb <- numeric(n)
  starts <- if(isTRUE(chunkwise)) seq.int(1L,n,by=chunk_size) else 1L
  for(i in seq_along(starts)) {
    st <- starts[[i]]
    idx <- st:min(n,st+chunk_size-1L)
    ppb[idx] <- stats::rgamma(
      length(idx),
      shape=gamma_shape,
      scale=ppb_mean[idx]/gamma_shape
    )
  }

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
  population <- integer(n)
  for(i in seq_along(starts)) {
    st <- starts[[i]]
    idx <- st:min(n,st+chunk_size-1L)
    population[idx] <- stats::rnbinom(
      length(idx),size=12,mu=mu_count[idx]
    )
  }

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
    verbose,
    "Synthetic tutorial data ready",
    paste0(
      "rows=",nrow(out),
      ", hierarchy=",include_hierarchy,
      ", total_buildings=",sum(out$buildings)
    )
  )

  out
}

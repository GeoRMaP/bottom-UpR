#' Derive data-adaptive SPDE mesh parameters
#'
#' Derives mesh distances from the spatial footprint of the supplied projected
#' coordinates. This is useful when applying bottom-UpR to study areas with very
#' different geographic extents. The returned distances are in the same units
#' as the coordinate columns.
#'
#' @param data Data frame containing projected coordinates.
#' @param coords Two coordinate column names.
#' @param inner_edge_fraction Fraction of the study-area diagonal used for the
#'   inner maximum edge.
#' @param outer_edge_fraction Fraction of the study-area diagonal used for the
#'   outer maximum edge.
#' @param inner_offset_fraction Fraction of the study-area diagonal used for the
#'   inner mesh extension.
#' @param outer_offset_fraction Fraction of the study-area diagonal used for the
#'   outer mesh extension.
#' @param cutoff_fraction Fraction of the study-area diagonal used as cutoff.
#' @param min_cutoff Optional lower bound for cutoff in coordinate units.
#' @param verbose Logical; print progress messages.
#' @return A list with code{max_edge}, code{offset}, code{cutoff}, spatial
#'   ranges, diagonal length, and coordinate centre.
#' @export
bottom_up_mesh_parameters <- function(
    data,
    coords=c("x","y"),
    inner_edge_fraction=0.08,
    outer_edge_fraction=0.16,
    inner_offset_fraction=0.05,
    outer_offset_fraction=0.10,
    cutoff_fraction=0.005,
    min_cutoff=0,
    verbose=TRUE) {

  .bottom_up_progress(verbose,"Deriving adaptive mesh parameters")
  .validate_columns(data,coords)
  .warn_coordinate_scale(data,coords)
  xy <- as.matrix(data[,coords,drop=FALSE])
  if (ncol(xy) != 2L)
    stop("coords must identify exactly two coordinate columns.",call.=FALSE)

  ok <- is.finite(xy[,1]) & is.finite(xy[,2])
  xy <- xy[ok,,drop=FALSE]
  if (nrow(xy) < 3L)
    stop("At least three finite coordinate pairs are required.",call.=FALSE)

  xr <- range(xy[,1])
  yr <- range(xy[,2])
  xspan <- diff(xr)
  yspan <- diff(yr)
  diagonal <- sqrt(xspan^2 + yspan^2)

  if (!is.finite(diagonal) || diagonal <= 0)
    stop("Spatial coordinates must span a non-zero study area.",call.=FALSE)

  max_edge <- diagonal*c(inner_edge_fraction,outer_edge_fraction)
  offset <- diagonal*c(inner_offset_fraction,outer_offset_fraction)
  cutoff <- max(diagonal*cutoff_fraction,min_cutoff)

  if(any(max_edge <= 0) || any(offset < 0) || cutoff < 0)
    stop("Derived mesh parameters must be non-negative, with positive max_edge.",call.=FALSE)
  if(max_edge[1] >= max_edge[2])
    .bottom_up_warn("Adaptive inner max_edge is not smaller than outer max_edge.")
  .bottom_up_progress(
    verbose,
    "Adaptive mesh parameters ready",
    paste0("diagonal=",signif(diagonal,5),
           ", max_edge=",paste(signif(max_edge,5),collapse="/"),
           ", offset=",paste(signif(offset,5),collapse="/"),
           ", cutoff=",signif(cutoff,5))
  )

  list(
    max_edge=max_edge,
    offset=offset,
    cutoff=cutoff,
    diagonal=diagonal,
    x_range=xr,
    y_range=yr,
    centre=c(x=mean(xr),y=mean(yr)),
    fractions=list(
      inner_edge=inner_edge_fraction,
      outer_edge=outer_edge_fraction,
      inner_offset=inner_offset_fraction,
      outer_offset=outer_offset_fraction,
      cutoff=cutoff_fraction
    )
  )
}

#' Build an SPDE mesh
#'
#' Mesh parameters can either be supplied explicitly or derived automatically
#' from the geographic footprint of the input locations. Automatic mode scales
#' edge lengths, extensions, and cutoff to the diagonal extent of the supplied
#' projected coordinates, so the defaults adapt to differently sized study
#' areas while preserving the coordinate units.
#'
#' @param data Data frame containing projected coordinates.
#' @param coords Two coordinate column names.
#' @param adaptive Logical; derive unspecified mesh parameters from the spatial
#'   footprint of code{data}.
#' @param max_edge Inner/outer maximum triangle edges. If code{NULL} and
#'   code{adaptive=TRUE}, derived from code{data}.
#' @param offset Inner/outer mesh extensions. If code{NULL} and
#'   code{adaptive=TRUE}, derived from code{data}.
#' @param cutoff Minimum point separation. If code{NULL} and
#'   code{adaptive=TRUE}, derived from code{data}.
#' @param adaptive_args Named list passed to
#'   code{bottom_up_mesh_parameters()}.
#' @return An INLA mesh with the parameters used stored in the
#'   code{"bottom_up_mesh_parameters"} attribute.
#' @export
make_bottom_up_mesh <- function(
    data,
    coords=c("x","y"),
    adaptive=FALSE,
    max_edge=NULL,
    offset=NULL,
    cutoff=NULL,
    adaptive_args=list(),
    verbose=TRUE) {

  .bottom_up_progress(verbose,"Building SPDE mesh")
  .require_inla()
  .validate_columns(data,coords)
  .warn_coordinate_scale(data,coords)
  loc <- as.matrix(data[,coords,drop=FALSE])
  if(any(!is.finite(loc)))
    stop("Mesh coordinates must be finite.",call.=FALSE)

  if (adaptive) {
    auto <- do.call(
      bottom_up_mesh_parameters,
      c(list(data=data,coords=coords,verbose=verbose),adaptive_args)
    )
    if (is.null(max_edge)) max_edge <- auto$max_edge
    if (is.null(offset)) offset <- auto$offset
    if (is.null(cutoff)) cutoff <- auto$cutoff
  } else {
    if (is.null(max_edge)) max_edge <- c(50,100)
    if (is.null(offset)) offset <- c(50,100)
    if (is.null(cutoff)) cutoff <- 30
    auto <- NULL
  }

  if(length(max_edge)!=2L || length(offset)!=2L)
    stop("max_edge and offset must each have length 2.",call.=FALSE)
  if(any(max_edge<=0) || any(offset<0) || cutoff<0)
    stop("Mesh distances must be non-negative and max_edge strictly positive.",call.=FALSE)
  if(max_edge[1] >= max_edge[2])
    .bottom_up_warn("Inner max_edge is not smaller than outer max_edge; check mesh settings.")
  if(cutoff > max_edge[1])
    .bottom_up_warn("cutoff exceeds inner max_edge; this may over-thin input locations.")

  .bottom_up_progress(
    verbose,"Creating triangulation",
    paste0("n_locations=",nrow(loc),
           ", max_edge=",paste(signif(max_edge,5),collapse="/"),
           ", offset=",paste(signif(offset,5),collapse="/"),
           ", cutoff=",signif(cutoff,5))
  )

  mesh <- INLA::inla.mesh.2d(
    loc=loc,
    max.edge=max_edge,
    offset=offset,
    cutoff=cutoff
  )

  .bottom_up_progress(verbose,"SPDE mesh complete",
                      paste0("mesh nodes=",mesh$n))
  attr(mesh,"bottom_up_mesh_parameters") <- list(
    adaptive=adaptive,
    max_edge=max_edge,
    offset=offset,
    cutoff=cutoff,
    data_location=if(is.null(auto)) NULL else list(
      x_range=auto$x_range,
      y_range=auto$y_range,
      centre=auto$centre,
      diagonal=auto$diagonal
    )
  )

  mesh
}

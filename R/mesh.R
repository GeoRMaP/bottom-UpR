#' Build a default SPDE mesh
#'
#' @param data Data frame containing projected coordinates.
#' @param coords Two coordinate column names.
#' @param max_edge Inner/outer maximum triangle edge.
#' @param offset Inner/outer extension.
#' @param cutoff Minimum point separation.
#' @return An INLA mesh.
#' @export
make_bottom_up_mesh <- function(data,coords=c("x","y"),
                                max_edge=c(50,100),
                                offset=c(50,100),cutoff=30) {
  .require_inla()
  .validate_columns(data,coords)
  loc <- as.matrix(data[,coords,drop=FALSE])
  INLA::inla.mesh.2d(loc=loc,max.edge=max_edge,offset=offset,cutoff=cutoff)
}

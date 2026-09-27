#' Launch a bottom-UpR Shiny demonstration
#'
#' Opens an interactive Shiny application for exploring gridded posterior
#' population results. Two demos are bundled: a single-result grid explorer and
#' a two-model comparison app.
#'
#' @param app Which demo to launch: code{"grid"} or code{"compare"}.
#' @param launch.browser Passed to code{shiny::runApp()}.
#' @param verbose Logical; print launch progress.
#' @param ... Additional arguments passed to code{shiny::runApp()}.
#' @return Invisibly returns the value from code{shiny::runApp()}.
#' @export
run_bottom_up_shiny <- function(app=c("grid","compare"),
                                launch.browser=interactive(),
                                verbose=TRUE,...) {
  app <- match.arg(app)
  .bottom_up_progress(verbose,"Launching Shiny demo",paste0("app=",app))
  if(!requireNamespace("shiny",quietly=TRUE))
    stop("Package 'shiny' is required. Install it with install.packages('shiny').",
         call.=FALSE)
  if(!requireNamespace("leaflet",quietly=TRUE))
    stop("Package 'leaflet' is required. Install it with install.packages('leaflet').",
         call.=FALSE)

  app_dir <- system.file(
    "shiny",
    if(app=="grid") "grid-results" else "model-comparison",
    package="bottom.UpR"
  )
  if(!nzchar(app_dir))
    stop("Bundled Shiny application was not found in the installed package.",
         call.=FALSE)

  shiny::runApp(
    appDir=app_dir,
    launch.browser=launch.browser,
    ...
  )
}

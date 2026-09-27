\name{run_bottom_up_shiny}
\alias{run_bottom_up_shiny}
\title{Launch a bottom-UpR Shiny demonstration}
\usage{
run_bottom_up_shiny(app = c("grid", "compare"),
  launch.browser = interactive(), verbose = TRUE, ...)
}
\arguments{
\item{app}{Bundled demonstration to launch: \code{"grid"} for a gridded posterior explorer or \code{"compare"} for a two-model comparison app.}
\item{launch.browser}{Passed to \code{shiny::runApp}.}
\item{verbose}{Logical; print application-launch progress.}
\item{...}{Additional arguments passed to \code{shiny::runApp}.}
}
\value{The value returned by \code{shiny::runApp}, invisibly.}
\description{Launches bundled interactive demonstrations for exploring synthetic gridded posterior population outputs.}

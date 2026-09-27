\name{plot_hierarchical_effects}
\alias{plot_hierarchical_effects}
\title{Posterior hierarchical-effect plot}
\usage{plot_hierarchical_effects(object, effect = NULL, main = NULL)}
\arguments{
\item{object}{A fitted \code{bottom_up_fit}.}
\item{effect}{Effect column, label, or internal name; first effect if omitted.}
\item{main}{Optional plot title.}
}
\value{Posterior effect summaries with original levels, invisibly.}
\description{Plots posterior means and 95 percent credible intervals for a selected hierarchical effect.}

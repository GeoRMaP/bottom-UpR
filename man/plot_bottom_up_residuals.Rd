\name{plot_bottom_up_residuals}
\alias{plot_bottom_up_residuals}
\title{Residual diagnostic plots}
\usage{plot_bottom_up_residuals(object, draws = 500)}
\arguments{
\item{object}{A fitted \code{bottom_up_fit}.}
\item{draws}{Number of posterior draws.}
}
\value{Observation-level diagnostics, invisibly.}
\description{Plots residuals versus fitted and observed values, a histogram and a Q-Q plot.}

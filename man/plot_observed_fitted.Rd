\name{plot_observed_fitted}
\alias{plot_observed_fitted}
\title{Observed versus fitted or cross-validated predictions}
\usage{plot_observed_fitted(object, draws = 500, log_scale = FALSE, main = NULL, ...)}
\arguments{
\item{object}{A \code{bottom_up_fit} or \code{bottom_up_cv} object.}
\item{draws}{Posterior draws for fitted-model predictions.}
\item{log_scale}{Use log1p axes.}
\item{main}{Optional plot title.}
\item{...}{Additional graphical arguments passed to \code{plot}.}
}
\value{Plotted data, invisibly.}
\description{Scatter plot with identity line and fitted trend.}

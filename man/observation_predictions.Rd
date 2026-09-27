\name{observation_predictions}
\alias{observation_predictions}
\title{Observation-level posterior predictions}
\usage{observation_predictions(object, draws = 500, seed = 123, verbose = TRUE)}
\arguments{
\item{object}{A fitted \code{bottom_up_fit}.}
\item{draws}{Number of joint posterior draws.}
\item{seed}{Random seed.}
\item{verbose}{Logical; print observation-prediction progress.}
}
\value{Observation-level observed, fitted, interval, residual and coverage results.}
\description{Generates joint-posterior results at modelled observation locations.}

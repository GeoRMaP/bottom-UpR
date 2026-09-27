\name{aggregate_bottom_up}
\alias{aggregate_bottom_up}
\title{Aggregate posterior population draws}
\usage{aggregate_bottom_up(draws, group)}
\arguments{
\item{draws}{Numeric matrix with cells/units in rows and posterior draws in columns.}
\item{group}{Grouping value for each row of \code{draws}.}
}
\value{A data frame of draw-wise aggregate posterior summaries.}
\description{Sums posterior draws by group before computing means, uncertainty intervals and coefficients of variation.}

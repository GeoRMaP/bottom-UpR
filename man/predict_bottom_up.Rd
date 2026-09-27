\name{predict_bottom_up}
\alias{predict_bottom_up}
\title{Predict from a bottom-UpR model}
\usage{
predict_bottom_up(object, newdata = NULL, population = TRUE,
  draws = 500, chunkwise = TRUE, chunk_size = 50000,
  return_draws = FALSE,
  aggregate_by = NULL, omit_effects = NULL,
  marginalize_effects = NULL, seed = 123, verbose = TRUE)
}
\arguments{
\item{object}{A fitted \code{bottom_up_fit}.}
\item{newdata}{Optional prediction data frame; fitted-support summaries if absent.}
\item{population}{For PPB models, convert PPB to population using mapped buildings.}
\item{draws}{Number of joint posterior samples.}
\item{chunkwise}{If \code{TRUE}, process prediction rows in chunks. If \code{FALSE}, process all rows in one block.}
\item{chunk_size}{Rows processed per fine-grid chunk when \code{chunkwise=TRUE}.}
\item{return_draws}{Retain the cell-by-draw matrix.}
\item{aggregate_by}{Optional grouping columns for draw-wise aggregation.}
\item{omit_effects}{Hierarchical effects whose contributions are set to zero.}
\item{marginalize_effects}{IID hierarchical effects integrated over their sampled variance.}
\item{seed}{Random seed used for posterior sampling.}
\item{verbose}{Logical; print posterior-sampling and chunk progress plus warnings.}
}
\value{A fitted-support data frame or a \code{bottom_up_prediction} object.}
\description{Performs synchronized joint-posterior fine-grid prediction, structural-zero handling and coherent administrative aggregation.}

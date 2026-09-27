\name{predict_bottom_up}
\alias{predict_bottom_up}
\title{Predict from a bottom-UpR model}
\usage{
predict_bottom_up(object, newdata = NULL, population = TRUE,
  draws = 500, chunk_size = 50000, return_draws = FALSE,
  aggregate_by = NULL, omit_effects = NULL, seed = 123)
}
\description{Returns fitted-support summaries or chunked joint-posterior newdata predictions with coherent draw-wise administrative aggregation.}

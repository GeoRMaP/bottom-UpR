test_that("non-spatial count fit and joint posterior prediction work", {
  skip_if_not_installed("INLA")
  set.seed(11)

  n <- 30L
  d <- data.frame(
    y = NA_integer_,
    buildings = sample(5:20, n, replace = TRUE),
    x1 = rnorm(n),
    x = seq_len(n),
    yy = seq_len(n),
    source = rep(c("A", "B", "C"), length.out = n),
    admin = rep(c("north", "south"), each = n / 2)
  )
  mu <- d$buildings * exp(0.2 + 0.25 * d$x1)
  d$y <- stats::rpois(n, mu)

  fit <- bottom_up(
    data = d,
    response = "y",
    buildings = "buildings",
    response_type = "count",
    likelihood = "poisson",
    covariates = "x1",
    coords = c("x", "yy"),
    hierarchical_effects = list("source"),
    spatial = FALSE,
    config = TRUE
  )

  expect_s3_class(fit, "bottom_up_fit")
  expect_equal(fit$likelihood, "poisson")

  pr <- predict_bottom_up(
    fit,
    newdata = d,
    draws = 10,
    chunk_size = 10,
    aggregate_by = "admin",
    seed = 21
  )

  expect_s3_class(pr, "bottom_up_prediction")
  expect_equal(nrow(pr$summary), n)
  expect_true(all(is.finite(pr$summary$mean)))
  expect_true("admin" %in% names(pr$aggregates))
  expect_equal(nrow(pr$aggregates$admin), 2)
})

test_that("spatial count fit can project joint posterior samples", {
  skip_if_not_installed("INLA")
  set.seed(12)

  g <- expand.grid(x = 0:4, y = 0:4)
  n <- nrow(g)
  d <- data.frame(
    y_count = NA_integer_,
    buildings = sample(5:15, n, replace = TRUE),
    cov = rnorm(n),
    x = g$x,
    y = g$y
  )
  d$y_count <- stats::rpois(
    n,
    d$buildings * exp(0.1 + 0.15 * d$cov)
  )

  mesh <- make_bottom_up_mesh(
    d,
    coords = c("x", "y"),
    max_edge = c(1.5, 3),
    offset = c(1, 2),
    cutoff = 0.1
  )

  fit <- bottom_up(
    data = d,
    response = "y_count",
    buildings = "buildings",
    response_type = "count",
    likelihood = "poisson",
    covariates = "cov",
    coords = c("x", "y"),
    spatial = TRUE,
    mesh = mesh,
    prior_range = c(2, 0.05),
    prior_sigma = c(1, 0.05),
    config = TRUE
  )

  expect_s3_class(fit, "bottom_up_fit")

  pr <- predict_bottom_up(
    fit,
    newdata = d,
    draws = 5,
    chunk_size = 10,
    seed = 22
  )

  expect_equal(nrow(pr$summary), n)
  expect_true(all(is.finite(pr$summary$mean)))
})


test_that("chunk-wise and single-block posterior prediction agree", {
  skip_if_not_installed("INLA")
  set.seed(31)

  n <- 24L
  d <- data.frame(
    y = NA_integer_,
    buildings = sample(5:15, n, replace = TRUE),
    x1 = rnorm(n),
    x = seq_len(n),
    yy = seq_len(n)
  )
  d$y <- stats::rpois(
    n,
    d$buildings * exp(0.2 + 0.2 * d$x1)
  )

  fit <- bottom_up(
    data = d,
    response = "y",
    buildings = "buildings",
    response_type = "count",
    likelihood = "poisson",
    covariates = "x1",
    coords = c("x", "yy"),
    spatial = FALSE,
    config = TRUE
  )

  a <- predict_bottom_up(
    fit,
    newdata = d,
    draws = 8,
    chunkwise = TRUE,
    chunk_size = 5,
    seed = 77
  )

  b <- predict_bottom_up(
    fit,
    newdata = d,
    draws = 8,
    chunkwise = FALSE,
    seed = 77
  )

  expect_equal(a$summary, b$summary, tolerance = 1e-10)
  expect_true(a$chunkwise)
  expect_false(b$chunkwise)
})


test_that("automatic likelihood selection reuses selected fit", {
  skip_if_not_installed("INLA")
  set.seed(91)

  n <- 28L
  d <- data.frame(
    y = rpois(n, 8),
    buildings = rep(5L, n),
    x1 = rnorm(n),
    x = seq_len(n),
    yy = seq_len(n)
  )

  fit <- bottom_up(
    data = d,
    response = "y",
    buildings = "buildings",
    response_type = "count",
    likelihood = "auto",
    covariates = "x1",
    coords = c("x", "yy"),
    spatial = FALSE,
    config = FALSE,
    verbose = FALSE
  )

  expect_s3_class(fit, "bottom_up_fit")
  expect_true(fit$likelihood %in% c("poisson", "nbinomial"))
  expect_true(fit$likelihood %in% names(fit$selection$fits))
  expect_identical(
    fit$inla,
    fit$selection$fits[[fit$likelihood]]$inla
  )
})

test_that("separate-selection CV does not call supplied selector", {
  skip_if_not_installed("INLA")
  set.seed(92)

  n <- 30L
  d <- data.frame(
    y = rpois(n, 10),
    buildings = rep(5L, n),
    x1 = rnorm(n),
    x2 = rnorm(n),
    x = seq_len(n),
    yy = seq_len(n)
  )

  selector_that_must_not_run <- function(...) {
    stop("selector should not be called")
  }

  expect_warning(
    cv <- bottom_up_cv(
      data = d,
      response = "y",
      buildings = "buildings",
      response_type = "count",
      likelihood = "poisson",
      covariates = c("x1", "x2"),
      fixed_covariates = "x1",
      selection_mode = "separate",
      covariate_selector = selector_that_must_not_run,
      coords = c("x", "yy"),
      spatial = FALSE,
      method = "random",
      folds = 2,
      prediction_draws = 5,
      verbose = FALSE,
      seed = 5
    ),
    "selection_mode='separate'"
  )

  expect_s3_class(cv, "bottom_up_cv")
  expect_equal(cv$selection_mode, "separate")
  expect_equal(cv$fixed_covariates, "x1")
  expect_true(all(vapply(cv$selected_covariates, identical, logical(1), "x1")))
})

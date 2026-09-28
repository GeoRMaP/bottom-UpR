test_that("hierarchical effects are optional and adaptable", {
  d <- data.frame(
    y=c(1,2,3,4),
    x=1:4,
    yy=1:4,
    source=c("A","A","B","B")
  )

  s0 <- bottom_up_spec(
    d,response="y",coords=c("x","yy"),
    response_type="count",spatial=FALSE
  )
  expect_length(s0$hierarchical_effects,0)

  s1 <- bottom_up_spec(
    d,response="y",coords=c("x","yy"),
    response_type="count",spatial=FALSE,
    hierarchical_effects=list(
      list(column="source",model="iid",prediction="known")
    )
  )
  expect_length(s1$hierarchical_effects,1)
  expect_equal(s1$hierarchical_effects[[1]]$column,"source")
})

test_that("random and spatial fold constructors return valid fold IDs", {
  f <- bottom.UpR:::.make_random_folds(20,5,seed=1)
  expect_length(f,20)
  expect_equal(sort(unique(f)),1:5)

  d <- data.frame(x=rep(seq(0,300,100),each=4),
                  y=rep(seq(0,300,100),4))
  sf <- bottom.UpR:::.make_spatial_folds(
    d,c("x","y"),block_size=100,k=4,seed=1
  )
  expect_length(sf,nrow(d))
  expect_true(all(sf %in% 1:4))
})

test_that("effect omission can be requested by column or label", {
  e <- list(column="source",label="survey source",
            internal=".bottom_re1",prediction="known")
  expect_false(bottom.UpR:::.effect_omitted(e,NULL))
  expect_true(bottom.UpR:::.effect_omitted(e,"source"))
  expect_true(bottom.UpR:::.effect_omitted(e,"survey source"))
})


test_that("adaptive mesh parameters scale with spatial footprint", {
  d1 <- data.frame(
    x = c(0, 100, 0, 100),
    y = c(0, 0, 100, 100)
  )
  d2 <- data.frame(
    x = d1$x * 10,
    y = d1$y * 10
  )

  p1 <- bottom_up_mesh_parameters(d1, coords = c("x", "y"))
  p2 <- bottom_up_mesh_parameters(d2, coords = c("x", "y"))

  expect_equal(p2$diagonal / p1$diagonal, 10)
  expect_equal(p2$max_edge / p1$max_edge, c(10, 10))
  expect_equal(p2$offset / p1$offset, c(10, 10))
  expect_equal(p2$cutoff / p1$cutoff, 10)
})


test_that("bundled Shiny demo directories are installed", {
  grid_app <- system.file(
    "shiny", "grid-results", "app.R",
    package = "bottom.UpR"
  )
  compare_app <- system.file(
    "shiny", "model-comparison", "app.R",
    package = "bottom.UpR"
  )

  expect_true(nzchar(grid_app))
  expect_true(nzchar(compare_app))
  expect_true(file.exists(grid_app))
  expect_true(file.exists(compare_app))
})


test_that("progress helper obeys verbose and warnings are prefixed", {
  expect_message(
    bottom.UpR:::.bottom_up_progress(TRUE, "Testing", "progress"),
    "\\[bottom-UpR\\] Testing: progress"
  )

  expect_silent(
    bottom.UpR:::.bottom_up_progress(FALSE, "Testing", "progress")
  )

  expect_warning(
    bottom.UpR:::.bottom_up_warn("test warning"),
    "\\[bottom-UpR\\] test warning"
  )
})

test_that("synthetic simulator supports quiet mode", {
  expect_silent(
    d <- simulate_bottom_up_data(
      n = 50,
      seed = 1,
      verbose = FALSE
    )
  )
  expect_equal(nrow(d), 50)
})


test_that("full and chunk-wise synthetic simulation agree", {
  a <- simulate_bottom_up_data(
    n = 80,
    seed = 42,
    chunkwise = FALSE,
    verbose = FALSE
  )

  b <- simulate_bottom_up_data(
    n = 80,
    seed = 42,
    chunkwise = TRUE,
    chunk_size = 17,
    verbose = FALSE
  )

  expect_identical(a, b)
})

test_that("warning catalogue contains required fields", {
  x <- bottom_up_warning_catalogue()

  expect_true(is.data.frame(x))
  expect_true(all(
    c("category", "message", "trigger", "suggested_action") %in% names(x)
  ))
  expect_gt(nrow(x), 10)
  expect_true(all(nzchar(x$message)))
})

test_that("synthetic warning demo is installed", {
  p <- system.file(
    "examples",
    "synthetic_warning_demo.R",
    package = "bottom.UpR"
  )

  expect_true(nzchar(p))
  expect_true(file.exists(p))
})

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

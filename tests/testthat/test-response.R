test_that("response support is detected", {
  expect_equal(detect_response_type(c(0,1,2,5)), "count")
  expect_equal(detect_response_type(c(.5,1.2,3.4)), "ppb")
  expect_error(detect_response_type(c(-1,2)))
})

test_that("likelihood compatibility is checked", {
  d <- data.frame(y=c(1,2,3),x=1:3,yy=1:3)
  expect_error(bottom_up_spec(d,"y",response_type="count",likelihood="gamma",
                              coords=c("x","yy"),spatial=FALSE))
})

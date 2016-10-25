context("sfsR matrix input")

test_that("the input is a matrix", {
  exmtx=matrix(sample(2,40,replace=T),8,5)[,1]
  expect_error(sfsR(exmtx))
})

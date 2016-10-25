context("nucdivR input")

test_that("input of 2 matrices", {
  exmtx=matrix(sample(2,80,replace=T),ncol=10,nrow=8)
  exmtx1=matrix(sample(2,80,replace=T),ncol=10,nrow=8)
  expect_error(nucdivR(exmtx))
  expect_error(nucdivR(exmtx,exmtx1[,1]))
})

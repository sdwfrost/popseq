context("ehhRs calculation")

test_that("focal site is within range and outputdimesion", {
  exmtx=matrix(sample(2,40,replace=T),8,5)
  expect_error(ehhsR(exmtx,6))
  expect_equal(dim(ehhsR(exmtx,2)),NULL)
})


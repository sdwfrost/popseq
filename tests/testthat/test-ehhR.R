context("ehhR calculation")

test_that("the focal site is within the haplotype matrix", {
  exmtx=matrix(sample(2,40,replace=T),8,5)
  expect_error(ehhR(exmtx,6))
  expect_equal(dim(ehhR(exmtx,2))[1],2)
  })


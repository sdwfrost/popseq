context("flogliksingle input")

test_that("rownames are required", {
  freq=cbind(D2=c(0,0,57,0),D3=c(0,0,44,0),D4=c(3,0,45,0))
  expect_error(flogliksingle(c(0,1,1,1,1,1),freq,1/3*10^(-5),1))
})

test_that("6 parameters", {
  freq=cbind(D2=c(0,0,57,0),D3=c(0,0,44,0),D4=c(3,0,45,0))
  rownames(freq)=c("a","c","g","t")
  expect_error(flogliksingle(c(0,1,1,1,1),freq,1/3*10^(-5),1))
})

test_that("approxfun_rcpp works in the same manner as approxfun", {
  x <- seq(0, 1, length.out = 11)
  y <- 2*seq(0, 1, length.out = 11)

  expect_equal(approxfun_rcpp(x, y, 0.55)[1], 1.1)
})

test_that("approxfun_rcpp has rule = 2 behaviour as per approxfun at the boundaries", {
  x <- seq(0, 1, length.out = 11)
  y <- 2*seq(0, 1, length.out = 11)

  expect_equal(approxfun_rcpp(x, y, c(1.1, 1.2))[1:2], c(2, 2))
  expect_equal(approxfun_rcpp(x, y, c(-0.1, -0.2))[1:2], c(0, 0))
})

test_that("approxfun_rcpp throws an error with unsorted x", {
  x <- seq(0, 1, length.out = 11)
  x <- x * -1
  y <- x * -2

  expect_error(approxfun_rcpp(x, y, 0.55), "Input vector x is not sorted in ascending order")
})

test_that("approxExtrap_rcpp works as Hmisc::approxExtrap does", {
  x <- seq(0, 1, length.out = 11)
  y <- 2*seq(0, 1, length.out = 11)

  expect_equal(approxExtrap_rcpp(x, y, c(-0.1, 0.5, 1.1))[1:3], Hmisc::approxExtrap(x, y, c(-0.1, 0.5, 1.1))$y)
})

test_that("binary_cfdr produces expected result on whole simulated data set", {
  set.seed(2)
  n <- 1000
  n1p <- 50
  zp <- c(rnorm(n1p, sd=5), rnorm(n-n1p, sd=1))
  p <- 2*pnorm(-abs(zp))

  # generate q
  q <- rbinom(n, 1, 0.1)

  group <- c(rep("A", n/2), rep("B", n/2)) 

  expect_equal(digest::digest(binary_cfdr(p, q, group)), "0a7a84c07f88c3747eaa2f113d882a27")
})

# NB: would normally be overkill to test this, the evaluation of the if(length(ind) > length(p)*0.5) expression but I'd like to check my fix of the original error works before I submit a PR
test_that("binary_cfdr executes with low p-q correlation", {
  set.seed(2)
  n <- 10000
  n1p <- 500
  zp <- c(rnorm(n1p, sd=5), rnorm(n-n1p, sd=1))
  p <- 2*pnorm(-abs(zp))

  # generate q
  q <- rbinom(n, 1, 0.1)

  group <- c(rep("A", n/2), rep("B", n/2)) 

  expect_equal(digest::digest(binary_cfdr(p, q, group)), "1ceb99eb0fbcf560d4a340bd7fc67455")
})

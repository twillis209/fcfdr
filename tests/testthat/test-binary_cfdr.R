test_that("approxfun_rcpp works in the same manner as approxfun", {
  x <- seq(0, 1, length.out = 11)
  y <- 2*seq(0, 1, length.out = 11)

  expect_equal(approxfun_rcpp(x, y, 0.55)[1], 1.1)
})

test_that("approxfun_rcpp has rule = 2 behaviour as per approxfun at the boundaries", {
  x <- seq(0, 1, length.out = 11)
  y <- 2*seq(0, 1, length.out = 11)

  expect_equal(approxfun_rcpp(x, y, c(0, 1))[1:2], c(0, 2))
  expect_equal(approxfun_rcpp(x, y, c(1.1, 1.2))[1:2], c(2, 2))
  expect_equal(approxfun_rcpp(x, y, c(-0.1, -0.2))[1:2], c(0, 0))
})

test_that("approxExtrap_rcpp works as Hmisc::approxExtrap does", {
  x <- seq(0, 1, length.out = 11)
  y <- 2*seq(0, 1, length.out = 11)

  expect_equal(approxExtrap_rcpp(x, y, c(-0.1, 0, 0.5, 1.1))[1:4], Hmisc::approxExtrap(x, y, c(-0.1, 0, 0.5, 1.1))$y)
})

test_that("approxExtrap_rcpp works as Hmisc::approxExtrap does on a shuffled data set", {
  x <- seq(0, 1, length.out = 11)
  y <- 2*seq(0, 1, length.out = 11)

  set.seed(2)
  reidx <- sample(1:11)

  x <- x[reidx]
  y <- y[reidx]

  expect_equal(approxExtrap_rcpp(x, y, c(-0.1, 0, 0.5, 1.1))[1:4], Hmisc::approxExtrap(x, y, c(-0.1, 0, 0.5, 1.1))$y)
})

# TODO unsorted inputs for approxExtrap_rcpp

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

test_that("per_group_binary_cfdr_rcpp works like per_group_binary_cfdr", {
  set.seed(2)
  n <- 100
  n1p <- 5
  zp <- c(rnorm(n1p, sd=5), rnorm(n-n1p, sd=1))
  p <- 2*pnorm(-abs(zp))

  # generate q
  q <- rbinom(n, 1, 0.1)

  group <- c(rep("A", n/2), rep("B", n/2))

  logx=seq(log10(min(p)),log10(max(p)),length.out=1000)
  x=c(exp(logx),1)

  expect_equal(per_group_binary_cfdr_rcpp(p[1:(n/2)], q[1:(n/2)], p[((n/2)+1):n], q[((n/2)+1):n], x)[1:5000], per_group_binary_cfdr(p[1:(n/2)], q[1:(n/2)], p[((n/2)+1):n], q[((n/2)+1):n], x))
})

test_that("binary_cfdr_rcpp works like binary_cfdr", {
  set.seed(2)
  n <- 10000
  n1p <- 500
  zp <- c(rnorm(n1p, sd=5), rnorm(n-n1p, sd=1))
  p <- 2*pnorm(-abs(zp))

  # generate q
  q <- rbinom(n, 1, 0.1)

  group <- c(rep("A", n/2), rep("B", n/2)) 

  expect_equal(binary_cfdr_cpp(p, q, group), binary_cfdr(p, q, group))
})

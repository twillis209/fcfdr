test_that("approxfun_rcpp works in the same manner as approxfun with sorted input", {
  x <- seq(0, 1, length.out = 11)
  y <- 2*seq(0, 1, length.out = 11)

  xout <- c(0.2, 0.55, 0.7, 0.8)

  expect_equal(approxfun_rcpp(x, y, xout)[1:4], approxfun(x, y)(xout))
})

test_that("approxfun_rcpp works in the same manner as approxfun with unsorted input", {
  x <- seq(0, 1, length.out = 11)
  y <- 2*seq(0, 1, length.out = 11)

  xout <- c(0.2, 0.55, 0.7, 0.8)

  set.seed(2)

  reidx <- sample(1:11)
  x <- x[reidx]
  y <- y[reidx]

  expect_equal(approxfun_rcpp(x, y, xout)[1:4], approxfun(x, y)(xout))
})

test_that("approxfun_rcpp has rule = 2 behaviour as per approxfun at the boundaries", {
  x <- seq(0, 1, length.out = 11)
  y <- 2*seq(0, 1, length.out = 11)

  expect_equal(approxfun_rcpp(x, y, c(0, 1))[1:2], c(0, 2))
  expect_equal(approxfun_rcpp(x, y, c(1.1, 1.2))[1:2], c(2, 2))
  expect_equal(approxfun_rcpp(x, y, c(-0.1, -0.2))[1:2], c(0, 0))
})

test_that("approxExtrap_rcpp works as Hmisc::approxExtrap does when interpolating", {
  x <- seq(0, 1, length.out = 11)
  y <- 2*seq(0, 1, length.out = 11)

  xout <- c(0, 1, 0.55, 1)

  expect_equal(approxExtrap_rcpp(x, y, xout)[1:4], Hmisc::approxExtrap(x, y, xout)$y)
})

test_that("approxExtrap_rcpp works as Hmisc::approxExtrap does when interpolating with shuffled data", {
  x <- seq(0, 1, length.out = 11)
  y <- 2*seq(0, 1, length.out = 11)

  set.seed(2)
  reidx <- sample(1:11)

  x <- x[reidx]
  y <- y[reidx]

  xout <- c(0, 1, 0.55, 1)

  expect_equal(approxExtrap_rcpp(x, y, xout)[1:4], Hmisc::approxExtrap(x, y, xout)$y)
})

test_that("approxExtrap_rcpp works as Hmisc::approxExtrap does", {
  x <- seq(0, 1, length.out = 11)
  y <- 2*seq(0, 1, length.out = 11)

  expect_equal(approxExtrap_rcpp(x, y, c(-0.1, 0, 0.5, 1.1))[1:4], Hmisc::approxExtrap(x, y, c(-0.1, 0, 0.5, 1.1))$y)
})

test_that("approxExtrap_rcpp works as Hmisc::approxExtrap does on unordered data set", {
  x <- seq(0, 1, length.out = 11)
  y <- 2*seq(0, 1, length.out = 11)

  set.seed(2)
  reidx <- sample(1:11)

  x <- x[reidx]
  y <- y[reidx]
  xout <- c(-0.1, 0, 0.5, 1.1)

  expect_equal(approxExtrap_rcpp(x, y, xout)[1:4], Hmisc::approxExtrap(x, y, xout)$y)
})

test_that("binary_cfdr produces expected result on whole simulated data set", {
  set.seed(2)
  n <- 10000
  n1p <- 500
  zp <- c(rnorm(n1p, sd=5), rnorm(n-n1p, sd=1))
  p <- 2*pnorm(-abs(zp))

  # generate q
  q <- rbinom(n, 1, 0.1)

  group <- c(rep("A", n/2), rep("B", n/2))

  expected_daf <- readRDS(system.file('testdata', 'binary_cfdr_test_data_set.rds', package = 'fcfdr', mustWork = T))

  expect_equal(binary_cfdr(p, q, group), expected_daf)
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
  
  p_loo <- p[1:(n/2)]
  q_loo <- q[1:(n/2)]
  ps <- p[((n/2)+1):n]
  qs <- q[((n/2)+1):n]

  expect_equal(per_group_binary_cfdr_rcpp(p_loo, q_loo, ps, qs, x)[1:50], per_group_binary_cfdr(p_loo, q_loo, ps, qs, x), tolerance = 1e-5)
})

test_that("per_group_binary_cfdr_with_ecdf works like per_group_binary_cfdr", {
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

  p_loo <- p[1:(n/2)]
  q_loo <- q[1:(n/2)]
  ps <- p[((n/2)+1):n]
  qs <- q[((n/2)+1):n]

  expect_equal(per_group_binary_cfdr_with_ecdf(p_loo, q_loo, ps, qs, x)[1:50], per_group_binary_cfdr(p_loo, q_loo, ps, qs, x), tolerance = 1e-5)
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

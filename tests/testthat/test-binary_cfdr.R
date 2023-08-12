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

# NB: would normally be overkill to test this warning condition but I'd like to check my fix of the original error works before I submit a
test_that("binary_cfdr produces expected warning for low p-q correlation and high proportion of problematic v-values", {
  set.seed(2)
  n <- 10000
  n1p <- 500
  zp <- c(rnorm(n1p, sd=5), rnorm(n-n1p, sd=1))
  p <- 2*pnorm(-abs(zp))

                                        # generate q
  q <- rbinom(n, 1, 0.1)

  group <- c(rep("A", n/2), rep("B", n/2)) 

  expect_warning(binary_cfdr(p, q, group), "p,q have low correlation and >50% of v-values may be problematic - check results")
})

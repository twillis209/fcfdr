devtools::load_all()

set.seed(2)
n <- 1000
n1p <- 5
zp <- c(rnorm(n1p, sd=5), rnorm(n-n1p, sd=1))
p <- 2*pnorm(-abs(zp))

                                        # generate q
q <- rbinom(n, 1, 0.1)

group <- c(rep("A", n/2), rep("B", n/2))

logx=seq(log10(min(p)),log10(max(p)),length.out=1000)
x=c(exp(logx),1)

res <- per_group_binary_cfdr(p[1:(n/2)], q[1:(n/2)], p[((n/2)+1):n], q[((n/2)+1):n], x, 'res.RData')
res2 <- per_group_binary_cfdr_rcpp(p[1:(n/2)], q[1:(n/2)], p[((n/2)+1):n], q[((n/2)+1):n], x, T)

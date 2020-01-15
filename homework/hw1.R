
#3.3
n <- 1000
u <- runif(n)
x <- 2/u^(1/2)
hist(x, prob = TRUE, breaks = 300, xlim = c(0,30), main = "")
y <- seq(2, 30, length.out = n)
lines(y, 8/y^3, lwd=1.5, col="red")

#3.5
n <- 1000
u <- runif(n)
x <- as.integer(u>0.1) + as.integer(u>0.3) + as.integer(u>0.5) + as.integer(u>0.7)
p <- c(.1, .2, .2, .2, .3)
y <- sample(c(0,1,2,3,4), size = n, replace = TRUE, prob = p)
round(rbind(table(x)/n, table(y)/n, p), 3)

#3.6
#proof see rmd

#3.7
nBetaAR <- function(n, a, b){
  k <- 0 #counter for accepted
  y <- numeric(n)
  c <- beta(a, b)
  while (k < n){
    u <- runif(1)
    x <- runif(1) #random variate from g
    if (x^(a-1) * (1-x)^(b-1) > u){
      #we accept x
      k <- k + 1
      y[k] <- x
    }
  }
  return (y)
}

n <- 1000
x <- nBetaAR(n, 3, 2)
hist(x, prob = TRUE, breaks = 100, main = "", xlim = c(0, 1))
y <- seq(0, 1, length.out = n)
lines(y, y^2*(1-y)/beta(3,2), lwd=1.5, col="red")

#3.12
n <- 1000
r <- 4
beta <- 2
lambda <- rgamma(n, r, beta)
y <- rexp(n, lambda)
hist(y, prob = TRUE, breaks = 100, main = "", xlim = c(0, 10))

#3.13
yy <- seq(0, 10, length.out = n)
lines(yy, 64/(2+y)^5, lwd=1.5, col="red")

#3.14
n <- 200
mu <- c(0, 1, 2)
Sigma <- matrix(c(1, -.5, .5, -.5, 1, -.5, .5, -.5, 1), nrow = 3, ncol = 3)
rmvn.Choleski <-
  function(n, mu, Sigma) {
    # generate n random vectors from MVN(mu, Sigma)
    # dimension is inferred from mu and Sigma
    d <- length(mu)
    Q <- chol(Sigma) # Choleski factorization of Sigma
    Z <- matrix(rnorm(n*d), nrow=n, ncol=d)
    X <- Z %*% Q + matrix(mu, n, d, byrow=TRUE)
    X
  }
# generate the sample
X <- rmvn.Choleski(n, mu, Sigma)
pairs(X)
print(cor(X))

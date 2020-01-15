cat("\014");rm(list=ls())

## 1-9.1
drayleigh <- function(x, sigma){
  if (any(x < 0)) return (0)
  stopifnot(sigma > 0)
  return ((x/sigma^2) * exp(-x^2 / (2*sigma^2)))
}

rrayleigh.mh <- function(sigma, m=10000){
  x <- numeric(m)
  x[1] <- rchisq(1, df=1)
  k <- 0
  u <- runif(m)
  for (i in 2:m){
    xt <- x[i-1]
    y <- rchisq(1, df = xt)
    num <- drayleigh(y, sigma) * dchisq(xt, df = y)
    den <- drayleigh(xt, sigma) * dchisq(y, df = xt)
    if (u[i] <= num/den) x[i] <- y else{
      x[i] <- xt
      k <- k+1
    }
  }
  return (list(x=x, k=k))
}

r1 <- rrayleigh.mh(sigma = 4)
r2 <- rrayleigh.mh(sigma = 2)
r1$k; r2$k

index <- 5000:5500
y1 <- r1$x[index]
plot(index, y1, type="l", main="sigma=4", ylab="x")
y2 <- r2$x[index]
plot(index, y2, type="l", main="sigma=2", ylab="x", pin=c(1,2))

## 1-9.3
m <- 10000
x <- numeric(m)
x[1] <- 0
k <- 0
u <- runif(m)
for (i in 2:m){
  xt <- x[i-1]
  y <- rnorm(1, mean=xt, sd=1)
  num <- dcauchy(y, location=0, scale=1) * dnorm(xt, mean=y, sd=1)
  den <- dcauchy(xt, location=0, scale=1) * dnorm(y, mean=xt, sd=1)
  if (u[i] <= num/den) x[i] <- y else{
    x[i] <- xt
    k <- k+1
  }
}
x1 <- x[1001:m]
plot(x1,type='l')
k
quantile(x1, probs=seq(0, 1, 0.1))
qcauchy(seq(0, 1, 0.1), location=0, scale=1)

## 1-9.4
dLaplace <- function(x){
  return (0.5 * exp(-abs(x)))
}
rw.Metropolis <- function(sigma, x0, N){
  x <- numeric(N)
  x[1] <- x0
  u <- runif(N)
  k <- 0
  for (i in 2:N) {
    y <- rnorm(1, x[i-1], sigma)
    if (u[i] <= (dLaplace(y) / dLaplace(x[i-1])))
      x[i] <- y else {
        x[i] <- x[i-1]
        k <- k + 1
      }
  }
  return(list(x=x, k=k))
}

N <- 2000
sigma <- c(.05, .5, 2, 16)

x0 <- 25
rw1 <- rw.Metropolis(sigma[1], x0, N)
rw2 <- rw.Metropolis(sigma[2], x0, N)
rw3 <- rw.Metropolis(sigma[3], x0, N)
rw4 <- rw.Metropolis(sigma[4], x0, N)
#number of candidate points rejected
print(c(rw1$k, rw2$k, rw3$k, rw4$k))
#acceptance rate
r <- (N - c(rw1$k, rw2$k, rw3$k, rw4$k))/N
print(rbind(sigma, r))
#plot
par(mfrow=c(2,2))
plot(rw1$x, type='l', main="sigma=0.05", ylab="x")
plot(rw2$x, type='l', main="sigma=0.5", ylab="x")
plot(rw3$x, type='l', main="sigma=2", ylab="x")
plot(rw4$x, type='l', main="sigma=16", ylab="x")
dev.off()

## 1-9.7
N <- 5000 #length of chain
burn<- 1000 #burn-in length
X <- matrix(0, N, 2) #the chain, a bivariate sample
rho <- 0.9 #correlation
mu1 <- mu2 <- 0
sigma1 <- sigma2 <- 1
s1 <- sqrt(1-rho^2)*sigma1
s2 <- sqrt(1-rho^2)*sigma2

# generate the chain
X[1, ] <- c(mu1, mu2) #initialize

for (i in 2:N) {
  x2 <- X[i-1, 2]
  m1 <- mu1 + rho * (x2 - mu2) * sigma1/sigma2
  X[i, 1] <- rnorm(1, m1, s1)
  x1 <- X[i, 1]
  m2 <- mu2 + rho * (x1 - mu1) * sigma2/sigma1
  X[i, 2] <- rnorm(1, m2, s2)
}

b <- burn + 1
x <- X[b:N, ]

# compare sample statistics to parameters and plot
colMeans(x); cov(x); cor(x)
par(mfrow=c(1,2))
plot(x, main="generated data", cex=.5, xlab=bquote(X[1]),
     ylab=bquote(X[2]), ylim=range(x[,2]))

# linear regression
yy <- x[,2]
xx <- x[,1]
lm1 <- lm(yy ~ xx)
summary(lm1)
res <- lm1$residuals
hist(res, main = "residuals of linear model")
dev.off()

## 2

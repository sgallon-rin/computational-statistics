#1-5.15
M <- 10000
theta.hat <- se <- numeric(3)
g <- function(x){
  exp(-x - log(1+x^2)) * (x>0) * (x<1)
}
f <- function(x){
  exp(-x) / (1 - exp(-1))
}
F.cdf <- function(x){
  (1 - exp(-x)) / (1 - exp(-1))
}
F.quantile <- function(u){
  - log(1- u * (1 - exp(-1)))
}
#using importance sampling
u <- runif(M)
x <- F.quantile(u) #inverse transform method
fg <- g(x) / f(x)
theta.hat[1] <- mean(fg)
se[1] <- sd(fg)
#using stratified importance sampling with subinterval Ia
k <- 5
m <- M/k
T1 <- var1 <- numeric(5)
a <- seq(0, 1, length.out = 6)
for (j in 1:5){
  ai <- a[j] #a_{i-1}
  aii <- a[j+1] #a_i
  fj <- function(x){f(x) / (F.cdf(j/5) - F.cdf((j-1)/5))}
  cnt <- 0
  fg <- numeric(m)
  while (cnt < m){
    u <- runif(1)
    x <- F.quantile(u)
    if ((ai < x) & (x < aii)){
      cnt <- cnt+1
      fg[cnt] <- g(x)/fj(x)
    }
  }
  T1[j] <- mean(fg)
  var1[j] <- var(fg)
}
theta.hat[2] <- sum(T1)
se[2] <- sqrt(1/m * sum(var1))
#using stratified importance sampling with subinterval Ib
T2 <- var2 <- numeric(5)
dim(a) <- length(a)
b <- apply(a, 1, F.quantile)
fj <- function(x){5 * f(x)}
for (j in 1:5){
  ai <- b[j] #a_{i-1}
  aii <- b[j+1] #a_i
  cnt <- 0
  fg <- numeric(m)
  while (cnt < m){
    u <- runif(1)
    x <- F.quantile(u)
    if ((ai < x) & (x < aii)){
      cnt <- cnt+1
      fg[cnt] <- g(x)/fj(x)
    }
  }
  T2[j] <- mean(fg)
  var2[j] <- var(fg)
}
theta.hat[3] <- sum(T2)
se[3] <- sqrt(1/m * sum(var2))
theta.hat
se

#2-6.1
n <- 20
K <- n/2 - 1
m <- 1000
mse <- matrix(0, K, 2)
trimmed.mse <- function(n, m, k){
  #MC est for k-level trimmed mean of standard Cauchy distribution
  tmean <- numeric(m)
  for (i in 1:m){
    x <- sort(rcauchy(n))
    tmean[i] <- sum(x[(k+1):(n-k)]) / (n-2*k)
  }
  theta <- mean(tmean)
  mse.hat <- mean((tmean - theta)^2)
  se.mse <- sqrt(mean((tmean - theta)^2))/sqrt(m)
  return (c(mse.hat, se.mse))
}
for (k in 1:K){
  mse[k, 1:2] <- trimmed.mse(n=n, m=m, k=k)
}
mse.table <- data.frame(seq(1,9), mse)
names(mse.table) <- c('k', 'mse of k-level trimmed mean', 'standard error')
mse.table

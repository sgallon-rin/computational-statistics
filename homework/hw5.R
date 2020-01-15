cat("\014");rm(list=ls())

## 1-6.4
n <- 20
alpha <- 0.05
calCI <- function(n, alpha){
  x <- rlnorm(n, meanlog = 0, sdlog = 1)
  y <- log(x)
  y.bar <- mean(y)
  y.sd <- sqrt(n/(n-1)*var(y))
  z <- y.sd / sqrt(n) * qt(alpha/2, df = n-1) #Lower α/2 quantile
  return (c(y.bar+z, y.bar-z))
}
m <- 1000
CI <- matrix(0, ncol = 2, nrow = m)
cnt <- 0
for (i in 1:m){
  CI[i,] <- calCI(n, alpha)
  if ((CI[i,1]<0) & (CI[i,2]>0)){
    cnt <- cnt + 1
  }
}
cnt/m

## 1-6.9
gini <- function(x){
  x <- sort(x)
  mu <- mean(x)
  n <- length(x)
  g <- 0
  for (i in 1:n){
    g <- g + (2*i-n-1) * x[i]
  }
  g <- g / (n^2 * mu)
  return (g)
}
n <- 100
m <- 1000 #replicates
#standard lognormal
G.hat <- replicate(m, expr = {
  x <- rlnorm(n, meanlog = 0, sdlog = 1)
  gini(x)
})
mean(G.hat)
median(G.hat)
quantile(G.hat, seq(0.1,1,0.1)) #deciles
hist(G.hat, main = "Gini ratio estimation for standard lognormal dist(m=1000)")
#uniform
G.hat <- replicate(m, expr = {
  x <- runif(n)
  gini(x)
})
mean(G.hat)
median(G.hat)
quantile(G.hat, seq(0.1,1,0.1)) #deciles
hist(G.hat, main = "Gini ratio estimation for uniform dist(m=1000)")
#Bernouli(0.1)
G.hat <- replicate(m, expr = {
  x <- sample(c(0,1), size=n, replace=TRUE, prob=c(0.9,0.1))
  gini(x)
})
mean(G.hat)
median(G.hat)
quantile(G.hat, seq(0.1,1,0.1)) #deciles
hist(G.hat, main = "Gini ratio estimation for Bernouli(0.1) dist(m=1000)")

x <- seq(0,1,0.01)
g0 <- gini(x)
g <- 0
n <- length(x)
for (i in 1:n){
  for (j in 1:n){
    g <- g + abs(x[i]-x[j])
  }
}
g <- g/(2*n^2*mean(x))
g
g0

## 1-6.10
#h <- function(u1, u2){
#  #as.numeric(u2<=u1)*u1 + as.numeric(u1<=u2)*u2
#  pmax(u1, u2)
#}
F_n <- function(u, x){
  n <- length(x)
  mean(x<=u)
}
h1.hat<- function(u, x){
  u*F_n(u, x) + mean(x*(x>=u))
}
set.seed(123)
G <- gini(rlnorm(1000000, meanlog = 0, sdlog = 1)) #as nominal value
n <- 50
alpha <- 0.05
calCI.G <- function(n, alpha){
  x <- rlnorm(n, meanlog = 0, sdlog = 1)
  mu.hat <- mean(x)
  G.hat <- gini(x)
  u1 <- numeric(n)
  for (i in 1:n){
    u1[i] <- 2*h1.hat(x[i], x) - (G.hat+1)*x[i]
  }
  sigma1 <- var(u1)*n/(n-1)
  sigma1 <- sqrt(sigma1)/mu.hat
  z <- qnorm(1-alpha/2) #Upper alpha/2 quantile
  z <- z*sigma1/sqrt(n)
  return (c(G.hat-z, G.hat+z))
}
m <- 1000
CI <- matrix(0, ncol = 2, nrow = m)
cnt <- 0
for (i in 1:m){
  CI[i,] <- calCI.G(n, alpha)
  if ((CI[i,1]<G) & (CI[i,2]>G)){
    cnt <- cnt + 1
  }
}
cnt/m

## 3
#(2)
#type 1
set.seed(123)
n <- 20
alpha <- 0.05
m <- 1000
calCI1 <- function(n, alpha){
  x <- rnorm(n, mean=0, sd=2)
  return ((n-1) * var(x)/qchisq(alpha, df = n-1))
}
UCL <- replicate(m, expr = {
  calCI1(n=n, alpha=alpha)
})
p1 <- mean(UCL > 4) #empirical coverage probability
w1 <- mean(UCL) #average conficence interval width
#type 2
calCI2 <- function(n, alpha){
  x <- rnorm(n, mean=0, sd=2)
  z1 <- (n-1) * var(x)/qchisq(alpha/2, df = n-1) #Lower α/2 quantile
  z2 <- (n-1) * var(x)/qchisq(1-alpha/2, df = n-1)
  return (c(z2, z1))
}
CI <- matrix(0, ncol = 2, nrow = m)
cnt <- 0
for (i in 1:m){
  CI[i,] <- calCI2(n, alpha)
  if ((CI[i,1]<4) & (CI[i,2]>4)){
    cnt <- cnt + 1
  }
}
p2 <- cnt/m #empirical coverage probability
w2 <- mean(CI[,2]-CI[,1]) #average conficence interval width
cbind(c(p1,p2),c(w1,w2))

#(3)
set.seed(123)
n <- 20
alpha <- 0.05
m <- 1000
calCI1 <- function(n, alpha){
  x <- rchisq(n, df=2)
  return ((n-1) * var(x)/qchisq(alpha, df = n-1))
}
UCL <- replicate(m, expr = {
  calCI1(n=n, alpha=alpha)
})
p1 <- mean(UCL > 4) #empirical coverage probability
w1 <- mean(UCL) #average conficence interval width
#type 2
calCI2 <- function(n, alpha){
  x <- rchisq(n, df=2)
  z1 <- (n-1) * var(x)/qchisq(alpha/2, df = n-1) #Lower α/2 quantile
  z2 <- (n-1) * var(x)/qchisq(1-alpha/2, df = n-1)
  return (c(z2, z1))
}
CI <- matrix(0, ncol = 2, nrow = m)
cnt <- 0
for (i in 1:m){
  CI[i,] <- calCI2(n, alpha)
  if ((CI[i,1]<4) & (CI[i,2]>4)){
    cnt <- cnt + 1
  }
}
p2 <- cnt/m #empirical coverage probability
w2 <- mean(CI[,2]-CI[,1]) #average conficence interval width
cbind(c(p1,p2),c(w1,w2))
cat("\014");rm(list=ls())

x <- seq(.1, 2.5, length = 10)
m <- 10000
u <- runif(m)
cdf <- numeric(length(x))
for (i in 1:length(x)) {
  g <- x[i] * exp(-(u * x[i])^2 / 2)
  cdf[i] <- mean(g) / sqrt(2 * pi) + 0.5
}
Phi <- pnorm(x)
print(round(rbind(x, cdf, Phi), 3))

x <- seq(.1, 2.5, length = 10)
m <- 10000
z <- rnorm(m)
dim(x) <- length(x)
p <- apply(x, MARGIN = 1,
           FUN = function(x, z) {mean(z <= x)}, z = z)
Phi <- pnorm(x)
print(round(rbind(x, p, Phi), 3))

#5.3
m <- 10000
x <- runif(m, min = 0, max = 0.5)
g <- exp(-x)
theta.hat <- mean(g) * 0.5
print(theta.hat)
print(1 - exp(-0.5))
#v.hat <- mean((g - mean(g))^2) / m
v.hat <- var(g) / m
print(v.hat)
z <- rexp(m, rate = 1) #sampling from exponential dist
g <- (z<=0.5) #the indicator function
theta.star <- mean(g)
print(theta.star)
#v.star <- mean((g - mean(g))^2) / m
v.star <- var(g) / m
print(v.star)

#5.4
m <- 10000
x <- seq(0.1, 0.9, 0.1)
mc.beta.cdf <- function(x, a, b, m=10000){
  z <- rbeta(m, 3, 3)
  dim(x) <- length(x)
  p <- apply(x, MARGIN = 1, FUN = function(x, z){mean(z <= x)}, z = z)
  return(p)
}
p <- mc.beta.cdf(x, 3, 3, m)
cdf <- pbeta(x, 3, 3)
print(round(rbind(p, cdf), 3))

#5.5
x <- seq(.1, 2.5, length = 10)
m <- 10000
# sample mean Monte Carlo method
u <- runif(m)
cdf1 <- numeric(length(x))
var1 <- numeric(length(x))
for (i in 1:length(x)) {
  g <- x[i] * exp(-(u * x[i])^2 / 2)
  cdf1[i] <- mean(g) / sqrt(2 * pi) + 0.5
  var1[i] <- var(g) / m
}
# "hit-or-miss" method
z <- rnorm(m)
cdf2 <- numeric(length(x))
var2 <- numeric(length(x))
for (i in 1:length(x)) {
  g <- (z <= x[i])
  cdf2[i] <- mean(g)
  var2[i] <- var(g) / m
}
#print(rbind(cdf1, cdf2))
#print(rbind(var1, var2))
# efficiency comparison
effi <- character(length(x))
improve <- numeric(length(x))
for (i in 1:length(x)){
  if (var1[i] >= var2[i]){
    effi[i] <- "theta2"
    improve[i] <- (var1[i] - var2[i]) / var1[i]
  }else{
    effi[i] <- "theta1"
    improve[i] <- (var2[i] - var1[i]) / var2[i]
  }
}
improve <- paste0(as.character(round(improve*100, 2)), "%")
Phi <- pnorm(x)
cdfs <- round(rbind(x, cdf1, cdf2, Phi), 3)
vars <- round(rbind(var1, var2), 13)
print(rbind(cdfs, vars, effi, improve))


#5.9
MC.Rayleigh <- function(x, sigma, R = 10000, antithetic = TRUE) {
  u <- runif(R/2)
  if (!antithetic) v <- runif(R/2) else
    v <- 1 - u
  u <- c(u, v)
  cdf <- numeric(length(x))
  for (i in 1:length(x)) {
    g <- x[i]^2 * u * exp(-(u * x[i])^2 / (2 * sigma^2))
    cdf[i] <- mean(g) / (sigma^2)
  }
  cdf
}
#generation by antithetic variables
x <- seq(0.1, 4, length = 10)
sigma <- 1
MC <- MC.Rayleigh(x, sigma)
print(round(rbind(x, MC), 5))
#draw a graph of Rayleigh(1) cdf
x <- seq(0, 5, length = 500)
sigma <- 1
MC <- MC.Rayleigh(x, sigma)
plot(x, MC, type = "l")
#comparison
m <- 1000
MC1 <- MC2 <- numeric(m)
x <- 0.5
for (i in 1:m) {
  MC1[i] <- MC.Rayleigh(x, sigma, R = m, antithetic = FALSE)
  MC2[i] <- MC.Rayleigh(x, sigma, R = m)
}
print(sd(MC1))
print(sd(MC2))
improve = (var(MC1) - var(MC2))/var(MC1)
str <- paste0(as.character(round(improve*100, 2)), "%")
print(paste0("The percentage reduction in variance of antithetic variables is ", str))


#5.10
h <- function(x){
  return(exp(-x) / (1+x^2))
}
func <- function(R = 10000, antithetic = TRUE){
  u <- runif(R/2)
  if (!antithetic) v <- runif(R/2) else
    v <- 1 - u
  u <- c(u, v)
  g <- h(u)
  return(mean(g))
}
#calculate with antithetic method
set.seed(123)
theta.hat <- func()
theo <- integrate(h,0,1)
print(c(theta.hat, theo$value))
#comparison of variance
m <- 1000
MC1 <- MC2 <- numeric(m)
for (i in 1:m) {
  MC1[i] <- func(R = m, antithetic = FALSE)
  MC2[i] <- func(R = m)
}
print(sd(MC1))
print(sd(MC2))
improve = (var(MC1) - var(MC2))/var(MC1)
str <- paste0(as.character(round(improve*100, 2)), "%")
print(paste0("The the approximate reduction in varianceis ", str))

#2
v.hypersphere <- function(n, r=1){
  v <- pi^(n/2) * r^n / gamma(n/2 + 1)
  return(v)
}
pi.hypersphere1 <- function(n, m=10000){
  x <- matrix(nrow = m, ncol = n)
  #r <- numeric(m)
  for (i in 1:n){
    x[, i] <- runif(m)
  }
  #for (j in 1:m){
  #  r[j] <- sum(x[j, ]^2)
  #}
  r <- apply(x^2, 1, sum)
  pi.hat <- (sum(r <= 1)/m * 2^n * gamma(n/2 + 1)) ^ (2/n)
  return(pi.hat)
}
pi.hypersphere <- function(n, m=10000){
  #calculate pi.hat from n-dimension hypersphere with snmple size m
  x <- matrix(runif(m*n), nrow = m)
  g <- apply(x^2, 1, sum)
  theta.hat <- mean(g <= 1)
  pi.hat <- (theta.hat * 2^n * gamma(n/2 + 1)) ^ (2/n)
  return(pi.hat)
}
pi.digit1 <- function(d, m=1e6, digit=4){
  set.seed(123)
  pi.hat <- pi.hypersphere(d, m)
  print(c(m, pi.hat))
  r <- abs(round(pi.hat-pi, digit))
  while(r >= 0.001){
    set.seed(123)
    m <- m + 1e6
    pi.hat <- pi.hypersphere(d, m)
    print(c(m, pi.hat))
    r <- abs(round(pi.hat-pi, digit))
  }
  #return(c(m, pi.hat))
}
pi.digit <- function(d, m=1e7, digit=4){
  #calculate the sample size needed to approximate pi to its digit-th digit,
  #starting with smaple size m
  set.seed(123)
  pi.hat <- pi.hypersphere(d, m)
  #print(c(m, pi.hat))
  r <- abs(round(pi.hat-pi, digit-1))
  while(r >= 0.001){
    set.seed(123)
    m <- m + 1e7
    pi.hat <- pi.hypersphere(d, m)
    #print(c(m, pi.hat))
    r <- abs(round(pi.hat-pi, digit-1))
  }
  return(c(m, pi.hat))
}
dmin <- 2
dmax <- 8
x <- seq(dmin, dmax)
n <- dmax - dmin + 1
#pi.hat <- matrix(nrow = n, ncol = 2)
for (d in dmin:dmax){
  sprintf("dimension: %d", d)
  pi.digit(d)
  #pi.hat[d-dmin+1,] <- pi.digit(d)
}
#pi.hat

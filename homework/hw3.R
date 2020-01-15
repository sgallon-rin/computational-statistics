#5.14
dRayleigh <- function(x, sigma){
  #density of Rayleigh(sigma) at x
  x/sigma^2 * exp(-x^2/(2*sigma^2))
}

rRayleigh <- function(x, sigma, R = 10000, antithetic = TRUE) {
  #random sample of Rayleigh(sigma) at interval x
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

s <- 1.5 #sigma of Rayleigh
m <- 1.5 #mean of normal
n <- 10000
x <- seq(0, 10, length.out = n+1)
x <- x[2:length(x)]
g0 <- function(x){x^2 / sqrt(2*pi) * exp(-x^2/2)}
plot(x, dRayleigh(x, sigma = s), type = 'l', col = 'blue')
lines(x, g0(x), col = 'red')
lines(x, dnorm(x, mean = m), col = 'black')

plot(x, g0(x)/dRayleigh(x, sigma = s), type = 'l', col = 'blue')
lines(x, g0(x)/dnorm(x, mean = m), col = 'red')

g <- function(x){x^2 / sqrt(2*pi) * exp(-x^2/2) * (x>1)}
f1 <- function(x){
  dRayleigh(x, sigma = s) * (x>1)
}
f2 <- function(x){
  dnorm(x, mean = m) * (x>1)
}

rf1 <- function(x){
  rRayleigh(x, sigma = s)
}

rf2 <- function(n){
  rnorm(n, mean = m)
}

#calculate by importance sampling
g <- function(x){x^2 / sqrt(2*pi) * exp(-x^2/2) * (x>1)}
f1 <- function(x){dnorm(x, mean = m)}
f2 <- function(x){dRayleigh(x, sigma = s)}
x1 <- rnorm(n, mean = m)
x2 <- rRayleigh(x, sigma = s)
theta1 <- mean(g(x1)/f1(x1))
theta2 <- mean(g(x2)/f2(x2))
theta1
theta2

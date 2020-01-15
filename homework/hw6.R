cat("\014");rm(list=ls())
library(boot)
library(bootstrap)

## EX-7.3
boot.t.ci <-
  function(x, B = 500, R = 100, level = .95, statistic){
    #compute the bootstrap t CI
    x <- as.matrix(x);  n <- nrow(x)
    stat <- numeric(B); se <- numeric(B)
    boot.se <- function(x, R, f) {
      #local function to compute the bootstrap
      #estimate of standard error for statistic f(x)
      x <- as.matrix(x); m <- nrow(x)
      th <- replicate(R, expr = {
        i <- sample(1:m, size = m, replace = TRUE)
        f(x[i, ])
      })
      return(sd(th))
    }
    for (b in 1:B) {
      j <- sample(1:n, size = n, replace = TRUE)
      y <- x[j, ]
      stat[b] <- statistic(y)
      se[b] <- boot.se(y, R = R, f = statistic)
    }
    stat0 <- statistic(x)
    t.stats <- (stat - stat0) / se
    se0 <- sd(stat)
    alpha <- 1 - level
    Qt <- quantile(t.stats, c(alpha/2, 1-alpha/2), type = 1)
    names(Qt) <- rev(names(Qt))
    CI <- rev(stat0 - Qt * se0)
  }

data(law, package = "bootstrap")
theta <- cor(law82$LSAT, law82$GPA)
theta.hat <- cor(law$LSAT, law$GPA)
stat <- function(dat) {
  cor(dat[, 1], dat[, 2]) }
ci <- boot.t.ci(law, statistic = stat, B=2000, R=200)
print(list(theta=theta, theta.hat=theta.hat, ci=ci))

## Ex-7.4
data(aircondit, package = "boot")
#MLE estimate for lambda
x <- aircondit$hours
n <- length(x)
logL <- function(lambda, n, x){
  #objective function: log likelihood
  n*log(lambda) - lambda*sum(x)
}
f <- function(lambda){logL(lambda,n,x)}
MLE <- optimize(f, c(0,1), maximum = TRUE)
lambda.hat <- MLE$maximum
#set up the bootstrap, sample size is n
B <- 200            #number of replicates
lambda.b <- numeric(B)     #storage for replicates
#bootstrap estimate of bias and standard error
for (b in 1:B) {
  #randomly select the indices
  i <- sample(1:n, size = n, replace = TRUE)
  x.star <- sample(x, replace = TRUE)
  f <- function(lambda){logL(lambda, n, x.star)}
  lambda.b[b] <- optimize(f, c(0,1), maximum = TRUE)$maximum
}
bias.hat <- mean(lambda.b - lambda.hat)
se.hat <- sd(lambda.b)
print(list(est=lambda.hat, bias=bias.hat, se=se.hat))

## Ex-7.6
data(scor, package = "bootstrap")
plot(scor)
#Sigma <- cor(scor)
Sigma <- cov(scor)
Sigma
#bootstrap estimate of standard error of correlation
r <- function(x, i, c1, c2) {
  #want correlation of columns c1 and c2
  cor(x[i,c1], x[i,c2])
}
boot.sd <- function(c1, c2){
  rc1c2 <- function(x, i){r(x, i, c1, c2)}
  obj <- boot(data = scor, statistic = rc1c2, R = 200)
  y <- obj$t
  sd(y)
}
se12 <- boot.sd(1, 2)
se34 <- boot.sd(3, 4)
se35 <- boot.sd(3, 5)
se45 <- boot.sd(4, 5)
print(rbind(cbind(Sigma[1,2], Sigma[3,4], Sigma[3,5], Sigma[4,5]),
            cbind(se12, se34, se35, se45)))

## Ex-7.7
data(scor, package = "bootstrap")
Sigma <- cor(scor)
e0 <- eigen(Sigma)$values
theta.hat <- e0[1]/sum(e0)
#bootstrap estimate of standard error of theta.hat
r <- function(x, i) {
  e <- eigen(cor(x[i,]))$values
  theta.b <- e[1]/sum(e)
  return(theta.b)
}
obj <- boot(data = scor, statistic = r, R = 200)
theta.b <- obj$t
bias.hat <- mean(theta.b - theta.hat)
se.hat <- sd(theta.b)
print(list(est=theta.hat, bias=bias.hat, se=se.hat))

## Proj-6D
#(Power comparison of tests of normality)
#only one loop, for epsilon=0.1, was shown in the text
#the simulation below takes several minutes to run
library(energy)
library(mlbench)

kurtosis <- function(x, normal = TRUE) {
  #computes the sample's multivariate kurtosis statistic
  #x is a matrix, with each row being an observation
  #if normal is TRUE, return
  #(b-d(d+2))/sqrt(8d(d+2)/n)
  #to test normality
  n <- nrow(x) #sample size
  d <- ncol(x) #dimension
  sigma <- cov(x) #covariance matrix
  sigma <- solve(sigma) #compute the inverse
  xbar <- apply(x, MARGIN = 2, FUN = mean)
  f <- function(v){
    #local func to compute (x_i-xbar)^T*S^{-1}*(x_i-xbar)
    v <- v - xbar
    return((v %*% sigma %*% v)^2)
  }
  b <- mean(apply(x, MARGIN = 1, FUN = f))
  if(normal){
    b <- (b-d*(d+2)) / sqrt(8*d*(d+2)/n)
  }
  return(b)
}

# initialize input and output
alpha <- .1
D <- seq(2,5) #dimension
n <- 30 #sample size
m <- 100 #try small m for a trial run
test1 <- test2 <- numeric(m)

#critical value for Mardia's skewness and kurtosis tests
cv <- qnorm(1-alpha/2)
sim <- matrix(0, 9, 3)

# estimate power
for (d in 2:10) {
  #different dimension
  for (j in 1:m) {
    x <- mlbench.twonorm(n, d=d)$x
    test1[j] <- as.integer(abs(kurtosis(x)) >= cv)
    test2[j] <- as.integer(
      mvnorm.etest(x, R=200)$p.value <= alpha)
  }
  print(c(d, mean(test1), mean(test2)))  
  sim[d-1, ] <- c(d, mean(test1), mean(test2))
}

# plot the empirical estimates of power
plot(sim[,1], sim[,2], ylim = c(0, 1), type = "l",
     xlab = "dimension", ylab = "power")
lines(sim[,1], sim[,3], lty = 2, col=2)
#lines(sim[,1], sim[,4], lty = 4, col=3)
abline(h = alpha, lty = 3)
legend("right", 1, c("kurtosis", "energy"),
       lty = c(1,2), col=1:2, inset = .02)


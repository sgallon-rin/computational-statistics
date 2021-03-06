cat("\014");rm(list=ls())

#test
a <- seq(10, 40, .1)     #sequence for plotting fits

par(mfrow = c(2, 2))    #layout for graphs

L1 <- lm(magnetic ~ chemical)
plot(chemical, magnetic, main="Linear", pch=16)
yhat1 <- L1$coef[1] + L1$coef[2] * a
lines(a, yhat1, lwd=2)

L2 <- lm(magnetic ~ chemical + I(chemical^2))
plot(chemical, magnetic, main="Quadratic", pch=16)
yhat2 <- L2$coef[1] + L2$coef[2] * a + L2$coef[3] * a^2
lines(a, yhat2, lwd=2)

L3 <- lm(log(magnetic) ~ chemical)
plot(chemical, magnetic, main="Exponential", pch=16)
logyhat3 <- L3$coef[1] + L3$coef[2] * a
yhat3 <- exp(logyhat3)
lines(a, yhat3, lwd=2)

L4 <- lm(magnetic ~ chemical + I(chemical^2) + I(chemical^3))
plot(chemical, magnetic, main="Cubic", pch=16)
yhat4 <- L4$coef[1] + L4$coef[2] * a + L4$coef[3] * a^2 + L4$coef[4] * a^3
lines(a, yhat4, lwd=2)

dev.off()


## 7.10
library(DAAG); attach(ironslag)

n <- length(magnetic)   #in DAAG ironslag
e1 <- e2 <- e3 <- e4 <- numeric(n)

# for n-fold cross validation
# fit models on leave-one-out samples
for (k in 1:n) {
  y <- magnetic[-k]
  x <- chemical[-k]
  
  J1 <- lm(y ~ x)
  yhat1 <- J1$coef[1] + J1$coef[2] * chemical[k]
  e1[k] <- magnetic[k] - yhat1
  
  J2 <- lm(y ~ x + I(x^2))
  yhat2 <- J2$coef[1] + J2$coef[2] * chemical[k] +
    J2$coef[3] * chemical[k]^2
  e2[k] <- magnetic[k] - yhat2
  
  J3 <- lm(log(y) ~ x)
  logyhat3 <- J3$coef[1] + J3$coef[2] * chemical[k]
  yhat3 <- exp(logyhat3)
  e3[k] <- magnetic[k] - yhat3
  
  J4 <- lm(y ~ x + I(x^2) + I(x^3))
  yhat4 <- J4$coef[1] + J4$coef[2] * chemical[k] +
    J4$coef[3] * chemical[k]^2 +
    J4$coef[4] * chemical[k]^3
  e4[k] <- magnetic[k] - yhat4
}
c(mean(e1^2), mean(e2^2), mean(e3^2), mean(e4^2))

# fit models on the whole sample
L1 <- lm(magnetic ~ chemical)
r1 <- summary(L1)$adj.r.squared

L2 <- lm(magnetic ~ chemical + I(chemical^2))
r2 <- summary(L2)$adj.r.squared

L3 <- lm(log(magnetic) ~ chemical)
r3 <- summary(L3)$adj.r.squared

L4 <- lm(magnetic ~ chemical + I(chemical^2) + I(chemical^3))
r4 <- summary(L4)$adj.r.squared

c(r1, r2, r3, r4)

## 7.11
library(DAAG); attach(ironslag)

n <- length(magnetic)   #in DAAG ironslag
e1 <- e2 <- e3 <- e4 <- numeric(n*(n-1))

# for (n-1)-fold cross validation
# fit models on leave-two-out samples
idx <- 1
for (k in 1:n) {
  y <- magnetic[-k]
  x <- chemical[-k]
  for (i in 1:(n-1)){
    yy <- y[-i]
    xx <- x[-i]
    
    J1 <- lm(yy ~ xx)
    yhat1 <- J1$coef[1] + J1$coef[2] * chemical[k]
    yyhat1 <- J1$coef[1] + J1$coef[2] * x[i]
    e1[idx] <- magnetic[k] - yhat1
    e1[idx+1] <- y[i] - yyhat1
    
    J2 <- lm(yy ~ xx + I(xx^2))
    yhat2 <- J2$coef[1] + J2$coef[2] * chemical[k] +
      J2$coef[3] * chemical[k]^2
    yyhat2 <- J2$coef[1] + J2$coef[2] * x[i] +
      J2$coef[3] * x[i]^2
    e2[idx] <- magnetic[k] - yhat2
    e2[idx+1] <- y[i] - yyhat2
    
    J3 <- lm(log(yy) ~ xx)
    logyhat3 <- J3$coef[1] + J3$coef[2] * chemical[k]
    yhat3 <- exp(logyhat3)
    logyyhat3 <- J3$coef[1] + J3$coef[2] * x[i]
    yyhat3 <- exp(logyyhat3)
    e3[idx] <- magnetic[k] - yhat3
    e3[idx+1] <- y[i] - yyhat3
    
    J4 <- lm(log(yy) ~ log(xx))
    logyhat4 <- J4$coef[1] + J4$coef[2] * log(chemical[k])
    yhat4 <- exp(logyhat4)
    logyyhat4 <- J4$coef[1] + J4$coef[2] * log(x[i])
    yyhat4 <- exp(logyyhat4)
    e4[idx] <- magnetic[k] - yhat4
    e4[idx+1] <- y[i] - yyhat4
    
    idx <- idx + 2
  }
}

c(mean(e1^2), mean(e2^2), mean(e3^2), mean(e4^2))

## 7.A
mu <- 20
sigma <- 10
n <- 1000 #sample size
x <- rnorm(n, mean = mu, sd = sigma)
mu.hat <- mean(x) #sample mean

# bootstrap samples for sample mean
B <- 200
mu.b <- numeric(B)
for (b in 1:B){
  i <- sample(1:n, replace = TRUE)
  x.b <- x[i]
  mu.b[b] <- mean(x.b)
}
se.hat <- sd(mu.b)

# compute CI's
alpha <- 0.05 #significance level

# standard normal bootstrap CI
z <- qnorm(alpha/2)
ci.sn <- c(mu.hat + z, mu.hat - z)
names(ci.sn) <- c(paste0(as.character(100*alpha/2),"%"), 
                  paste0(as.character(100*(1-alpha/2)),"%"))
ci.sn

# basic bootstrap CI
ci.basic <- 2*mu.hat - rev(quantile(mu.b, c(alpha/2, 1-alpha/2)))
names(ci.basic) <- rev(names(ci.basic))
ci.basic

# percentile bootstrap CI
ci.percentile <- quantile(mu.b, c(alpha/2, 1-alpha/2))
ci.percentile

# MC study
m <- 100 #size
M <- 10000 #replicates
mc.study = data.frame(rbind(ci.sn, ci.basic, ci.percentile))
names(mc.study) <- c(paste0(as.character(100*alpha/2),"%"), 
                     paste0(as.character(100*(1-alpha/2)),"%"))
mc.study$miss.left <- numeric(3)
mc.study$miss.right <- numeric(3)

for (i in 1:M){
  i <- sample(1:n, size = m, replace = TRUE)
  mu.sample <- mean(x[i])
  for (c in 1:3){
    lower <- mc.study[c,1]
    upper <- mc.study[c,2]
    if (mu.sample < lower){
      mc.study[c,3] <- mc.study[c,3] + 1
    }else if (mu.sample > upper){
      mc.study[c,4] <- mc.study[c,4] + 1
    }
  }
}

mc.study$miss.left <- mc.study$miss.left/M
mc.study$miss.right <- mc.study$miss.right/M
mc.study

## 7.B
sk <- function(x) {
  #computes the sample skewness coeff.
  xbar <- mean(x)
  m3 <- mean((x - xbar)^3)
  m2 <- mean((x - xbar)^2)
  return( m3 / m2^1.5 )
}

mc.sk <- function(x, B=200, alpha=0.05, m=100, M=10000){
  theta.hat <- sk(x)
  
  # bootstrap samples for sample skewness
  theta.b <- numeric(B)
  for (b in 1:B){
    i <- sample(1:n, replace = TRUE)
    x.b <- x[i]
    theta.b[b] <- sk(x.b)
  }
  
  # standard normal bootstrap CI
  z <- qnorm(alpha/2)
  ci.sn <- c(theta.hat + z, theta.hat - z)
  
  # basic bootstrap CI
  ci.basic <- 2*theta.hat - rev(quantile(theta.b, c(alpha/2, 1-alpha/2)))
  
  # percentile bootstrap CI
  ci.percentile <- quantile(theta.b, c(alpha/2, 1-alpha/2))
  
  # MC study
  mc.study = data.frame(rbind(ci.sn, ci.basic, ci.percentile))
  names(mc.study) <- c(paste0(as.character(100*alpha/2),"%"), 
                       paste0(as.character(100*(1-alpha/2)),"%"))
  mc.study$miss.left <- numeric(3)
  mc.study$miss.right <- numeric(3)
  
  for (i in 1:M){
    i <- sample(1:n, size = m, replace = TRUE)
    theta.sample <- sk(x[i])
    for (c in 1:3){
      lower <- mc.study[c,1]
      upper <- mc.study[c,2]
      if (theta.sample < lower){
        mc.study[c,3] <- mc.study[c,3] + 1
      }else if (theta.sample > upper){
        mc.study[c,4] <- mc.study[c,4] + 1
      }
    }
  }
  
  mc.study$miss.left <- mc.study$miss.left/M
  mc.study$miss.right <- mc.study$miss.right/M
  mc.study
}

#test on normal and chi-square
mu <- 20
sigma <- 10
n <- 1000 #sample size
x1 <- rnorm(n, mean = mu, sd = sigma)
mc.sk(x1)

df <- 5
x2 <- rchisq(n, df = df)
mc.sk(x2)

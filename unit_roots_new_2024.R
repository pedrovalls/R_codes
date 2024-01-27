# clear workspace
rm(list=ls())
# set local directory to where data set is 
# setwd("C:/Users/Pedro/Dropbox/EcoIII2021/Lecture4_Unit_Root")
setwd("~/Dropbox/EcoIII2021/Lecture4_Unit_Root")

# Load package using a function load_package-----------------------------------------------------------------
load_package<-function(x){
  x<-as.character(match.call()[[2]])
  if (!require(x,character.only=TRUE)){
    install.packages(pkgs=x,repos="http://cran.r-project.org")
    require(x,character.only=TRUE)
  }
}

load_package("tseries")

# sample size
n <- 1000

# repetitions
num_simulations <- 1000

# fix random seed
set.seed(123456)


# Pre-allocate vectors for storing results
rphi1 <- rep(0, num_simulations)
rt1 <- rep(0, num_simulations)
rphi2 <- rep(0, num_simulations)
rt2 <- rep(0, num_simulations)
rphi3 <- rep(0, num_simulations)
rt3 <- rep(0, num_simulations)
y1 <- numeric(num_simulations)
y2 <- numeric(num_simulations)
y3 <- numeric(num_simulations)



# First Process ----------------------------------------------------------

# generate random walk 
rphi1 <- 0
rt1 <- 0 
for(j in 1:num_simulations){
  y1[1] <- rnorm(n)[1]
  for(i in 2:n){
    y1[i] <- y1[i-1] + rnorm(n)[i] 
  }
  y1_1 <- c(NA, y1[1:length(y1)-1]) # compute the lagged series
  eq1 <- lm(y1 ~ y1_1 +0) # OLS estimation no intercept
  rphi1[j] <- summary(eq1)$coefficients[1]  # storing the slope
  rt1[j] <- (summary(eq1)$coefficients[1]-1)/summary(eq1)$coefficients[2]  # storing the t statistic 
}
options(digits=12)

# histogram of the coefficients

par(mfrow=c(1,1))
# plot histogram
hist(rphi1, breaks=30, freq =F, xlab = '', ylab='', main='')
##
# declare function dist the plot the normal distribution with the same mean and standard deviation as the series gwn
##
dist<- function(n) { 
  dnorm(n,mean(rphi1),sd(rphi1))
}

curve(dist, add=T, col='red')

d <- density(rphi1)
lines(d, col='blue')
legend('topleft', legend=c('rphi', 'Normal', 'Kernel'), 
       col=c(1,2,"blue"), pch=15)

summary(rphi1)



# histogram of coefficients t-statistic 


par(mfrow=c(1,1))
# plot histogram
hist(rt1, breaks=30, freq =F, xlab = '', ylab='', main='')
##
# declare function dist the plot the normal distribution with the same mean and standard deviation as the series gwn
##
dist<- function(n) { 
  dnorm(n,mean(rt1),sd(rt1))
}

curve(dist, add=T, col='red')

d <- density(rt1)
lines(d, col='blue')
legend('topleft', legend=c('rt1', 'Normal', 'Kernel'), 
       col=c(1,2,"blue"), pch=15)

summary(rt1)


# Second Process ----------------------------------------------------------

# generate random walk with drift
rphi2 <- 0
rt2 <- 0 
for(j in 1:num_simulations){
  y2[1] <- 0.5++ rnorm(n)[1] 
  for(i in 2:n){
    y2[i] <- 0.5 + y2[i-1] + rnorm(n)[i] 
  }
  y2_1 <- c(NA,y2[1:length(y2)-1]) # compute the lagged series
  eq2 <- lm(y2 ~ y2_1) # OLS estimation 
  rphi2[j] <- summary(eq2)$coefficients[2]  # storing the slope
  rt2[j] <- (summary(eq2)$coefficients[2]-1)/summary(eq2)$coefficients[4]  # storing the t statistic 
}

# histogram of the coefficients

par(mfrow=c(1,1))
# plot histogram
hist(rphi2, breaks=30, freq =F, xlab = '', ylab='', main='')
##
# declare function dist the plot the normal distribution with the same mean and standard deviation as the series gwn
##
dist<- function(n) { 
  dnorm(n,mean(rphi2),sd(rphi2))
}

curve(dist, add=T, col='red')

d <- density(rphi2)
lines(d, col='blue')
legend('topleft', legend=c('rphi2', 'Normal', 'Kernel'), 
       col=c(1,2,"blue"), pch=15)

summary(rphi2)



# histogram of coefficients t-statistic 


par(mfrow=c(1,1))
# plot histogram
hist(rt2, breaks=30, freq =F, xlab = '', ylab='', main='')
##
# declare function dist the plot the normal distribution with the same mean and standard deviation as the series gwn
##
dist<- function(n) { 
  dnorm(n,mean(rt2),sd(rt2))
}

curve(dist, add=T, col='red')

d <- density(rt2)
lines(d, col='blue')
legend('topleft', legend=c('rt2', 'Normal', 'Kernel'), 
       col=c(1,2,"blue"), pch=15)

summary(rt2)



# Third Process -----------------------------------------------------------

# generate random walk with linear trend and drift
trend <- seq(1:n)
rphi3 <- 0
rt3 <- 0 
for(j in 1:num_simulations){
  y3[1] <- 0.1 + 0.5 + rnorm(n)[1]
  for(i in 2:n){
    y3[i] <- 0.1 + 0.5*trend[i-1] + y3[i-1] + rnorm(n)[i] 
  }
  y3_1 <- c(NA,y3[1:length(y3)-1]) # compute the lagged series
  eq3 <- lm(y3 ~ trend + y3_1) # OLS estimation 
  rphi3[j] <- summary(eq3)$coefficients[3]# storing the slope
  rt3[j] <- (summary(eq3)$coefficients[3]-1)/summary(eq3)$coefficients[6]  # storing the t statistic 
}


# histogram of the coefficients

par(mfrow=c(1,1))
# plot histogram
hist(rphi3, breaks=30, freq =F, xlab = '', ylab='', main='')
##
# declare function dist the plot the normal distribution with the same mean and standard deviation as the series gwn
##
dist<- function(n) { 
  dnorm(n,mean(rphi3),sd(rphi3))
}

curve(dist, add=T, col='red')

d <- density(rphi3)
lines(d, col='blue')
legend('topleft', legend=c('rphi3', 'Normal', 'Kernel'), 
       col=c(1,2,"blue"), pch=15)

summary(rphi3)



# histogram of coefficients t-statistic 


par(mfrow=c(1,1))
# plot histogram
hist(rt3, breaks=30, freq =F, xlab = '', ylab='', main='')
##
# declare function dist the plot the normal distribution with the same mean and standard deviation as the series gwn
##
dist<- function(n) { 
  dnorm(n,mean(rt3),sd(rt3))
}

curve(dist, add=T, col='red')

d <- density(rt3)
lines(d, col='blue')
legend('topleft', legend=c('rt3', 'Normal', 'Kernel'), 
       col=c(1,2,"blue"), pch=15)

summary(rt3)



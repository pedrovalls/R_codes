# set local directory to where data set is 
setwd("C:/Users/Pedro/Dropbox/EcoIII2021/Lecture7_var_vec/VEC/R_script")

# Load package using a function load_package-----------------------------------------------------------------
load_package<-function(x){
  x<-as.character(match.call()[[2]])
  if (!require(x,character.only=TRUE)){
    install.packages(pkgs=x,repos="http://cran.r-project.org")
    require(x,character.only=TRUE)
  }
}

load_package('tseries')
load_package('urca')
load_package('dlm')
load_package('openxlsx')
load_package('xts')
load_package('class')
load_package('zoo')
load_package('fBasics')
load_package('qrmtools')
load_package('stats')
load_package('MTS')
load_package('vars')
load_package('graphics')
load_package('MASS')
load_package('urca')


# sample size
n <- 1000

##
# Fix the random seed 
##
set.seed(123456)

##
# generate a multivariade normal distribuition with mean = (0,0), Sigma = diag(0.5,0.5) 
# using mvrnorm
##

Sigma <- matrix(c(0.5,0.0,
                  0.0,0.5),2)
mgwn <- mvrnorm(n, rep(0,2), Sigma, tol = 1e-6, empirical = TRUE)

gwn1=ts(mgwn[,1])
gwn2=ts(mgwn[,2])

##
# check the correlation matrix of gwn1 and gwn2
##
cor(mgwn)
#
#
# simulated bivariate cointegrated system with beta=(1,-1)
# y1(t)=y2(t)+u(t)
# y2(t)=y2(t-1)+gwn1(t)
# u(t)=0.75*u(t-1)+gwm2(t)
##

y1 <- 0
y2 <- 0 
u <- 0
u[1]=gwn2[1]
y2[1]=gwn1[1]
y1[1]=y2[1]+u[1]
  for(i in 2:n){
    u[i] <- 0.75*u[i-1]+gwn2[i]
    y2[i]<- y2[i-1]+gwn1[i]
    y1[i] <- y2[i] + u[i] 
  }

##
#plot the two series in separate plots
##
par(mfrow=c(1,1))
plot(y1,type='l', col='blue', main = 'First Component y1')
plot(y2,type='l', col='red', main = 'Second Component y2')

##
# test unit roots in y1 and y2 using adf,test in tseries
##
y1_adf=adf.test(y1)
y1_adf

y2_adf=adf.test(y2)
y2_adf


##
# generate the EqCM = y1 - y2
##

EqCM=y1-y2


##
# plot EqCM series
##
plot(EqCM,type='l', col='darkblue', main = 'EqCM = y1 - y2')

##
# ACF for EqCm
#
par(mfrow=c(2,1))
acf(EqCM, lag.max=12)
pacf(EqCM, lag.max=12)

##
# test unit roots in EqCM
##
EqCM_adf = adf.test(EqCM)
EqCM_adf


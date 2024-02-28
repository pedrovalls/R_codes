##
# SV Model  QMLE using Kalman Filter
#

# clear workspace
rm(list=ls())
# set local directory to where data set is 
# setwd("C:/Users/Pedro/Dropbox/EcoIII2020/Lecture4_Unit_Root")

# Load package using a function load_package-----------------------------------------------------------------
load_package<-function(x){
  x<-as.character(match.call()[[2]])
  if (!require(x,character.only=TRUE)){
    install.packages(pkgs=x,repos="http://cran.r-project.org")
    require(x,character.only=TRUE)
  }
}

load_package("tseries")
load_package("astsa")
load_package("forecast")
load_package("urca")
load_package("xtable")
load_package("readxl")
load_package("stats")
library(xtable)
library(forecast)
library(tseries) 
library(urca)
library(astsa)
library(stats)
bolsa_ibv <- read_excel("bolsa_ibv.xls")
IBOV_ts=ts(bolsa_ibv$ibvf)
# 
# return of Log IBOV
#
LIBOV=log(IBOV_ts)
LIBOV_1<- stats::lag(LIBOV,-1)
u = ts.intersect(LIBOV,LIBOV_1) 
RLIBOV=u[,1]-u[,2]

#plot the series
par(mfrow=c(1,1))

tsplot(RLIBOV, type ='l',col=4, pch=19, ylab=expression(y[~t]), main="Log Return of IBOV")

##
# in order to avoid log of zero center RLIBOV
#
RLIBOVC=RLIBOV-mean(RLIBOV)

##
# SV using SVfilter
#
SV_RLIBOVC = SV.mle(RLIBOVC,feedback=FALSE)

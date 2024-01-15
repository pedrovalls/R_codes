# clear workspace
rm(list=ls())
# set local directory to where data set is 
# setwd("C:/Users/Pedro/Dropbox/EcoIII2020/Lecture1")
setwd("~/Dropbox/ecoiii2021/lecture3")
# Load package using a function load_package-----------------------------------------------------------------
load_package<-function(x){
  x<-as.character(match.call()[[2]])
  if (!require(x,character.only=TRUE)){
    install.packages(pkgs=x,repos="http://cran.r-project.org")
    require(x,character.only=TRUE)
  }
}

load_package('tseries')
load_package('stats')
load_package('dplyr')
load_package('forecast')
load_package('lmtest')
load_package('TSA')
load_package("readxl")
load_package("data.table")
load_package("xtable")
library(forecast)
library(stats)
library(lmtest)
library(tidyverse)
library(modelr)
library(broom)
library(TSA)
library(data.table)
library(readxl)
library(xtable)
bolsa_cambio <- read_excel("~/Dropbox/EcoIII2021/Lecture3/bolsa_cambio.xls")
View(bolsa_cambio)    

class(bolsa_cambio)

# Assuming dldolar is your time series data
# You should replace 'dldolar' with your actual time series object

dldolar = diff(log(bolsa_cambio$Cambiof),lag=1)
dldolar <-ts(dldolar)

# ACF for AR(2)
par(mfrow=c(2,1))
Acf(dldolar, main = "" ,lag.max=12)
Pacf(dldolar, main = "" ,lag.max=12)


# Fit the ARMA model
# using a Loop and the arima and MLE and with and without mean
#

#
# without mean 
#

start.time_no_mean <- Sys.time()
result_bic_no_mean <- matrix(data=0, nrow=13, ncol=13)
result_bic_mean <- matrix(data=0, nrow=13, ncol=13)
final.bic_no_mean <- Inf
final.bic_mean <- Inf
final.order_no_mean <- c(0,0,0)
final.order_mean <- c(0,0,0)
for (i in 0:12) for (j in 0:12) {
      current.bic_no_mean <- BIC(arima(dldolar, order=c(i, 0, j), include.mean = FALSE, method = "ML"))
      print(i)
      print(j)
      print(current.bic_no_mean)
      result_bic_no_mean[i+1,j+1] = current.bic_no_mean
      print(result_bic_no_mean[i+1,j+1])
  if (current.bic_no_mean < final.bic_no_mean) {
         final.bic_no_mean <- current.bic_no_mean
         final.order_no_mean <- c(i, 0, j)
         final.arma_no_mean <- arima(dldolar, order=final.order_no_mean,include.mean = FALSE, method = "ML")
      }
  }


print(result_bic_no_mean)
print(final.order_no_mean)
print(final.arma_no_mean)

result_bic_no_mean_latex <- xtable(result_bic_no_mean, caption=NULL, label = NULL, align = "c", digits = 4, display = NULL, auto = FALSE)
print(result_bic_no_mean_latex)

end.time_no_mean <- Sys.time()
time.taken_no_mean <- end.time_no_mean - start.time_no_mean
time.taken_no_mean

#
# with mean 
#
start.time <- Sys.time()


result_bic_mean <- matrix(data=0, nrow=13, ncol=13)
final.bic_mean <- Inf
final.order_mean <- c(0,0,0)
for (i in 0:12) for (j in 0:12) {
  current.bic_mean <- BIC(arima(dldolar, order=c(i, 0, j), include.mean = TRUE, method = "ML"))
  print(i)
  print(j)
  print(current.bic_mean)
  result_bic_mean[i+1,j+1] = current.bic_mean
  print(result_bic_mean[i+1,j+1])
  if (current.bic_mean < final.bic_mean) {
    final.bic_mean <- current.bic_mean
    final.order_mean <- c(i, 0, j)
    final.arma_mean <- arima(dldolar, order=final.order_mean,include.mean = TRUE, method = "ML")
  }
}
print(result_bic_mean)
print(final.order_mean)
print(final.arma_mean)

result_bic_mean_latex <- xtable(result_bic_mean, caption=NULL, label = NULL, align = "c", digits = 4, display = NULL, auto = FALSE)
print(result_bic_mean_latex)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken
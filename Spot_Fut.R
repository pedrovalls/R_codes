#' ---
#' title: "Spot and Fut VEC modelling"
#' author: 
#'  - "Pedro Valls"
#'  
#' date: "September 23th, 2020"
#' ---








#' set local directory to where data set is 
setwd("C:/Users/Pedro/Dropbox/ecoiii2021/Lecture7_var_vec/VEC/R_script")

#' Load package using a function load_package-----------------------------------------------------------------
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
load_package('readxl')
load_package("tidyverse")
load_package("tsDyn")
load_package("ggplot2")

library(ggplot2)
library(tseries)
library(urca)
library(dlm)
library(openxlsx)
library(xts)
library(class)
library(zoo)
library(fBasics)
library(qrmtools)
library(stats)
library(MTS)

library(vars)
library(graphics)
library(readxl)
library(tidyverse)
library(tsDyn)
#'
#'
#' ## Read S&P500 spot and fut already transformed in logs
#'

sp500 <- read_excel("C:/Users/Pedro/Dropbox/ecoiii2020/Lecture7_var_vec/VEC/sp500.xls")

#'
#' transform in time series 
#'
 sp500_spot=ts(log(sp500$SPOT))
 sp500_fut=ts(log(sp500$FUT))
 
#'
#' plot the two series
#'
 par(mfrow=c(1,2))
 plot(sp500_spot,type='l', col='blue', main = 'Spot of log(SP500)')
 plot(sp500_fut,type='l', col='red', main = 'Fut of log(SP500)')
 


 
 
 
 # plotting in the same plot
 par(mfrow=c(1,1))
 
 # Plot the two series in the same plot
 plot(sp500_spot, type='l', col='blue', main='Log(SP500_SPOT) and Log(SP500_Future)', ylab='Value', xlab='Time')
 lines(sp500_fut, col='red')
 
 # Add a legend to distinguish the series
 legend("topright", legend=c("Spot", "Future"), col=c("blue", "red"), lty=1)
 
 
 
 par(mfrow=c(1,1))
 plot.default(sp500_spot, sp500_fut, col='darkblue')
 regressao <- lm(sp500_spot ~ sp500_fut)
  summary(regressao)
 abline(coef(regressao), col='red')
 legend('topleft', legend=c('sp500_spot=0.0613024+0.9661331*sp500_fut'), col=c('red'), pch=15)
 title('Scatter Plot - least squared fitted') 

 
#'
#' ## Define the Equilibrium correction by $sp500spot_{t}-\alpha - \beta*sp500fut_{t}$
#'
 EqCM=sp500_spot-regressao$coefficients[1]-regressao$coefficients[2]*sp500_fut
 
#'
#' plot of EqCM
#'
 plot(EqCM,type='l', col='steelblue', main = 'Equilibrium correction using OLS')
 
 
#'
#' ACF for EqCm
#'
 par(mfrow=c(2,1))
 acf(EqCM, lag.max=12)
 pacf(EqCM, lag.max=12)
 
#'
#' test unit roots in EqCM
#'
 EqCM_adf = adf.test(EqCM)
 EqCM_adf
 
 
#'
#' ## Generate EqCM with $\beta = (1,-1)$ the procedure is not necessarily optimum
#'
 EqCM11 = sp500_spot-sp500_fut
 
#'
#' plot of EqCM
#'
 par(mfrow=c(1,1))
 plot(EqCM11,type='l', col='steelblue', main = 'Equilibrium correction using beta = (1,-1)')
 
 
#'
#' ACF for EqCm
#'
 par(mfrow=c(2,1))
 acf(EqCM11, lag.max=12)
 pacf(EqCM11, lag.max=12)
 
#'
#' test unit roots in EqCM11
#'
 EqCM11_adf = adf.test(EqCM11)
 EqCM11_adf
 
#'
#' ## EqCM as the residual for the static regression of sp500spot into constant and sp500fut this is the two-step procedure of Engle and Granger(1987)
#'
 EqCM_regressao <- lm(sp500_spot ~ sp500_fut)
 summary(EqCM_regressao)
 EqCM_res=ts(EqCM_regressao$residuals)
 
#'
#' plot of EqCM_res
#'
 par(mfrow=c(1,1))
 plot(EqCM_res,type='l', col='steelblue', main = 'Equilibrium correction using beta from OLS')
 
 
#'
#' ACF and PACFfor EqCm_res
#'
 par(mfrow=c(2,1))
 acf(EqCM_res, lag.max=12)
 pacf(EqCM_res, lag.max=12)
 
 
#'
#' test unit roots in EqCM_res
#'
 EqCM_res_adf = adf.test(EqCM_res)
 EqCM_res_adf
 

 
#'
#' lag one period transforme sp500_spot and sp500_fut 
#'
sp500_spot_1=ts(c(NA,sp500_spot[1:length(sp500_spot)-1]))
sp500_fut_1=ts(c(NA,sp500_fut[1:length(sp500_fut)-1]))


#'
#' transforme sp500_spot and sp500_fut in first difference
#'
dsp500_spot = sp500_spot-sp500_spot_1
dsp500_fut = sp500_fut-sp500_fut_1

#'
#' lag one period  dsp500_spot and dsp500_fut 
#'
dsp500_spot_1=ts(c(NA,dsp500_spot[1:length(dsp500_spot)-1]))
dsp500_fut_1=ts(c(NA,dsp500_fut[1:length(dsp500_fut)-1]))

#'
#' lag one period EqCM_res
#'
EqCM_res_1=ts(c(NA,EqCM_res[1:length(EqCM_res)-1]))

#'
#' ## Equation for $dsp500spot_{t} = \alpha_{1} + \beta_{11}*EqCMres_{t-1} + \beta_{12}* dsp500spot_{t-1} + \beta_{13}*dsp500fut_{t-1}$
#'
eq_spot <- lm(dsp500_spot ~ EqCM_res_1 + dsp500_spot_1 +  dsp500_fut_1)
summary(eq_spot)


#'
#' ## Equation for $dsp500fut_{t} = \alpha_{2} + \beta_{21}*EqCMres_{t-1} + \beta_{22}* dsp500spot_{t-1} + \beta_{23}*dsp500fut{t-1}$
#'
eq_fut <- lm(dsp500_fut ~ EqCM_res_1 + dsp500_spot_1 +  dsp500_fut_1)
summary(eq_fut)



#'
#' ## Define the multivariate time series $(spot_{t}, fut_{t})$
#'
sp500_spot_fut=ts(sp500[,3:2])



#'
#' ## VAR lag length selection using the package VARS
#'

VARselect(sp500_spot_fut,lag.max=20, type=c("const")) 

#'
#' best model is VAR(20) using AIC, VAR(6) using BIC and VAR(8) using HQ
#'

#'
#' we are going to use VAR(20)
#'

VAR_spot_fut20 <-VAR(sp500_spot_fut , p=20, type = "const")

VAR_spot_fut20

#'
#' summary of equation spot
#'
summary(VAR_spot_fut20, equation = "SPOT")

#'
#' plot fit and residuals and ACF and PACF for residuals for SPOT
#'#'
plot(VAR_spot_fut20, names = "SPOT")

#'
#' summary of equation fut
#'
summary(VAR_spot_fut20, equation = "FUT")
#'
#' plot fit and residuals and ACF and PACF for residuals for FUT
#'
plot(VAR_spot_fut20, names = "FUT")


#'
#'  Portmanteau test for system 
#'
ser20 <- serial.test(VAR_spot_fut20, lags.pt = 20, type = "PT.asymptotic")
ser20

#'
#' JB Test for system
#'
JB20 <- normality.test(VAR_spot_fut20)
JB20$jb.mul

#'
#' ARCH test for system
#'
arch20 <- arch.test(VAR_spot_fut20, lags.multi = 20)
arch20$arch.mul


#'
#' diagnostic for spot equation
#'
plot(arch20,name= "SPOT")
plot(stability(VAR_spot_fut20),nc=2)


plot(arch20,name= "FUT")



#'
#' Causality test
#' 

#'
#'Granger causality H0: FUT do not Granger-cause SPOT
#'

causality(VAR_spot_fut20, cause=c("FUT"))

#'
#'Granger causality H0: SPOT  do not Granger-cause FUT
#'

causality(VAR_spot_fut20, cause=c("SPOT"))

#'
#' ## Johansen procedure for SPOT and FUT
#'

vec_spot_fut_eigen <- (ca.jo(sp500_spot_fut, type=c("eigen"), ecdet = c("const"), K=19, spec = c("transitory")))
summary(vec_spot_fut_eigen)

vec_spot_fut_trace <- (ca.jo(sp500_spot_fut, type=c("trace"), ecdet = c("const"), K=19, spec = c("transitory")))
summary(vec_spot_fut_trace)

vec_z<- ca.jo(sp500_spot_fut, type=c("eigen"), ecdet = c("const"), K=19, spec = c("longrun"))
summary(vec_z)


#' ## Specifying a VEqCM
#' 
#' We will specify a **restricted** VEqCM with constant inside the EqCM
#' but no deterministic component in the short run dynamic. 
#' Since our VAR model is of order 20, the VEC model
#' will be of order 19, and the long term relation.
#' 
#' This long-term relation has already been estimated by our VAR(20) 
#' regression, thus, the VEqCM specification:
#' 
#' $$\Delta y_t = \alpha\beta^\prime y_{t-1}+\varepsilon_t \quad\text{VEqCM}$$
#' 
#' is equivalent to the VAR(1) specification
#' 
#' $$y_t = \Pi y_{t-1}+\varepsilon_t \quad\text{VAR}$$
#' 
#'  where $\Pi = \mathbf{I}+\alpha\beta^\prime$. Therefore, from our VAR(1) model
#'  already estimated we have all information we need to have the VEqCM 
#'  representation. 
#'  
#'  The estimated matrix $\Pi$ is given by:
#'  
Bcoef(VAR_spot_fut20)


#'  We can use the function VECM from tsDyn package to estimate the VECM(19) with 
#'  one cointegrating vectors with constant in the EqCM (**LRinclude=c("const")**) 
#'  but not in the VAR (**include=c("none")**)
#'  
vec <- VECM(sp500_spot_fut, lag = 19, r=1, estim = "ML", include = c("none"), LRinclude = c("const"))
summary(vec)



#' where  EqCM = $\beta^{'}y_{t-1}$


#'
#' ## Test a restriction in $\alpha$ using alrtest $\alpha_{1}=0$ using the package urca weak exogeneity for spot in fut conditional model
#' 





DA_spot = matrix (c(0,1),c(2,1))
DA_spot
teste_alpha_spot <- alrtest(vec_spot_fut_eigen,A=DA_spot,r=1)
summary(teste_alpha_spot)



#'
#' ## Test a restriction in $\alpha$ using alrtest $\alpha_{2}=0$ using the package urca weak exogeneity for fut in spot conditional model 
#' 



DA_fut = matrix (c(1,0),c(2,1))
DA_fut
teste_alpha_fut <- alrtest(vec_spot_fut_eigen,A=DA_fut,r=1)
summary(teste_alpha_fut)


#'
#' ## Test a restriction in $\beta$ using blrtest $\beta_{1}=1$ and $\beta_{2}=-1$ 
#' using the package urca for an one-to-one in the long-run 
#' relationship between spot and fut 
#'  

HD_prop <- matrix(c(1,-1,0),c(3,1))
HD_prop
test_beta_prop <- blrtest(vec_spot_fut_eigen,H=HD_prop,r=1)
summary(test_beta_prop)

#'
#´## Test both restrictions $\beta_{1}=1$ and $\beta_{2}=-1$  and $\alpha_{2}=0$ 
#'
test_alpha_beta_fut <- ablrtest(vec_spot_fut_eigen,H=HD_prop,A=DA_fut,r=1)
summary(test_alpha_beta_fut)



#'
#´## Test both restrictions $\beta_{1}=1$ and $\beta_{2}=-1$  and $\alpha_{1}= 0$ 
#'
test_alpha_beta_spot <- ablrtest(vec_spot_fut_eigen,H=HD_prop,A=DA_spot,r=1)
summary(test_alpha_beta_spot)


#'
#' ## IRF
#' 


 vec2var_vec_z <- vec2var(vec_z,r=1)


irf_vec_z_ort <- irf(vec2var_vec_z,n.ahead=10,ortho=TRUE)

plot(irf_vec_z_ort)

irf_vec_z_cumul <- irf(vec2var_vec_z,n.ahead=10,ortho=FALSE, cumulative = TRUE)

plot(irf_vec_z_cumul)

irf_vec_z <- irf(vec2var_vec_z,n.ahead=10,ortho=FALSE, cumulative = FALSE)

plot(irf_vec_z)


#'
#' ## Variance Decomposition
#' 

fevd(vec2var_vec_z, n.ahead=10)

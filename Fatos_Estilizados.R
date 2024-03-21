##
# Stylized Facts in Finance
##

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
load_package("lmtest")
load_package("ggplot2")
load_package("sandwich")
load_package("zoo")
load_package("car")
load_package("gets")

library(xtable)
library(forecast)
library(tseries) 
library(urca)
library(astsa)
library(stats)
library(lmtest) # for linear regression diagnostics
library(ggplot2) # for plotting
library(sandwich)
library(zoo)
library(car)
library(gets)


#setwd("C:/Users/Pedro/Dropbox/EcoIII2021/Lecture6_adl/R_code")
DADOS_BOLSA <- read_excel("C:/Users/Pedro/Dropbox/EcoIII2021/Lecture_Volatilidade_Univariada/Vol_R/DADOS_BR.xlsx")
##
# trasnform IBOVC in time series
##
IBOVC_ts = ts(DADOS_BOLSA$IBOVC)
##
# num is the length of the time series
##
num = length(IBOVC_ts)

##
# transition equation is local level x_{t} = x_{t-1} +omega_{t}
# observational equation is y_{t} = A_{t} x_{t} + v_{t} 
# where A_{t} is 1 when y_{t} is observed and A_{t}= 0  when y_{t} not observed
##
y = as.numeric(IBOVC_ts)
num = length(IBOVC_ts)
##
# define the A_{t} vector
##
A = array(0,dim=c(1,1,num))
##
# define  A_{t} is 1 when y_{t} is observed and A_{t}= 0  when y_{t} not observed
##
for(k in 1:num) if(!is.na(y[k])) A[,,k] = diag(1,1)
##
# Initial values
##
mu0    = matrix(0,1,1)
Sigma0 = diag(c(.1) ,1)
Phi    = diag(1, 1)
Q     = diag(c(1), 1) 
R     = diag(c(1),1) 
##
# Run EM
#
(em = EM(y, A, mu0, Sigma0, Phi, Q, R))

# Run smoother at the estimates (using the new Ksmooth script)
sQ = em$Q^.5
sR = sqrt(em$R)
ks  = Ksmooth(y, A, em$mu0, em$Sigma0, em$Phi, sQ, sR)
##
# Pull out the values
##
y1s = ks$Xs[1,,]
p1  = 2*sqrt(ks$Ps[1,1,])

##
# Irregular
##

Irregular = rep(0,num)

for(k in 1:num) if (!is.na(y[k])) {Irregular[k] = y[k] - y1s[k]}

##
# series without missing
##
y_no_missing= rep(0,num)
for(k in 1:num) y_no_missing[k] = y1s[k]+Irregular[k]


# plots
par(mfrow=c(1,1))
tsplot(y_no_missing, type='p', pch=19, ylim=c(8300,131000), col=6, lwd=2, cex=1)
lines(y_no_missing)
#xx = c(time(y_no_missing), rev(time(y_no_missing)))
xx = c(DADOS_BOLSA$Date, DADOS_BOLSA$Date)
yy = c(y_no_missing-p1, rev(y_no_missing+p1))
polygon(xx, yy, border=8, col=astsa.col(8, alpha = .1))

write.csv(y_no_missing, file="C:/Users/Pedro/Dropbox/EcoIII2021/Lecture_Volatilidade_Univariada/Vol_R/IBOVC_SB.csv" )

##
# Plot IBOV without Missing

tsplot(DADOS_BOLSA$Date,y_no_missing, type='l', pch=19, ylim=c(8300,131000), col=2, lwd=2, cex=1, ylab="IBOVC_SB",
       main = "IBOV close wihout missing") 

##
# Define the return of IBOVC
##
Ibov = y_no_missing

Ibov_ts = ts(Ibov)

u = ts.intersect(Ibov_ts,
                 stats::lag(Ibov_ts,-1))
RLIbov = (log(u[,1]) - log(u[,2]))*100

##
# plot the log return of Ibvoc
##
tsplot(DADOS_BOLSA$Date[2:num],RLIbov, type='l', pch=19, 
       #ylim=c(8300,131000),
       xlim=c(min(DADOS_BOLSA$Date), max(DADOS_BOLSA$Date)),
       col=2, lwd=2, cex=1, ylab="RLIbov",
       main = "DLIbovc") 
##
# squared returns
##
RLIbov_sq = RLIbov^2

##
# Plot Squared Returns and FAC and FACP of Squared Returns
##

par(mfrow=c(1,1))

tsplot(DADOS_BOLSA$Date[2:num],RLIbov_sq, type='l', pch=19, 
       #ylim=c(8300,131000),
       xlim=c(min(DADOS_BOLSA$Date), max(DADOS_BOLSA$Date)),
       col=2, lwd=2, cex=1, ylab="RLIbov_sq",
       main = "DLIbovc^2") 

par(mfrow=c(2,1))
Acf(RLIbov_sq, lag.max = 24)
Pacf(RLIbov_sq, lag.max = 24)

##
# tsdiag(series, gof = n)
##
LB_RLIbov_sq = Box.test(RLIbov_sq, lag = 50, type = "Ljung-Box", fitdf = 0)
LB_RLIbov_sq
##



u_engle = ts.intersect(RLIbov_sq,
                 stats::lag(RLIbov_sq,-1),
                 stats::lag(RLIbov_sq,-2),
                 stats::lag(RLIbov_sq,-3),
                 stats::lag(RLIbov_sq,-4),
                 stats::lag(RLIbov_sq,-5),
                 stats::lag(RLIbov_sq,-6),
                 stats::lag(RLIbov_sq,-7),
                 stats::lag(RLIbov_sq,-8),
                 stats::lag(RLIbov_sq,-9),
                 stats::lag(RLIbov_sq,-10),
                 stats::lag(RLIbov_sq,-11),
                 stats::lag(RLIbov_sq,-12))
Test_engle <- lm(u_engle[,1] ~ u_engle[,2]+
                   u_engle[,3]+
                   u_engle[,4]+
                   u_engle[,5]+
                   u_engle[,6]+
                   u_engle[,7]+
                   u_engle[,8]+
                   u_engle[,9]+
                   u_engle[,10]+
                   u_engle[,11]+
                   u_engle[,12]+
                   u_engle[,13])
summary(Test_engle)

num = length(RLIbov_sq)
num
AIC_Test_engle= AIC(Test_engle)
AIC_Test_engle
AIC_Test_engle_T=AIC_Test_engle/num
AIC_Test_engle_T
BIC_Test_engle = BIC(Test_engle)

BIC_Test_engle
BIC_Test_engle_T=BIC_Test_engle/num
BIC_Test_engle_T
npar_Test_engle = length(Test_engle$coefficients)
loglik_Test_engle = -(1/2)*(AIC_Test_engle - 2*npar_Test_engle)
loglik_Test_engle


##
# Sampling Variance
##


moving_sd_22 = runsd(RLIbov, 
                     22, 
                     center = runmean(RLIbov,22))
moving_sd_44 = runsd(RLIbov, 
                     44, 
                     center = runmean(RLIbov,44))
moving_sd_66 = runsd(RLIbov, 
                     66, 
                     center = runmean(RLIbov,66))
moving_sd_126 = runsd(RLIbov, 
                     126, 
                     center = runmean(RLIbov,126))
moving_sd_252 = runsd(RLIbov, 
                     252, 
                     center = runmean(RLIbov,252))
par(mfrow=c(1,1))

yy = cbind(moving_sd_22,moving_sd_44,moving_sd_66,moving_sd_126,moving_sd_252)



matplot(DADOS_BOLSA$Date[1:num], yy, type='l', lty = 1,
        xlim=c(min(DADOS_BOLSA$Date), max(DADOS_BOLSA$Date)),
        col = 2:6, lwd = 2, xlab = "Date", ylab = "Moving Standard Deviation",
        main = "Moving Standard Deviations using Different Window Sizes for RLIbov")
legend("topleft", legend = c("22 days", "44 days", "66 days", "126 days", "252 days"),
       col = 2:6, lty = 1, lwd = 2)

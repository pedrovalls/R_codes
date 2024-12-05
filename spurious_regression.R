##
# Regress\~{a}o Esp\'{u}ria
## 

# set local directory to where data set is 
# setwd("~/Users/pedrovallspereira/Dropbox/ecoiii2020/Lecture7_var_vec/VEC/R_script")
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




# Fix the random seed
set.seed(123456)
# matrix A of coeff of VAR(1)
p1 = matrix(c(1.0,0.0,
              0.0,1.0),2)
# matrix SIGMA var-cov of disturbances
sig=matrix(c(1.0, 0.0,
             0.0,1.0),2)
# simulate a VAR(1)
m1=VARMAsim(200,arlags=c(1),malags=NULL,phi=p1,theta=NULL,sigma=sig)
# rename series to yt
yt=m1$series

# transform yt in time series

yt_ts=ts(yt)

# rename the first component of yt to y1 and the second to y2
y1=yt_ts[,1]
y2=yt_ts[,2]
#plot the two series
par(mfrow=c(1,2))
plot(y1,type='l', col='blue', main = 'First RW Component')
plot(y2,type='l', col='red', main = 'Second RW Component')
##
# scatter plot
##
par(mfrow=c(1,1))
plot(y1,y2,pch='8', col='black', main = 'Scatter Plot of two RW')

par(mfrow=c(1,1))

regressao <- lm(y2 ~ y1)
summary(regressao)
abline(coef(regressao), col='red')

legend('topleft', legend=c('y2=18.93392 - 0.29451*y1'), col=c('red'), pch=15)



##
# correlation between the two series
##
cor(y1,y2)


##
# define the lagged series
##

y1_1=c(NA,y1[1:length(y1)-1])

##
# defined the difference series and called dy1
##

dy1=y1-y1_1


##
# define the lagged series
##

y2_1=c(NA,y2[1:length(y2)-1])

##
# defined the difference series and called dy2
##

dy2=y2-y2_1



##
# plot the first difference of the two series - dy1 and dy2
##


par(mfrow=c(1,2))
plot(dy1,type='l', col='blue', main = 'First Difference of y1')
plot(dy2,type='l', col='red', main = 'First Difference of y2')

##
# scatter plot
##

par(mfrow=c(1,1))
plot(dy1,dy2,pch='8', col='black', main = 'Scatter Plot of dy1 and dy2')


regressaod <- lm(dy2 ~ dy1)
summary(regressaod)
abline(coef(regressaod), col='red')

legend('topleft', legend=c('dy2=-0.07150+0.03865*dy1'), col=c('red'), pch=15)



##
# regression in levels
##
eq1_level<- lm(y2 ~ y1)

##
# summary of the regression in levels y2 = a +b*y1
##

summary.lm(eq1_level)



##
# regression in fisrt difference 
##
eq1_diff<- lm(dy2 ~ dy1)

##
# summary of the regression in difference dy2 = a +b*dy1
##

summary.lm(eq1_diff)



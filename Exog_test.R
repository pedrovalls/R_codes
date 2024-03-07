##
# Exogeneity Test Durbin-Wu_Hausmann
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
library(xtable)
library(forecast)
library(tseries) 
library(urca)
library(astsa)
library(stats)

setwd("C:/Users/Pedro/Dropbox/EcoIII2021/Lecture6_adl/R_code")
datagive1 <- read_excel("datagive1.xlsx")
inc = datagive1$inc[2:159]
cons = datagive1$cons[2:159]
cons1=datagive1$cons1[2:159]
eq1 <- lm(inc ~ cons+cons1)
res_eq1 = eq1$residuals
eq2 <- lm(cons ~ inc+res_eq1)
summary(eq2)

######################################################
############### Econometrics III #################
################### Valls Windows ####################
######################################################
# clean workplace
rm(list=ls())

# set local directory to where data set is 
setwd("C:/Users/Pedro/Dropbox/EcoIII2020/Lecture1")

# Load Packages fucntion -----------------------------------------------------------------
load_package<-function(x){
  x<-as.character(match.call()[[2]])
  if (!require(x,character.only=TRUE)){
    install.packages(pkgs=x,repos="http://cran.r-project.org")
    require(x,character.only=TRUE)
  }
}

load_package('readxl')
load_package('zoo')
###

ukhp<-read_xls('UKHP.xls')
# Change name to hp
hp<-zoo(ukhp$`Average House Price`, order.by=ukhp$Month)

# CI [m(hp)+1.96*sd(hp),m(hp)-1.96*sd(hp)] for the series in levels
mean_hp<-mean(hp[1:77])
stdev_hp<-sd(hp[1:77])
plot(hp)
abline(h=mean_hp+1.96*stdev_hp, col='red', lty=2)
abline(h=mean_hp-1.96*stdev_hp, col='red', lty=2)
# Define series of monthly percentage changes
dhp<-100*(hp-lag(hp,-1))/lag(hp,-1)
plot(dhp)
# CI [m(dhp)+1.96*sd(dhp),m(dhp)-1.96*sd(dhp)] for the series in first difference
mean_dhp<-mean(dhp[1:76])
stdev_dhp<-sd(dhp[1:76])
plot(dhp)
abline(h=mean_dhp+1.96*stdev_dhp, col='red', lty=2)
abline(h=mean_dhp-1.96*stdev_dhp, col='red', lty=2)

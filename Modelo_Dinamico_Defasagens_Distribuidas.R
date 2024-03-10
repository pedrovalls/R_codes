##
# Dynamic Models Typology - Distributed Lags 
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

setwd("C:/Users/Pedro/Dropbox/EcoIII2021/Lecture6_adl/R_code")

datagive1 <- read_excel("datagive1.xlsx")
class(datagive1)





##
# Distributed Lags 
##
# Perform linear regression using data from2 to 159 
model_DL<- lm(log(cons) ~ log(inc)+log(inc1), data = datagive1[2:159,])

# Summary of the model to display regression results
summary(model_DL)

# Using White correction
coeftest(model_DL, vcov = vcovHC(model_DL, type = "HC0"))
#using Newey-West Correction
coeftest(model_DL, vcov = vcovHAC(model_DL))

num = length(datagive1$cons)-1
num
AIC_DL= AIC(model_DL)
AIC_DL
AIC_DL_T=AIC_DL/num
AIC_DL_T
BIC_DL = BIC(model_DL)

BIC_DL
BIC_DL_T=BIC_DL/num
BIC_DL_T
npar_DL = length(model_DL$coefficients)
loglik_DL = -(1/2)*(AIC_DL - 2*npar_DL)
loglik_DL

# Create a data frame of predicted values and standard errors
predictions_DL <- data.frame(
  lcons_DL = log(datagive1$cons[2:159]),
  Predicted_DL = predict(model_DL, interval = "prediction")[,1],
  Lower_DL = predict(model_DL, interval = "prediction")[,2],
  Upper_DL = predict(model_DL, interval = "prediction")[,3]
)

# Save  residuals 
Residuals_DL <- residuals(model_DL)

# Plot observed vs predicted values with confidence intervals



ggplot(predictions_DL, aes(x = datagive1$Date[2:159])) +
  geom_line(aes(y = Predicted_DL, color = "Predicted"), size = 1) +  # Define color within aes for legend
  geom_ribbon(aes(ymin = Lower_DL, ymax = Upper_DL, fill = "Confidence Interval"), alpha = 0.2) +  # Define fill within aes for legend
  geom_point(aes(y = lcons_DL, color = "Observed"), size = 2) +  # Define color within aes for legend
  scale_color_manual(values = c("Predicted" = "blue", "Observed" = "red")) +  # Manual color values
  scale_fill_manual(values = c("Confidence Interval" = "grey")) +  # Manual fill values
  labs(title = "Observed vs Predicted", x = "Date", y = "Value", 
       color = "Legend", fill = "Legend") +  # Customize legend title
  theme_minimal()


# Plot residuals with 2*sigma^{2}
sigma_value_DL = sigma(model_DL)
sigma_value_DL

ggplot(predictions_DL, aes(x = datagive1$Date[2:159], y = Residuals_DL)) +
  geom_point(aes(color = "Residuals")) +  # Assign color within aes for legend
  geom_line(aes(color = "Residuals")) +  # Use same color as points to keep them in a single legend item
  geom_hline(aes(yintercept = 2 * sigma_value_DL, linetype = "2*sigma"), color = "blue") +  # Map linetype to aes
  geom_hline(aes(yintercept = -2 * sigma_value_DL, linetype = "-2*sigma"), color = "gray") +  # Map linetype to aes
  scale_color_manual(values = c("Residuals" = "red")) +  # Manual color values for points and lines
  scale_linetype_manual(values = c("2*sigma" = "dashed", "-2*sigma" = "dashed")) +  # Manual linetype values
  labs(title = "Residuals", x = "Date", y = "residuals", color = "Legend", linetype = "Legend") +  # Custom legend titles
  theme_minimal()





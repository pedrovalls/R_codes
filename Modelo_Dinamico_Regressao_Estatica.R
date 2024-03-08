##
# Dynamic Models Typology
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
# Static Regression
##
# Perform linear regression
model <- lm(cons ~ inc, data = datagive1)

# Summary of the model to display regression results
summary(model)

coeftest(model, vcov = vcovHC(model, type = "HC0"))

AIC(model)
BIC(model)

# Create a data frame of predicted values and standard errors
predictions <- data.frame(
  cons = datagive1$cons,
  Predicted = predict(model, interval = "prediction")[,1],
  Lower = predict(model, interval = "prediction")[,2],
  Upper = predict(model, interval = "prediction")[,3]
)

# Save  residuals 
Residuals <- residuals(model)

# Plot observed vs predicted values with confidence intervals

ggplot(predictions, aes(x = datagive1$Date)) +
  geom_line(aes(y = Predicted), color = "blue") +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.2) +
  geom_point(aes(y = cons), color = "red") +
  labs(title = "Observed vs Predicted", x = "Date", y = "Value") +
  theme_minimal()

# Plot residuals with 2*sigma^{2}

sigma_value <- sigma(model)  # Calculate or specify the sigma of the model

ggplot(predictions, aes(x = datagive1$Date, y = Residuals)) +
  geom_point() +
  geom_line() +
  geom_point(aes(y = Residuals), color = "red") +
  geom_hline(yintercept = 2 * sigma_value, linetype = "dashed", color = "blue") +  # Add line at 2*sigma
  geom_hline(yintercept = -2 * sigma_value, linetype = "dashed", color = "blue") + # Add line at -2*sigma
  labs(title = "Residuals", x = "Date", y = "residuals") +
  theme_minimal()

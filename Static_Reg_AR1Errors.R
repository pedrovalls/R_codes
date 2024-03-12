##
# Dynamic Models Typology - Static Regrssion with AR(1) error
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
# Static Regression with AR(1) Errors 
##
# Perform linear regression



lcons = log(datagive1$cons)
lcons = lcons[2:159]
lcons1=log(datagive1$cons1)
lcons1 = lcons1[2:159]
linc = log(datagive1$inc)
linc = linc[2:159]

#######

# Assuming 'data' is your dataframe containing the variables 'cons' and 'inc'
# Set the data sample - this is typically implicit in R, as operations are applied to the whole dataset by default

# Estimate the autoregressive model with a constant and an independent variable
model_statRegAR <- Arima(lcons, order=c(1,0,0), xreg=linc, include.mean=TRUE)
summary(model_statRegAR)
AIC_statRegAR = model_statRegAR$aic
AIC_statRegAR

BIC_statRegAR = model_statRegAR$bic
BIC_statRegAR

loglik_statRegAR = model_statRegAR$loglik
loglik_statRegAR

sigma2_statRegAR = model_statRegAR$sigma2
sigma2_statRegAR


# Extract the fitted values
lconshat <- fitted(model_statRegAR)

# Calculate the standard errors of predictions
se <- sqrt(summary(model_statRegAR)$sigma2)

# Generate the confidence intervals
lconsup <- lconshat + 2*se
lconslow <- lconshat - 2*se

# Plot the actual vs fitted values with confidence intervals
# Prepare a data frame for plotting
plot_data <- data.frame(
  Time = datagive1$Date[2:159],  # Adjust 'Data' to your actual time variable name
  lcons,
  lconshat,
  lconsup,
  lconslow
)

# Plot the actual vs fitted values with confidence intervals
ggplot(plot_data, aes(x = Time)) +
  geom_line(aes(y = lcons, color = "Actual"), size = 1) +
  geom_line(aes(y = lconshat, color = "Fitted"), size = 1) +
  geom_ribbon(aes(ymin = lconslow, ymax = lconsup, fill = "Confidence Interval"), alpha = 0.5) +
  scale_color_manual(values = c("Actual" = "blue", "Fitted" = "red")) +
  scale_fill_manual(values = c("Confidence Interval" = "grey80")) +
  labs(title = "Actual vs Fitted Log(Cons) with Confidence Intervals",
       y = "Log(Cons)",
       color = "Legend", # This labels the color legend
       fill = "") + # This ensures the fill legend does not have a title, adjust as needed
  theme_minimal()


# Plot residuals
lconsres <- residuals(model_statRegAR)
# Generate the confidence intervals for residuals
lconsresup <- lconsres + 2*se
lconsreslow <- lconsres - 2*se

ggplot(plot_data, aes(x = Time)) +
  geom_line(aes(y = lconsres, color = "Residuals"), size = 1) +
  geom_ribbon(aes(ymin = lconsreslow, ymax = lconsresup, fill = "Confidence Interval"), alpha = 0.5) +
  scale_color_manual(values = c("Residuals" = "blue")) +
  scale_fill_manual(values = c("Confidence Interval" = "grey80")) +
  labs(title = "Residuals with Confidence Intervals",
       y = "Residuals",
       color = "Legend", # This labels the color legend
       fill = "") + # This ensures the fill legend does not have a title, adjust as needed
  theme_minimal()
########








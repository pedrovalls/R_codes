##
# Dynamic Models Typology  - Partial Adjustment
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
model_adjpar<- lm(log(cons) ~ log(inc)+log(cons1), data = datagive1[2:159,])

# Summary of the model to display regression results
summary(model_adjpar)

# Using White correction
coeftest(model_adjpar, vcov = vcovHC(model_adjpar, type = "HC0"))
#using Newey-West Correction
coeftest(model_adjpar, vcov = vcovHAC(model_adjpar))

num = length(datagive1$cons)-1
num
AIC_adjpar= AIC(model_adjpar)
AIC_adjpar
AIC_adjpar_T=AIC_adjpar/num
AIC_adjpar_T
BIC_adjpar = BIC(model_adjpar)

BIC_adjpar
BIC_adjpar_T=BIC_adjpar/num
BIC_adjpar_T
npar_adjpar = length(model_adjpar$coefficients)
loglik_adjpar = -(1/2)*(AIC_adjpar - 2*npar_adjpar)
loglik_adjpar

# Create a data frame of predicted values and standard errors
predictions_adjpar <- data.frame(
  lcons_adjpar = log(datagive1$cons[2:159]),
  Predicted_adjpar = predict(model_adjpar, interval = "prediction")[,1],
  Lower_adjpar = predict(model_adjpar, interval = "prediction")[,2],
  Upper_adjpar = predict(model_adjpar, interval = "prediction")[,3]
)

# Save  residuals 
Residuals_adjpar <- residuals(model_adjpar)

# Plot observed vs predicted values with confidence intervals



ggplot(predictions_adjpar, aes(x = datagive1$Date[2:159])) +
  geom_line(aes(y = Predicted_adjpar, color = "Predicted"), size = 1) +  # Define color within aes for legend
  geom_ribbon(aes(ymin = Lower_adjpar, ymax = Upper_adjpar, fill = "Confidence Interval"), alpha = 0.2) +  # Define fill within aes for legend
  geom_point(aes(y = lcons_adjpar, color = "Observed"), size = 2) +  # Define color within aes for legend
  scale_color_manual(values = c("Predicted" = "blue", "Observed" = "red")) +  # Manual color values
  scale_fill_manual(values = c("Confidence Interval" = "grey")) +  # Manual fill values
  labs(title = "Observed vs Predicted", x = "Date", y = "Value", 
       color = "Legend", fill = "Legend") +  # Customize legend title
  theme_minimal()


# Plot residuals with 2*sigma^{2}
sigma_value_adjpar = sigma(model_adjpar)
sigma_value_adjpar

ggplot(predictions_adjpar, aes(x = datagive1$Date[2:159], y = Residuals_adjpar)) +
  geom_point(aes(color = "Residuals")) +  # Assign color within aes for legend
  geom_line(aes(color = "Residuals")) +  # Use same color as points to keep them in a single legend item
  geom_hline(aes(yintercept = 2 * sigma_value_adjpar, linetype = "2*sigma"), color = "blue") +  # Map linetype to aes
  geom_hline(aes(yintercept = -2 * sigma_value_adjpar, linetype = "-2*sigma"), color = "gray") +  # Map linetype to aes
  scale_color_manual(values = c("Residuals" = "red")) +  # Manual color values for points and lines
  scale_linetype_manual(values = c("2*sigma" = "dashed", "-2*sigma" = "dashed")) +  # Manual linetype values
  labs(title = "Residuals", x = "Date", y = "residuals", color = "Legend", linetype = "Legend") +  # Custom legend titles
  theme_minimal()





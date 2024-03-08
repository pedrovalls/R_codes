##
# Dynamic Models Typology -Rate of Change
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
# Leading Indicator
##
# Perform linear regression using data from2 to 159 
dlcons=log(datagive1$cons)-log(datagive1$cons1)
dlcons=dlcons[2:159]
dlinc=log(datagive1$inc)-log(datagive1$inc1)
dlinc=dlinc[2:159]

model_TxCre<- lm(dlcons ~ dlinc-1)

# Summary of the model to display regression results
summary(model_TxCre)

# Using White correction
coeftest(model_TxCre, vcov = vcovHC(model_TxCre, type = "HC0"))
#using Newey-West Correction
coeftest(model_TxCre, vcov = vcovHAC(model_TxCre))

num = length(datagive1$cons)-1
num
AIC_TxCre= AIC(model_TxCre)
AIC_TxCre
AIC_TxCre_T=AIC_TxCre/num
AIC_TxCre_T
BIC_TxCre = BIC(model_TxCre)

BIC_TxCre
BIC_TxCre_T=BIC_TxCre/num
BIC_TxCre_T
npar_TxCre = 2
loglik_TxCre = -(1/2)*(AIC_TxCre - 2*npar_TxCre)
loglik_TxCre

# Create a data frame of predicted values and standard errors
predictions_TxCre <- data.frame(
  dlcons,
  Predicted_TxCre = predict(model_TxCre, interval = "prediction")[,1],
  Lower_TxCre = predict(model_TxCre, interval = "prediction")[,2],
  Upper_TxCre = predict(model_TxCre, interval = "prediction")[,3]
)

# Save  residuals 
Residuals_TxCre <- residuals(model_TxCre)

# Plot observed vs predicted values with confidence intervals



ggplot(predictions_TxCre, aes(x = datagive1$Date[2:159])) +
  geom_line(aes(y = Predicted_TxCre, color = "Predicted"), size = 1) +  # Define color within aes for legend
  geom_ribbon(aes(ymin = Lower_TxCre, ymax = Upper_TxCre, fill = "Confidence Interval"), alpha = 0.2) +  # Define fill within aes for legend
  geom_point(aes(y = dlcons, color = "Observed"), size = 2) +  # Define color within aes for legend
  scale_color_manual(values = c("Predicted" = "blue", "Observed" = "red")) +  # Manual color values
  scale_fill_manual(values = c("Confidence Interval" = "grey")) +  # Manual fill values
  labs(title = "Observed vs Predicted", x = "Date", y = "Value", 
       color = "Legend", fill = "Legend") +  # Customize legend title
  theme_minimal()


# Plot residuals with 2*sigma^{2}
sigma_value_TxCre = sigma(model_TxCre)
sigma_value_TxCre

ggplot(predictions_TxCre, aes(x = datagive1$Date[2:159], y = Residuals_TxCre)) +
  geom_point(aes(color = "Residuals")) +  # Assign color within aes for legend
  geom_line(aes(color = "Residuals")) +  # Use same color as points to keep them in a single legend item
  geom_hline(aes(yintercept = 2 * sigma_value_TxCre, linetype = "2*sigma"), color = "blue") +  # Map linetype to aes
  geom_hline(aes(yintercept = -2 * sigma_value_TxCre, linetype = "-2*sigma"), color = "gray") +  # Map linetype to aes
  scale_color_manual(values = c("Residuals" = "red")) +  # Manual color values for points and lines
  scale_linetype_manual(values = c("2*sigma" = "dashed", "-2*sigma" = "dashed")) +  # Manual linetype values
  labs(title = "Residuals", x = "Date", y = "residuals", color = "Legend", linetype = "Legend") +  # Custom legend titles
  theme_minimal()


##
# Now estimate with log(cons_{t}) as dependent variable and log(cons_{t-1}) as a regressor 
# but with fixed coefficient to 1
##

model_new_TxCre <- lm(log(cons) ~ offset(log(cons1)) + I(log(inc) - log(inc1))-1, data = datagive1[2:159,])

summary(model_new_TxCre)


# Using White correction
coeftest(model_new_TxCre, vcov = vcovHC(model_new_TxCre, type = "HC0"))
#using Newey-West Correction
coeftest(model_new_TxCre, vcov = vcovHAC(model_new_TxCre))

num = length(datagive1$cons)-1
num
AIC_TxCre= AIC(model_new_TxCre)
AIC_TxCre
AIC_TxCre_T=AIC_TxCre/num
AIC_TxCre_T
BIC_TxCre = BIC(model_new_TxCre)

BIC_TxCre
BIC_TxCre_T=BIC_TxCre/num
BIC_TxCre_T
npar_TxCre = 2
loglik_TxCre = -(1/2)*(AIC_TxCre - 2*npar_TxCre)
loglik_TxCre


# Create a data frame of predicted values and standard errors
predictions_new_TxCre <- data.frame(
  lcons_new_TxCre = log(datagive1$cons[2:159]),
  Predicted_new_TxCre = predict(model_new_TxCre, interval = "prediction")[,1],
  Lower_new_TxCre = predict(model_new_TxCre, interval = "prediction")[,2],
  Upper_new_TxCre = predict(model_new_TxCre, interval = "prediction")[,3]
)

# Save  residuals 
Residuals_new_TxCre <- residuals(model_new_TxCre)

# Plot observed vs predicted values with confidence intervals



ggplot(predictions_new_TxCre, aes(x = datagive1$Date[2:159])) +
  geom_line(aes(y = Predicted_new_TxCre, color = "Predicted"), size = 1) +  # Define color within aes for legend
  geom_ribbon(aes(ymin = Lower_new_TxCre, ymax = Upper_new_TxCre, fill = "Confidence Interval"), alpha = 0.2) +  # Define fill within aes for legend
  geom_point(aes(y = lcons_new_TxCre, color = "Observed"), size = 2) +  # Define color within aes for legend
  scale_color_manual(values = c("Predicted" = "blue", "Observed" = "red")) +  # Manual color values
  scale_fill_manual(values = c("Confidence Interval" = "grey")) +  # Manual fill values
  labs(title = "Observed vs Predicted", x = "Date", y = "Value", 
       color = "Legend", fill = "Legend") +  # Customize legend title
  theme_minimal()


# Plot residuals with 2*sigma^{2}
sigma_value_new_TxCre = sigma(model_new_TxCre)
sigma_value_new_TxCre

ggplot(predictions_TxCre, aes(x = datagive1$Date[2:159], y = Residuals_new_TxCre)) +
  geom_point(aes(color = "Residuals")) +  # Assign color within aes for legend
  geom_line(aes(color = "Residuals")) +  # Use same color as points to keep them in a single legend item
  geom_hline(aes(yintercept = 2 * sigma_value_new_TxCre, linetype = "2*sigma"), color = "blue") +  # Map linetype to aes
  geom_hline(aes(yintercept = -2 * sigma_value_new_TxCre, linetype = "-2*sigma"), color = "gray") +  # Map linetype to aes
  scale_color_manual(values = c("Residuals" = "red")) +  # Manual color values for points and lines
  scale_linetype_manual(values = c("2*sigma" = "dashed", "-2*sigma" = "dashed")) +  # Manual linetype values
  labs(title = "Residuals", x = "Date", y = "residuals", color = "Legend", linetype = "Legend") +  # Custom legend titles
  theme_minimal()



##
# Dynamic Models Typology  - MCE with Inflta and Inflat1
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

setwd("C:/Users/Pedro/Dropbox/EcoIII2021/Lecture6_adl/R_code")

datagive1 <- read_excel("datagive1.xlsx")
class(datagive1)


# Plotting the beta coefficients
par(mfrow=c(1,3))
plot(datagive1$Date[2:159], datagive1$cons[2:159], type = "l", col= "red", main = "Cons", ylab="cons", xlab="Date")
plot(datagive1$Date[2:159], datagive1$inc[2:159], type = "l", col= "red", main = "Inc", ylab="inc", xlab="Date")
plot(datagive1$Date[2:159], datagive1$inflat[2:159], type = "l", col= "red", main = "Inflat", ylab="inflat", xlab="Date")



##
# Test Omitted Variables: Inflat and Inflat1
##

##
# Unrrestricted model with Inflta and Inflat1
##

model_ADL_EqCM<- lm(log(cons) ~ log(cons1)+log(inc)+log(inc1)+inflat+inflat1, data = datagive1[2:159,])

# Summary of the model to display regression results
summary(model_ADL_Inflat)



# Using White correction
coeftest(model_ADL_Inflat, vcov = vcovHC(model_ADL_Inflat, type = "HC0"))
#using Newey-West Correction
coeftest(model_ADL_Inflat, vcov = vcovHAC(model_ADL_Inflat))

num = length(datagive1$cons)-1
num
AIC_ADL_Inflat= AIC(model_ADL_Inflat)
AIC_ADL_Inflat
AIC_ADL_Inflat_T=AIC_ADL_Inflat/num
AIC_ADL_Inflat_T
BIC_ADL_Inflat = BIC(model_ADL_Inflat)

BIC_ADL_Inflat
BIC_ADL_Inflat_T=BIC_ADL_Inflat/num
BIC_ADL_Inflat_T
npar_ADL_Inflat = length(model_ADL_Inflat$coefficients)
loglik_ADL_Inflat = -(1/2)*(AIC_ADL_Inflat - 2*npar_ADL_Inflat)
loglik_ADL_Inflat



# Create a data frame of predicted values and standard errors
predictions_ADL_Inflat <- data.frame(
  lcons_ADL_Inflat = log(datagive1$cons[2:159]),
  Predicted_ADL_Inflat = predict(model_ADL_Inflat, interval = "prediction")[,1],
  Lower_ADL_Inflat = predict(model_ADL_Inflat, interval = "prediction")[,2],
  Upper_ADL_Inflat = predict(model_ADL_Inflat, interval = "prediction")[,3]
)

# Save  residuals 
Residuals_ADL_Inflat <- residuals(model_ADL_Inflat)

# Plot observed vs predicted values with confidence intervals



ggplot(predictions_ADL_Inflat, aes(x = datagive1$Date[2:159])) +
  geom_line(aes(y = Predicted_ADL_Inflat, color = "Predicted"), size = 1) +  # Define color within aes for legend
  geom_ribbon(aes(ymin = Lower_ADL_Inflat, ymax = Upper_ADL_Inflat, fill = "Confidence Interval"), alpha = 0.2) +  # Define fill within aes for legend
  geom_point(aes(y = lcons_ADL_Inflat, color = "Observed"), size = 2) +  # Define color within aes for legend
  scale_color_manual(values = c("Predicted" = "blue", "Observed" = "red")) +  # Manual color values
  scale_fill_manual(values = c("Confidence Interval" = "grey")) +  # Manual fill values
  labs(title = "Observed vs Predicted - Model with Inflat and Inflat1", x = "Date", y = "Value", 
       color = "Legend", fill = "Legend") +  # Customize legend title
  theme_minimal()


# Plot residuals with 2*sigma^{2}
sigma_value_ADL_Inflat = sigma(model_ADL_Inflat)
sigma_value_ADL_Inflat

ggplot(predictions_ADL_Inflat, aes(x = datagive1$Date[2:159], y = Residuals_ADL_Inflat)) +
  geom_point(aes(color = "Residuals")) +  # Assign color within aes for legend
  geom_line(aes(color = "Residuals")) +  # Use same color as points to keep them in a single legend item
  geom_hline(aes(yintercept = 2 * sigma_value_ADL_Inflat, linetype = "2*sigma"), color = "blue") +  # Map linetype to aes
  geom_hline(aes(yintercept = -2 * sigma_value_ADL_Inflat, linetype = "-2*sigma"), color = "gray") +  # Map linetype to aes
  scale_color_manual(values = c("Residuals" = "red")) +  # Manual color values for points and lines
  scale_linetype_manual(values = c("2*sigma" = "dashed", "-2*sigma" = "dashed")) +  # Manual linetype values
  labs(title = "Residuals Model with Inflat and Inflat1", x = "Date", y = "residuals", color = "Legend", linetype = "Legend") +  # Custom legend titles
  theme_minimal()







###
# Test Ommited Variables
##


test_ommited_chi <-linearHypothesis(model_ADL_Inflat,c("inflat = 0.0","inflat1 = 0"), test ="Chisq")
test_ommited_F <-linearHypothesis(model_ADL_Inflat,c("inflat = 0.0","inflat1 = 0"), test ="F")
test_ommited_chi
test_ommited_F


###
# Test long-run elasticity
##

test_log_run_elasticity_chi_Inflat <-linearHypothesis(model_ADL_Inflat,c("log(cons1)=1","inflat + inflat1 = 0"), test ="Chisq")
test_log_run_elasticity_F_Inflat <-linearHypothesis(model_ADL_Inflat,c("log(cons1)=1","inflat + inflat1 = 0"), test ="F")
test_log_run_elasticity_chi_Inflat
test_log_run_elasticity_F_Inflat
test_log_run_elasticity_chi_inc <-linearHypothesis(model_ADL_Inflat,c("log(cons1) = 1.0","log(inc) = -log(inc1)"), test ="Chisq")
test_log_run_elasticity_F_inc <-linearHypothesis(model_ADL_Inflat,c("log(cons1) = 1.0","log(inc) = -log(inc1)"), test ="F")
test_log_run_elasticity_chi_inc
test_log_run_elasticity_F_inc



##
# Solved static long-run 
##
one_minus_alpha = (1-model_ADL_Inflat$coefficients[2])
one_minus_alpha
new_intercept = model_ADL_Inflat$coefficients[1]/one_minus_alpha
new_intercept
new_beta = (model_ADL_Inflat$coefficients[3]+model_ADL_Inflat$coefficients[4])/one_minus_alpha
new_beta
new_gamma = (model_ADL_Inflat$coefficients[5]+model_ADL_Inflat$coefficients[6])/one_minus_alpha
new_gamma

EqCM = log(datagive1$cons[2:159])-new_intercept-new_beta*log(datagive1$inc[2:159])-new_gamma*datagive1$inflat[2:159]

par(mfrow=c(1,1))

plot(datagive1$Date[2:159], EqCM, type = "l", col= "red", main = "EqCM", ylab="ECM", xlab="Date")



##
# unit roots test for EqCM
##

adf_ECM <- ur.df(EqCM, type="trend", selectlags = "BIC")
summary(adf_ECM)


##
# EqCM - Equilibrium Correction Model
##

DLCons = log(datagive1$cons[2:159]) - log(datagive1$cons1[2:159])
DLCons1 = DLCons[2:158]
DLInc = log(datagive1$inc[2:159]) - log(datagive1$inc1[2:159])
DLInc1=DLInc[2:158]
DInflat = datagive1$inflat[2:159] - datagive1$inflat1[2:159]
DInflat1 = DInflat[2:158]
EqCM1=EqCM[2:158]


model_ADL_EqCM<- lm(DLCons1 ~ DLInc1+DInflat1+EqCM1-1)
summary(model_ADL_EqCM)



# Using White correction
coeftest(model_ADL_EqCM, vcov = vcovHC(model_ADL_EqCM, type = "HC0"))
#using Newey-West Correction
coeftest(model_ADL_EqCM, vcov = vcovHAC(model_ADL_EqCM))

num = length(datagive1$cons)-1
num
AIC_ADL_EqCM= AIC(model_ADL_EqCM)
AIC_ADL_EqCM
AIC_ADL_EqCM_T=AIC_ADL_EqCM/num
AIC_ADL_EqCM_T
BIC_ADL_EqCM = BIC(model_ADL_EqCM)

BIC_ADL_EqCM
BIC_ADL_EqCM_T=BIC_ADL_EqCM/num
BIC_ADL_EqCM_T
npar_ADL_EqCM = length(model_ADL_EqCM$coefficients)
loglik_ADL_EqCM = -(1/2)*(AIC_ADL_EqCM - 2*npar_ADL_EqCM)
loglik_ADL_EqCM


# Create a data frame of predicted values and standard errors
predictions_ADL_EqCM <- data.frame(
  Dlcons_ADL_EqCM = DLCons1,
  Predicted_ADL_EqCM = predict(model_ADL_EqCM, interval = "prediction")[,1],
  Lower_ADL_EqCM = predict(model_ADL_EqCM, interval = "prediction")[,2],
  Upper_ADL_EqCM = predict(model_ADL_EqCM, interval = "prediction")[,3]
)

# Save  residuals 
Residuals_ADL_EqCM <- residuals(model_ADL_EqCM)



# Plot observed vs predicted values with confidence intervals



ggplot(predictions_ADL_EqCM, aes(x = datagive1$Date[3:159])) +
  geom_line(aes(y = Predicted_ADL_EqCM, color = "Predicted"), size = 1) +  # Define color within aes for legend
  geom_ribbon(aes(ymin = Lower_ADL_EqCM, ymax = Upper_ADL_EqCM, fill = "Confidence Interval"), alpha = 0.2) +  # Define fill within aes for legend
  geom_point(aes(y = Dlcons_ADL_EqCM, color = "Observed"), size = 2) +  # Define color within aes for legend
  scale_color_manual(values = c("Predicted" = "blue", "Observed" = "red")) +  # Manual color values
  scale_fill_manual(values = c("Confidence Interval" = "grey")) +  # Manual fill values
  labs(title = "Observed vs Predicted - Model EqCM", x = "Date", y = "Value", 
       color = "Legend", fill = "Legend") +  # Customize legend title
  theme_minimal()


# Plot residuals with 2*sigma^{2}
sigma_value_ADL_EqCM = sigma(model_ADL_EqCM)
sigma_value_ADL_EqCM

ggplot(predictions_ADL_EqCM, aes(x = datagive1$Date[3:159], y = Residuals_ADL_EqCM)) +
  geom_point(aes(color = "Residuals")) +  # Assign color within aes for legend
  geom_line(aes(color = "Residuals")) +  # Use same color as points to keep them in a single legend item
  geom_hline(aes(yintercept = 2 * sigma_value_ADL_EqCM, linetype = "2*sigma"), color = "blue") +  # Map linetype to aes
  geom_hline(aes(yintercept = -2 * sigma_value_ADL_EqCM, linetype = "-2*sigma"), color = "gray") +  # Map linetype to aes
  scale_color_manual(values = c("Residuals" = "red")) +  # Manual color values for points and lines
  scale_linetype_manual(values = c("2*sigma" = "dashed", "-2*sigma" = "dashed")) +  # Manual linetype values
  labs(title = "Residuals Model EqCM", x = "Date", y = "residuals", color = "Legend", linetype = "Legend") +  # Custom legend titles
  theme_minimal()


##
# Test Breusch Godfrey 
## 

bgtest_ADLEqCM = bgtest(DLCons1 ~ DLInc1+DInflat1+EqCM1-1, order = 4)
bgtest_ADLEqCM

Res_ADLEqCM = ts(Residuals_ADL_EqCM)


u_bg = ts.intersect(Res_ADLEqCM, stats::lag(Res_ADLEqCM,-1), stats::lag(Res_ADLEqCM,-2), stats::lag(Res_ADLEqCM,-3), stats::lag(Res_ADLEqCM,-4)) 


bg_test_ADLEqCM = lm(u_bg[,1] ~ u_bg[,2]+u_bg[,3]+u_bg[,4]+u_bg[,5])
summary(bg_test_ADLEqCM)

summary(bg_test_ADLEqCM)$r.squared
(summary(bg_test_ADLEqCM)$r.squared)*153
t_stat_bg = (summary(bg_test_ADLEqCM)$r.squared)*153
1-pchisq(t_stat_bg,df = 4)


##
# Test Breusch `Godfrey `Pagan which is the same as White´s test
##

bptest_ADLEqCM = bptest(model_ADL_EqCM, ~  I(DLInc1^2)+I(DInflat1^2)+ I(EqCM1^2))
bptest_ADLEqCM       


##
# Jarque & Bera test for residulas
##

jbtest_ADLEqCM = jarque.bera.test(Residuals_ADL_EqCM )
jbtest_ADLEqCM

##
# RESET test
##
reset_ADLEqCM = resettest(model_ADL_EqCM, power =3, type = "regressor")
reset_ADLEqCM


##
# ARCH test
##

Res_SQ_ADLEqCM =ts((Res_SQ_ADLEqCM))
u = ts.intersect(Res_SQ_ADLEqCM, stats::lag(Res_SQ_ADLEqCM,-1), stats::lag(Res_SQ_ADLEqCM,-2), stats::lag(Res_SQ_ADLEqCM,-3), stats::lag(Res_SQ_ADLEqCM,-4)) 



arch_test_ADLEqCM = lm(u[,1] ~ u[,2]+u[,3]+u[,4]+u[,5])
summary(arch_test_ADLEqCM)

t_stat_arch = (summary(arch_test_ADLEqCM)$r.squared)*153
t_stat_arch
p_value_arch = 1-pchisq(t_stat_arch,df = 4)
p_value_arch

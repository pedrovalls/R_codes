##
# Modelos GARCH usando fGARCH
#


library(fGarch) # Loading the fGarch package for GARCH modeling
library(forecast)
# Load the data

library(readxl)
data <- read_excel("C:/Users/Pedro/Dropbox/EcoIII2021/Lecture_Volatilidade_Univariada/Vol_R/DADOS_BR.xlsx")
#data <- read_excel("C:/Users/pedro.valls/Dropbox/EcoIII2021/Lecture_Volatilidade_Univariada/Vol_R/DADOS_BR.xlsx")
# data <- read.csv("C:/Users/Pedro/Dropbox/Topicos_em_Financas_2022/garch_model/dados_br.csv") # Adjust file extension if needed

# Define the variable dlibovm as the difference between dlibovc_sb_per and its mean
data$dlibov = (as.numeric(data$DLIBOVC))*100 
data$dlibovm[2:5880] <- as.numeric(data$dlibov[2:5880]) - mean(as.numeric(data$dlibov[2:5880]))
data$dlibovm
# Estimate a NARCH(1) model
model_narch <- garchFit(formula = ~ garch(1, 0), data = data$dlibovm[2:5880], cond.dist = "norm")

# Print model summary
print(summary(model_narch))

# Save model summary to a text file
#write.table(summary(model), file = "C:/Users/pedro.valls/Dropbox/EcoIII2021/Lecture_Volatilidade_Univariada/Vol_R/narch.txt", quote = FALSE)

# Calculate volatility
data$volatilidade_n_arch[2:5880] <- sqrt(model_narch@h.t)
  

# Plot volatility
par(mfrow=c(1,1))
plot(data$Date[2:5880],data$volatilidade_n_arch[2:5880], type = "l", main = "Volatilidade para RCPIBOV usando N-ARCH", xlab = "Time", ylab = "Volatility")

# Check the adequacy of the model
residuals_narch <- residuals(model_narch, standardize = TRUE)

# Plot FAC and FACP of standardized residuals
par(mfrow=c(1,2))
Acf(residuals_narch, lag.max = 12, main = "ACF of Standardized Residuals")
Pacf(residuals_narch, lag.max = 12, main = "PACF of Standardized Residuals")

# Square the residuals
data$res_n_arch_st_sq[2:5880] <- residuals_narch^2

# Plot FAC and FACP of squared residuals
par(mfrow=c(1,2))
Acf(data$res_n_arch_st_sq[2:5880], lag.max = 12, main = "ACF of Squared Standardized Residuals")
Pacf(data$res_n_arch_st_sq[2:5880], lag.max = 12, main = "PACF of Squared Standardized Residuals")



# Estimate a tARCH(1) model
model_tarch <- garchFit(formula = ~ garch(1, 0), data = data$dlibovm[2:5880], cond.dist = "std")

# Print model summary
print(summary(model_tarch))

# Save model summary to a text file
#write.table(summary(model), file = "C:/Users/pedro.valls/Dropbox/EcoIII2021/Lecture_Volatilidade_Univariada/Vol_R/narch.txt", quote = FALSE)

# Calculate volatility
data$volatilidade_t_arch[2:5880] <- sqrt(model_tarch@h.t)


# Plot volatility
par(mfrow=c(1,1))
plot(data$Date[2:5880],data$volatilidade_t_arch[2:5880], type = "l", main = "Volatilidade para RCPIBOV usando t-ARCH", xlab = "Time", ylab = "Volatility")

# Check the adequacy of the model
residuals_tarch <- residuals(model_tarch, standardize = TRUE)

# Plot FAC and FACP of standardized residuals
par(mfrow=c(1,2))
Acf(residuals_tarch, lag.max = 12, main = "ACF of Standardized Residuals")
Pacf(residuals_tarch, lag.max = 12, main = "PACF of Standardized Residuals")

# Square the residuals
data$res_t_arch_st_sq[2:5880] <- residuals_tarch^2

# Plot FAC and FACP of squared residuals
par(mfrow=c(1,2))
Acf(data$res_t_arch_st_sq[2:5880], lag.max = 12, main = "ACF of Squared Standardized Residuals")
Pacf(data$res_t_arch_st_sq[2:5880], lag.max = 12, main = "PACF of Squared Standardized Residuals")



# Estimate a GEDARCH(1) model
model_gedarch <- garchFit(formula = ~ garch(1, 0), data = data$dlibovm[2:5880], cond.dist = "ged", algorithm = "lbfgsb+nm" )

# Print model summary
print(summary(model_gedarch))

# Save model summary to a text file
#write.table(summary(model), file = "C:/Users/pedro.valls/Dropbox/EcoIII2021/Lecture_Volatilidade_Univariada/Vol_R/narch.txt", quote = FALSE)

# Calculate volatility
data$volatilidade_ged_arch[2:5880] <- sqrt(model_gedarch@h.t)


# Plot volatility
par(mfrow=c(1,1))
plot(data$Date[2:5880],data$volatilidade_ged_arch[2:5880], type = "l", main = "Volatilidade para RCPIBOV usando ged-ARCH", xlab = "Time", ylab = "Volatility")

# Check the adequacy of the model
residuals_gedarch <- residuals(model_gedarch, standardize = TRUE)

# Plot FAC and FACP of standardized residuals
par(mfrow=c(1,2))
Acf(residuals_gedarch, lag.max = 12, main = "ACF of Standardized Residuals")
Pacf(residuals_gedarch, lag.max = 12, main = "PACF of Standardized Residuals")

# Square the residuals
data$res_ged_arch_st_sq[2:5880] <- residuals_gedarch^2

# Plot FAC and FACP of squared residuals
par(mfrow=c(1,2))
Acf(data$res_ged_arch_st_sq[2:5880], lag.max = 12, main = "ACF of Squared Standardized Residuals")
Pacf(data$res_ged_arch_st_sq[2:5880], lag.max = 12, main = "PACF of Squared Standardized Residuals")






# Estimate a SSTDARCH(1) model
model_sstdarch <- garchFit(formula = ~ garch(1, 0), data = data$dlibovm[2:5880], cond.dist = "sstd", algorithm = "lbfgsb+nm" )

# Print model summary
print(summary(model_sstdarch))

# Save model summary to a text file
#write.table(summary(model), file = "C:/Users/pedro.valls/Dropbox/EcoIII2021/Lecture_Volatilidade_Univariada/Vol_R/narch.txt", quote = FALSE)

# Calculate volatility
data$volatilidade_sstd_arch[2:5880] <- sqrt(model_sstdarch@h.t)


# Plot volatility
par(mfrow=c(1,1))
plot(data$Date[2:5880],data$volatilidade_sstd_arch[2:5880], type = "l", main = "Volatilidade para RCPIBOV usando sstd-ARCH", xlab = "Time", ylab = "Volatility")

# Check the adequacy of the model
residuals_sstdarch <- residuals(model_sstdarch, standardize = TRUE)

# Plot FAC and FACP of standardized residuals
par(mfrow=c(1,2))
Acf(residuals_sstdarch, lag.max = 12, main = "ACF of Standardized Residuals")
Pacf(residuals_sstdarch, lag.max = 12, main = "PACF of Standardized Residuals")

# Square the residuals
data$res_sstd_arch_st_sq[2:5880] <- residuals_sstdarch^2

# Plot FAC and FACP of squared residuals
par(mfrow=c(1,2))
Acf(data$res_sstd_arch_st_sq[2:5880], lag.max = 12, main = "ACF of Squared Standardized Residuals")
Pacf(data$res_sstd_arch_st_sq[2:5880], lag.max = 12, main = "PACF of Squared Standardized Residuals")



##
# Infornation Criteria
##

model_narch@fit$ics
model_tarch@fit$ics
model_gedarch@fit$ics
model_sstdarch@fit$ics



##
# Loglik
##

(-1)*model_narch@fit$llh
(-1)*model_tarch@fit$llh
(-1)*model_gedarch@fit$llh
(-1)*model_sstdarch@fit$llh


##
# plot vol for four models
##

# Plot volatility
par(mfrow=c(2,2))
plot(data$Date[2:5880],data$volatilidade_n_arch[2:5880], type = "l", main = "Volatilidade para RCPIBOV usando N-ARCH", xlab = "Time", ylab = "Volatility")
plot(data$Date[2:5880],data$volatilidade_t_arch[2:5880], type = "l", main = "Volatilidade para RCPIBOV usando t-ARCH", xlab = "Time", ylab = "Volatility")
plot(data$Date[2:5880],data$volatilidade_ged_arch[2:5880], type = "l", main = "Volatilidade para RCPIBOV usando ged-ARCH", xlab = "Time", ylab = "Volatility")
plot(data$Date[2:5880],data$volatilidade_sstd_arch[2:5880], type = "l", main = "Volatilidade para RCPIBOV usando sstd-ARCH", xlab = "Time", ylab = "Volatility")

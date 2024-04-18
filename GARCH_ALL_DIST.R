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
# Estimate a Ngarch(1) model
model_ngarch <- garchFit(formula = ~ garch(1, 1), data = data$dlibovm[2:5880], cond.dist = "norm")

# Print model summary
print(summary(model_ngarch))

# Save model summary to a text file
#write.table(summary(model), file = "C:/Users/pedro.valls/Dropbox/EcoIII2021/Lecture_Volatilidade_Univariada/Vol_R/ngarch.txt", quote = FALSE)

# Calculate volatility
data$volatilidade_n_garch[2:5880] <- sqrt(model_ngarch@h.t)


# Plot volatility
par(mfrow=c(1,1))
plot(data$Date[2:5880],data$volatilidade_n_garch[2:5880], type = "l", main = "Volatilidade para RCPIBOV usando N-garch", xlab = "Time", ylab = "Volatility")

# Check the adequacy of the model
residuals_ngarch <- residuals(model_ngarch, standardize = TRUE)

# Plot FAC and FACP of standardized residuals
par(mfrow=c(1,2))
Acf(residuals_ngarch, lag.max = 12, main = "ACF of Standardized Residuals")
Pacf(residuals_ngarch, lag.max = 12, main = "PACF of Standardized Residuals")

# Square the residuals
data$res_n_garch_st_sq[2:5880] <- residuals_ngarch^2

# Plot FAC and FACP of squared residuals
par(mfrow=c(1,2))
Acf(data$res_n_garch_st_sq[2:5880], lag.max = 12, main = "ACF of Squared Standardized Residuals")
Pacf(data$res_n_garch_st_sq[2:5880], lag.max = 12, main = "PACF of Squared Standardized Residuals")



# Estimate a tgarch(1) model
model_tgarch <- garchFit(formula = ~ garch(1, 1), data = data$dlibovm[2:5880], cond.dist = "std")

# Print model summary
print(summary(model_tgarch))

# Save model summary to a text file
#write.table(summary(model), file = "C:/Users/pedro.valls/Dropbox/EcoIII2021/Lecture_Volatilidade_Univariada/Vol_R/ngarch.txt", quote = FALSE)

# Calculate volatility
data$volatilidade_t_garch[2:5880] <- sqrt(model_tgarch@h.t)


# Plot volatility
par(mfrow=c(1,1))
plot(data$Date[2:5880],data$volatilidade_t_garch[2:5880], type = "l", main = "Volatilidade para RCPIBOV usando t-garch", xlab = "Time", ylab = "Volatility")

# Check the adequacy of the model
residuals_tgarch <- residuals(model_tgarch, standardize = TRUE)

# Plot FAC and FACP of standardized residuals
par(mfrow=c(1,2))
Acf(residuals_tgarch, lag.max = 12, main = "ACF of Standardized Residuals")
Pacf(residuals_tgarch, lag.max = 12, main = "PACF of Standardized Residuals")

# Square the residuals
data$res_t_garch_st_sq[2:5880] <- residuals_tgarch^2

# Plot FAC and FACP of squared residuals
par(mfrow=c(1,2))
Acf(data$res_t_garch_st_sq[2:5880], lag.max = 12, main = "ACF of Squared Standardized Residuals")
Pacf(data$res_t_garch_st_sq[2:5880], lag.max = 12, main = "PACF of Squared Standardized Residuals")



# Estimate a GEDgarch(1) model
model_gedgarch <- garchFit(formula = ~ garch(1, 1), data = data$dlibovm[2:5880], cond.dist = "ged", algorithm = "lbfgsb+nm" )

# Print model summary
print(summary(model_gedgarch))

# Save model summary to a text file
#write.table(summary(model), file = "C:/Users/pedro.valls/Dropbox/EcoIII2021/Lecture_Volatilidade_Univariada/Vol_R/ngarch.txt", quote = FALSE)

# Calculate volatility
data$volatilidade_ged_garch[2:5880] <- sqrt(model_gedgarch@h.t)


# Plot volatility
par(mfrow=c(1,1))
plot(data$Date[2:5880],data$volatilidade_ged_garch[2:5880], type = "l", main = "Volatilidade para RCPIBOV usando ged-garch", xlab = "Time", ylab = "Volatility")

# Check the adequacy of the model
residuals_gedgarch <- residuals(model_gedgarch, standardize = TRUE)

# Plot FAC and FACP of standardized residuals
par(mfrow=c(1,2))
Acf(residuals_gedgarch, lag.max = 12, main = "ACF of Standardized Residuals")
Pacf(residuals_gedgarch, lag.max = 12, main = "PACF of Standardized Residuals")

# Square the residuals
data$res_ged_garch_st_sq[2:5880] <- residuals_gedgarch^2

# Plot FAC and FACP of squared residuals
par(mfrow=c(1,2))
Acf(data$res_ged_garch_st_sq[2:5880], lag.max = 12, main = "ACF of Squared Standardized Residuals")
Pacf(data$res_ged_garch_st_sq[2:5880], lag.max = 12, main = "PACF of Squared Standardized Residuals")






# Estimate a SSTDgarch(1) model
model_sstdgarch <- garchFit(formula = ~ garch(1, 1), data = data$dlibovm[2:5880], cond.dist = "sstd", algorithm = "lbfgsb+nm" )

# Print model summary
print(summary(model_sstdgarch))

# Save model summary to a text file
#write.table(summary(model), file = "C:/Users/pedro.valls/Dropbox/EcoIII2021/Lecture_Volatilidade_Univariada/Vol_R/ngarch.txt", quote = FALSE)

# Calculate volatility
data$volatilidade_sstd_garch[2:5880] <- sqrt(model_sstdgarch@h.t)


# Plot volatility
par(mfrow=c(1,1))
plot(data$Date[2:5880],data$volatilidade_sstd_garch[2:5880], type = "l", main = "Volatilidade para RCPIBOV usando sstd-garch", xlab = "Time", ylab = "Volatility")

# Check the adequacy of the model
residuals_sstdgarch <- residuals(model_sstdgarch, standardize = TRUE)

# Plot FAC and FACP of standardized residuals
par(mfrow=c(1,2))
Acf(residuals_sstdgarch, lag.max = 12, main = "ACF of Standardized Residuals")
Pacf(residuals_sstdgarch, lag.max = 12, main = "PACF of Standardized Residuals")

# Square the residuals
data$res_sstd_garch_st_sq[2:5880] <- residuals_sstdgarch^2

# Plot FAC and FACP of squared residuals
par(mfrow=c(1,2))
Acf(data$res_sstd_garch_st_sq[2:5880], lag.max = 12, main = "ACF of Squared Standardized Residuals")
Pacf(data$res_sstd_garch_st_sq[2:5880], lag.max = 12, main = "PACF of Squared Standardized Residuals")



##
# Infornation Criteria
##

model_ngarch@fit$ics
model_tgarch@fit$ics
model_gedgarch@fit$ics
model_sstdgarch@fit$ics



##
# Loglik
##

(-1)*model_ngarch@fit$llh
(-1)*model_tgarch@fit$llh
(-1)*model_gedgarch@fit$llh
(-1)*model_sstdgarch@fit$llh


##
# plot vol for four models
##

# Plot volatility
par(mfrow=c(2,2))
plot(data$Date[2:5880],data$volatilidade_n_garch[2:5880], type = "l", main = "Volatilidade para RCPIBOV usando N-garch", xlab = "Time", ylab = "Volatility")
plot(data$Date[2:5880],data$volatilidade_t_garch[2:5880], type = "l", main = "Volatilidade para RCPIBOV usando t-garch", xlab = "Time", ylab = "Volatility")
plot(data$Date[2:5880],data$volatilidade_ged_garch[2:5880], type = "l", main = "Volatilidade para RCPIBOV usando ged-garch", xlab = "Time", ylab = "Volatility")
plot(data$Date[2:5880],data$volatilidade_sstd_garch[2:5880], type = "l", main = "Volatilidade para RCPIBOV usando sstd-garch", xlab = "Time", ylab = "Volatility")
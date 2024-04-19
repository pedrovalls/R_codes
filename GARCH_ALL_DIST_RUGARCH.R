##
# Modelos GARCH usando rugarch
#

library(rugarch)  # Loading the rugarch package for GARCH modeling
library(forecast)
library(readxl)

# Load the data
data <- read_excel("C:/Users/Pedro/Dropbox/EcoIII2021/Lecture_Volatilidade_Univariada/Vol_R/DADOS_BR.xlsx")

# Define the variable dlibovm as the difference between dlibovc_sb_per and its mean
data$dlibov = (as.numeric(data$DLIBOVC))*100 
data$dlibovm <- as.numeric(data$dlibov) - mean(as.numeric(data$dlibov), na.rm = TRUE)

# Specify GARCH(1,1) model using rugarch
spec_ngarch <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
                          mean.model = list(armaOrder = c(0, 0), include.mean = TRUE),
                          distribution.model = "norm")

# Fit the model
fit_ngarch <- ugarchfit(spec = spec_ngarch, data = data$dlibovm[2:5880], solver = "hybrid")

# Print model summary
show(fit_ngarch)

# Calculate volatility
data$volatilidade_n_garch[2:5880] <- sigma(fit_ngarch)

# Plot volatility
par(mfrow=c(1,1))
plot(data$Date[2:5880], data$volatilidade_n_garch[2:5880], type = "l", main = "Volatilidade para RCPIBOV usando N-garch", xlab = "Time", ylab = "Volatility")

# Check the adequacy of the model for residuals
residuals_ngarch <- residuals(fit_ngarch, standardize = TRUE)

# Plot FAC and FACP of standardized residuals
par(mfrow=c(1,2))
Acf(residuals_ngarch, lag.max = 12, main = "ACF of Standardized Residuals")
Pacf(residuals_ngarch, lag.max = 12, main = "PACF of Standardized Residuals")


# Check the adequacy of the model for squared residuals
residuals_ngarch_sq = residuals_ngarch^2

# Plot FAC and FACP of standardized Squared residuals
par(mfrow=c(1,2))
Acf(residuals_ngarch_sq, lag.max = 12, main = "ACF of Standardized Squared Residuals")
Pacf(residuals_ngarch_sq, lag.max = 12, main = "PACF of Standardized Squared Residuals")




# Other models (t-GARCH, GED-GARCH, and SSTD-GARCH) following a similar pattern:
# Specify the model t-GARCH 
spec_tgarch <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
                          mean.model = list(armaOrder = c(0, 0), include.mean = TRUE),
                          distribution.model = "std")  # "std" for Student's t

# Fit the model
fit_tgarch <- ugarchfit(spec = spec_tgarch, data = data$dlibovm[2:5880], solver = "hybrid")

# Print model summary
show(fit_tgarch)

# Calculate and plot volatility as before
data$volatilidade_t_garch[2:5880] <- sigma(fit_tgarch)

par(mfrow=c(1,1))
plot(data$Date[2:5880],data$volatilidade_t_garch[2:5880], type = "l", main = "Volatilidade para RCPIBOV usando t-garch", xlab = "Time", ylab = "Volatility")


# Check the adequacy of the model for residuals
residuals_tgarch <- residuals(fit_tgarch, standardize = TRUE)

# Plot FAC and FACP of standardized residuals
par(mfrow=c(1,2))
Acf(residuals_tgarch, lag.max = 12, main = "ACF of Standardized Residuals")
Pacf(residuals_tgarch, lag.max = 12, main = "PACF of Standardized Residuals")


# Check the adequacy of the model for squared residuals
residuals_tgarch_sq = residuals_tgarch^2

# Plot FAC and FACP of standardized residuals
par(mfrow=c(1,2))
Acf(residuals_tgarch_sq, lag.max = 12, main = "ACF of Standardized Squared Residuals")
Pacf(residuals_tgarch_sq, lag.max = 12, main = "PACF of Standardized Squared Residuals")







# Specify the model GED-GARCH 
spec_gedgarch <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
                          mean.model = list(armaOrder = c(0, 0), include.mean = TRUE),
                          distribution.model = "ged")  # "ged" for GED-GARCH

# Fit the model
fit_gedgarch <- ugarchfit(spec = spec_gedgarch, data = data$dlibovm[2:5880], solver = "hybrid")

# Print model summary
show(fit_gedgarch)

# Calculate and plot volatility as before
data$volatilidade_ged_garch[2:5880] <- sigma(fit_gedgarch)

par(mfrow=c(1,1))
plot(data$Date[2:5880],data$volatilidade_ged_garch[2:5880], type = "l", main = "Volatilidade para RCPIBOV usando ged-garch", xlab = "Time", ylab = "Volatility")


# Check the adequacy of the model for residuals
residuals_gedgarch <- residuals(fit_gedgarch, standardize = TRUE)

# Plot FAC and FACP of standardized residuals
par(mfrow=c(1,2))
Acf(residuals_gedgarch, lag.max = 12, main = "ACF of Standardized Residuals")
Pacf(residuals_gedgarch, lag.max = 12, main = "PACF of Standardized Residuals")


# Check the adequacy of the model for squared residuals
residuals_gedgarch_sq = residuals_gedgarch^2

# Plot FAC and FACP of standardized residuals
par(mfrow=c(1,2))
Acf(residuals_gedgarch_sq, lag.max = 12, main = "ACF of Standardized Squared Residuals")
Pacf(residuals_gedgarch_sq, lag.max = 12, main = "PACF of Standardized Squared Residuals")



# Specify the model SSTD-GARCH models using "sstd" in distribution.model
spec_sstdgarch <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
                            mean.model = list(armaOrder = c(0, 0), include.mean = TRUE),
                            distribution.model = "sstd")  # "sstd" for GED-GARCH

# Fit the model
fit_sstdgarch <- ugarchfit(spec = spec_sstdgarch, data = data$dlibovm[2:5880], solver = "hybrid")

# Print model summary
show(fit_sstdgarch)

# Calculate and plot volatility as before
data$volatilidade_sstd_garch[2:5880] <- sigma(fit_sstdgarch)

par(mfrow=c(1,1))
plot(data$Date[2:5880],data$volatilidade_sstd_garch[2:5880], type = "l", main = "SSTD-GARCH", xlab = "Time", ylab = "Volatility")


# Check the adequacy of the model for residuals
residuals_sstdgarch <- residuals(fit_sstdgarch, standardize = TRUE)

# Plot FAC and FACP of standardized residuals
par(mfrow=c(1,2))
Acf(residuals_sstdgarch, lag.max = 12, main = "ACF of Standardized Residuals")
Pacf(residuals_sstdgarch, lag.max = 12, main = "PACF of Standardized Residuals")


# Check the adequacy of the model for squared residuals
residuals_sstdgarch_sq = residuals_sstdgarch^2

# Plot FAC and FACP of standardized residuals
par(mfrow=c(1,2))
Acf(residuals_sstdgarch_sq, lag.max = 12, main = "ACF of Standardized Squared Residuals")
Pacf(residuals_sstdgarch_sq, lag.max = 12, main = "PACF of Standardized Squared Residuals")




# Information Criteria for all models
info_criteria_ngarch <- infocriteria(fit_ngarch)
info_criteria_ngarch
info_criteria_tgarch <- infocriteria(fit_tgarch)
info_criteria_tgarch
info_criteria_gedgarch <- infocriteria(fit_gedgarch )
info_criteria_gedgarch
info_criteria_sstdgarch <- infocriteria(fit_sstdgarch )
info_criteria_sstdgarch

# Log Likelihood for all models
loglik_ngarch <- fit_ngarch@fit$LLH
loglik_ngarch
loglik_tgarch <- fit_tgarch@fit$LLH
loglik_tgarch
loglik_gedgarch <- fit_gedgarch@fit$LLH
loglik_gedgarch
loglik_sstdgarch <- fit_sstdgarch@fit$LLH
loglik_sstdgarch

# Plot volatility comparisons
par(mfrow=c(2,2))
plot(data$Date[2:5880],data$volatilidade_n_garch[2:5880], type = "l", main = "N-GARCH", xlab = "Time", ylab = "Volatility")
plot(data$Date[2:5880],data$volatilidade_t_garch[2:5880], type = "l", main = "t-GARCH", xlab = "Time", ylab = "Volatility")
plot(data$Date[2:5880],data$volatilidade_ged_garch[2:5880], type = "l", main = "GED-GARCH", xlab = "Time", ylab = "Volatility")
plot(data$Date[2:5880],data$volatilidade_sstd_garch[2:5880], type = "l", main = "SSTD-GARCH", xlab = "Time", ylab = "Volatility")
# Load necessary libraries
library(vars)
library(MASS)
library(ggplot2)
library(gridExtra)


# Set seed for reproducibility
set.seed(123456789)

# Parameters
a11 <- 0.4
a12 <- 0.2
a21 <- 0.6
a22 <- 0.1
m1 <- 0
m2 <- 0

# Covariance matrix
cov <- matrix(c(1, 0.5, 0.5, 4), nrow = 2)

# Generate random errors
u <- mvrnorm(n = 200, mu = c(0, 0), Sigma = cov)

# Initialize series
y1 <- numeric(200)
y2 <- numeric(200)

# Generate data
for (i in 2:200) {
  y1[i] <- m1 + a11 * y1[i - 1] + a12 * y2[i - 1] + u[i, 1]
  y2[i] <- m2 + a21 * y1[i - 1] + a22 * y2[i - 1] + u[i, 2]
}

# Combine into a data frame
data <- data.frame(y1, y2)

# Fit VAR model
var_model <- VAR(data, p = 1, type = "const")

# Extract coefficients and residual covariance matrix
a <- coef(var_model)
c1 <- a$y1[3]
c2 <- a$y2[3]
vcv <- cov(residuals(var_model))

# Forecast
forecast_horizon <- 20
y1_forecast <- numeric(forecast_horizon)
y2_forecast <- numeric(forecast_horizon)
vy1 <- numeric(forecast_horizon)
vy2 <- numeric(forecast_horizon)
y1_forecast[1] <- y1[200]
y2_forecast[1] <- y2[200]

for (j in 2:forecast_horizon) {
  y1_forecast[j] <- c1 + a$y1[1] * y1_forecast[j - 1] + a$y1[2] * y2_forecast[j - 1]
  y2_forecast[j] <- c2 + a$y2[1] * y1_forecast[j - 1] + a$y2[2] * y2_forecast[j - 1]
}

# Variance calculation
al <- diag(2)
vcvs <- matrix(0, 2, 2)
fvcv <- vcv
vy1[1] <- vcv[1, 1]
vy2[1] <- vcv[2, 2]
coef_var_model <- matrix(0,2,2)
coef_var_model[1,1] = a$y1[1]
coef_var_model[1,2] = a$y1[2]
coef_var_model[2,1] = a$y2[1]
coef_var_model[2,2] = a$y2[2]

for (k in 2:forecast_horizon) {
  al <- al %*% coef_var_model
  vcvs <- vcvs + al %*% vcv %*% t(al)
  fvcv <- vcv + vcvs
  vy1[k] <- fvcv[1, 1]
  vy2[k] <- fvcv[2, 2]
}

# Standard deviation of forecast errors
sdfy1 <- sqrt(vy1)
sdfy2 <- sqrt(vy2)

# Forecast intervals
li1 <- y1_forecast - 2 * sdfy1
ls1 <- y1_forecast + 2 * sdfy1
li2 <- y2_forecast - 2 * sdfy2
ls2 <- y2_forecast + 2 * sdfy2

# Create data frame for plotting
forecast_data <- data.frame(
  time = 201:220,
  y1_forecast = y1_forecast,
  y2_forecast = y2_forecast,
  li1 = li1,
  ls1 = ls1,
  li2 = li2,
  ls2 = ls2
)

# Plot forecasts and intervals
ggplot(forecast_data, aes(x = time)) +
  geom_line(aes(y = y1_forecast), color = "blue") +
  geom_line(aes(y = li1), linetype = "dashed", color = "red") +
  geom_line(aes(y = ls1), linetype = "dashed", color = "red") +
  labs(title = "Forecast for y1 with 95% Confidence Interval",
       x = "Time", y = "y1") +
  theme_minimal()

ggplot(forecast_data, aes(x = time)) +
  geom_line(aes(y = y2_forecast), color = "blue") +
  geom_line(aes(y = li2), linetype = "dashed", color = "red") +
  geom_line(aes(y = ls2), linetype = "dashed", color = "red") +
  labs(title = "Forecast for y2 with 95% Confidence Interval",
       x = "Time", y = "y2") +
  theme_minimal()


##
# side by side
##

# Plot forecasts and intervals
p1 <- ggplot(forecast_data, aes(x = time)) +
  geom_line(aes(y = y1_forecast), color = "blue") +
  geom_line(aes(y = li1), linetype = "dashed", color = "red") +
  geom_line(aes(y = ls1), linetype = "dashed", color = "red") +
  labs(title = "Forecast for y1 with 95% Confidence Interval",
       x = "Time", y = "y1") +
  theme_minimal()

p2 <- ggplot(forecast_data, aes(x = time)) +
  geom_line(aes(y = y2_forecast), color = "blue") +
  geom_line(aes(y = li2), linetype = "dashed", color = "red") +
  geom_line(aes(y = ls2), linetype = "dashed", color = "red") +
  labs(title = "Forecast for y2 with 95% Confidence Interval",
       x = "Time", y = "y2") +
  theme_minimal()

# Arrange plots side by side
grid.arrange(p1, p2, ncol = 2)
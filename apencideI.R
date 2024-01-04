# Set the number of observations
n <- 1000

# Initialize series
u <- rep(0, n)
du <- rep(0, n)
d2u <- rep(0, n)
ds2u <- rep(0, n)

# Generate normally distributed random numbers for 'u'
u <- rnorm(n)

# First difference of u
du <- diff(u)

# Second difference of u
d2u <- diff(u, differences = 2)

# Seasonal difference of u with a seasonal lag of 12 

seasonal_lag <- 12
dsu <- diff(u, lag = seasonal_lag)

# Plotting the series and differences, if needed
par(mfrow=c(2,2))
plot(u, type = "l", main = "Series U", xlab = "Time", ylab = "U")
plot(du, type = "l", main = "First Difference of U", xlab = "Time", ylab = "DU")
plot(d2u, type = "l", main = "Second Difference of U", xlab = "Time", ylab = "D2U")
plot(dsu, type = "l", main = "Seasonal Difference of U", xlab = "Time", ylab = "DS2U")

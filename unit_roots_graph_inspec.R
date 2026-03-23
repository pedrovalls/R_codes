############################
## Unit Roots - Graphical Inspection
############################

# Clear workspace
rm(list = ls())

# -----------------------------------
# Load packages
# -----------------------------------
load_package <- function(x) {
  x <- as.character(match.call()[[2]])
  if (!require(x, character.only = TRUE)) {
    install.packages(pkgs = x, repos = "https://cloud.r-project.org")
    require(x, character.only = TRUE)
  }
}

load_package("tseries")
load_package("forecast")

library(forecast)
library(tseries)

set.seed(123456)

# -----------------------------------
# Simulate the series
# -----------------------------------
n <- 1000
time_series <- 1:n

u   <- rep(0, n)
x1  <- rep(0, n)
x2  <- rep(0, n)
x3  <- rep(0, n)
x4  <- rep(0, n)
x5  <- rep(0, n)
x51 <- rep(0, n)

psi1 <- 0.95
psi2 <- 1.00
psi3 <- 1.05
psi4 <- 0.50

m      <- 1
delta0 <- 10
delta1 <- 0.5

x4[1]  <- delta0 + delta1
x5[1]  <- delta0 + delta1
x51[1] <- delta0 + delta1

for (t in 2:n) {
  u[t] <- rnorm(1)

  # Stationary around mean m
  x1[t] <- m + psi1 * (x1[t - 1] - m) + u[t]

  # Unit root
  x2[t] <- m + psi2 * (x2[t - 1] - m) + u[t]

  # Explosive
  x3[t] <- m + psi3 * (x3[t - 1] - m) + u[t]

  # Trend-stationary
  x4[t] <- delta0 + delta1 * t +
    psi4 * (x4[t - 1] - delta0 - delta1 * (t - 1)) + u[t]

  # Hybrid / badly behaved trending process
  x5[t] <- delta0 + delta1 * t + psi2 * x5[t - 1] + u[t]

  # Unit root around deterministic trend
  x51[t] <- delta0 + delta1 * t +
    psi2 * (x51[t - 1] - delta0 - delta1 * (t - 1)) + u[t]
}

# -----------------------------------
# Helper function: find crossings
# -----------------------------------
find_crossings <- function(series, level) {
  idx <- which((series[-1] - level) * (series[-length(series)] - level) < 0)
  return(idx + 1)
}

# Levels for horizontal lines
h1 <- -5.7
h2 <-  7.4

# Crossings for x1 and x2
cross_x1_h1 <- find_crossings(x1, h1)
cross_x1_h2 <- find_crossings(x1, h2)

cross_x2_h1 <- find_crossings(x2, h1)
cross_x2_h2 <- find_crossings(x2, h2)

# -----------------------------------
# Plot all series
# -----------------------------------
par(mfrow = c(3, 2), mar = c(3, 3, 2, 1))

plot(time_series, x1, type = "l", main = "x1: Stationary",
     xlab = "Time", ylab = "x1")
abline(h = m, lty = 2)

plot(time_series, x2, type = "l", main = "x2: Unit Root",
     xlab = "Time", ylab = "x2")

plot(time_series, x3, type = "l", main = "x3: Explosive",
     xlab = "Time", ylab = "x3")

plot(time_series, x4, type = "l", main = "x4: Trend-Stationary",
     xlab = "Time", ylab = "x4")
abline(a = delta0, b = delta1, lty = 2)

plot(time_series, x5, type = "l", main = "x5: Hybrid",
     xlab = "Time", ylab = "x5")

plot(time_series, x51, type = "l", main = "x51: Unit Root + Trend",
     xlab = "Time", ylab = "x51")
abline(a = delta0, b = delta1, lty = 2)

# -----------------------------------
# Special plot for x1 with horizontal lines and crossings
# -----------------------------------
par(mfrow = c(1, 1), mar = c(4, 4, 2, 1))

plot(time_series, x1, type = "l", lwd = 1.5,
     main = "x1 with Horizontal Lines and Crossings",
     xlab = "Time", ylab = "x1")

abline(h = h1, col = "red",  lty = 2, lwd = 2)
abline(h = h2, col = "blue", lty = 2, lwd = 2)

points(time_series[cross_x1_h1], x1[cross_x1_h1], col = "red",  pch = 16)
points(time_series[cross_x1_h2], x1[cross_x1_h2], col = "blue", pch = 16)

text(x = 900, y = h1, labels = "-5.7", col = "red",  pos = 3)
text(x = 900, y = h2, labels = "7.4",  col = "blue", pos = 3)

legend("topleft",
       legend = c("x1", "y = -5.7", "y = 7.4", "Crossings at -5.7", "Crossings at 7.4"),
       col    = c("black", "red", "blue", "red", "blue"),
       lty    = c(1, 2, 2, NA, NA),
       pch    = c(NA, NA, NA, 16, 16),
       bty    = "n")

# -----------------------------------
# Comparison plot: x1 versus x2
# -----------------------------------
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

plot(time_series, x1, type = "l", lwd = 1.5,
     main = "x1: Stationary",
     xlab = "Time", ylab = "x1")
abline(h = h1, col = "red",  lty = 2, lwd = 2)
abline(h = h2, col = "blue", lty = 2, lwd = 2)
points(time_series[cross_x1_h1], x1[cross_x1_h1], col = "red",  pch = 16)
points(time_series[cross_x1_h2], x1[cross_x1_h2], col = "blue", pch = 16)
legend("topleft",
       legend = c(
         paste("Crossings at -5.7:", length(cross_x1_h1)),
         paste("Crossings at 7.4:",  length(cross_x1_h2))
       ),
       bty = "n")

plot(time_series, x2, type = "l", lwd = 1.5,
     main = "x2: Unit Root",
     xlab = "Time", ylab = "x2")
abline(h = h1, col = "red",  lty = 2, lwd = 2)
abline(h = h2, col = "blue", lty = 2, lwd = 2)
points(time_series[cross_x2_h1], x2[cross_x2_h1], col = "red",  pch = 16)
points(time_series[cross_x2_h2], x2[cross_x2_h2], col = "blue", pch = 16)
legend("topleft",
       legend = c(
         paste("Crossings at -5.7:", length(cross_x2_h1)),
         paste("Crossings at 7.4:",  length(cross_x2_h2))
       ),
       bty = "n")

# -----------------------------------
# Print crossing counts
# -----------------------------------
cat("\n=============================\n")
cat("Crossing counts\n")
cat("=============================\n")
cat("x1 crossings at -5.7:", length(cross_x1_h1), "\n")
cat("x1 crossings at  7.4:", length(cross_x1_h2), "\n")
cat("x2 crossings at -5.7:", length(cross_x2_h1), "\n")
cat("x2 crossings at  7.4:", length(cross_x2_h2), "\n")

# -----------------------------------
# ACF and PACF for x2 and x3
# -----------------------------------
par(mfrow = c(2, 2))
Acf(x2, lag.max = 12, main = "Acf x2")
Pacf(x2, lag.max = 12, main = "Pacf x2")
Acf(x3, lag.max = 12, main = "Acf x3")
Pacf(x3, lag.max = 12, main = "Pacf x3")

# -----------------------------------
# ACF and PACF for x4
# -----------------------------------
par(mfrow = c(1, 2))
Acf(x4, lag.max = 12, main = "Acf x4")
Pacf(x4, lag.max = 12, main = "Pacf x4")

# -----------------------------------
# ACF and PACF for x5
# -----------------------------------
par(mfrow = c(1, 2))
Acf(x5, lag.max = 12, main = "Acf x5")
Pacf(x5, lag.max = 12, main = "Pacf x5")

# -----------------------------------
# ACF and PACF for x51
# -----------------------------------
par(mfrow = c(1, 2))
Acf(x51, lag.max = 12, main = "Acf x51")
Pacf(x51, lag.max = 12, main = "Pacf x51")

# -----------------------------------
# Detrending x4
# -----------------------------------
x4_dt <- lm(x4 ~ time_series)
summary(x4_dt)
x4_detrend <- x4_dt$residuals

par(mfrow = c(1, 1))
plot(time_series, x4_detrend, type = "l", main = "",
     xlab = "Time", ylab = "x4 detrended")

# -----------------------------------
# Detrending x2
# -----------------------------------
x2_dt <- lm(x2 ~ time_series)
summary(x2_dt)
x2_detrend <- x2_dt$residuals

par(mfrow = c(1, 1))
plot(time_series, x2_detrend, type = "l", main = "",
     xlab = "Time", ylab = "x2 detrended")

# -----------------------------------
# First-differencing x2
# -----------------------------------
dx2 <- diff(x2, lag = 1)
par(mfrow = c(1, 1))
plot(time_series[2:n], dx2, type = "l", main = "",
     xlab = "Time", ylab = "dx2")

# -----------------------------------
# First-differencing x4
# -----------------------------------
dx4 <- diff(x4, lag = 1)
par(mfrow = c(1, 2))
plot(time_series[2:n], dx4, type = "l", main = "",
     xlab = "Time", ylab = "dx4")
plot(time_series[2:n], x4_detrend[2:n], type = "l", main = "",
     xlab = "Time", ylab = "x4_detrend")

# -----------------------------------
# ACF and PACF for dx4 and x4_detrend
# -----------------------------------
par(mfrow = c(2, 2))
Acf(dx4, lag.max = 12, main = "Acf dx4")
Pacf(dx4, lag.max = 12, main = "Pacf dx4")
Acf(x4_detrend, lag.max = 12, main = "Acf x4hat")
Pacf(x4_detrend, lag.max = 12, main = "Pacf x4hat")

# -----------------------------------
# ACF and PACF for dx2, x2_detrend and x2
# -----------------------------------
par(mfrow = c(3, 2))
Acf(dx2, lag.max = 12, main = "Acf dx2")
Pacf(dx2, lag.max = 12, main = "Pacf dx2")
Acf(x2_detrend, lag.max = 12, main = "Acf x2hat")
Pacf(x2_detrend, lag.max = 12, main = "Pacf x2hat")
Acf(x2, lag.max = 12, main = "Acf x2")
Pacf(x2, lag.max = 12, main = "Pacf x2")

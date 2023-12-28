# Load necessary libraries
library(tseries)

# Set parameters
num_realizations <- 1000
num_points <- 100

# Initialize matrices
ma_error <- matrix(nrow = num_points, ncol = num_realizations)
mma <- matrix(nrow = num_points, ncol = num_realizations)

# Populate ma_error with random normal variables
set.seed(123) # for reproducibility
ma_error <- matrix(rnorm(num_points * num_realizations), num_points, num_realizations)

# Generating the MA(1) process
for (j in 1:num_realizations) {
  for (i in 2:num_points) {
    mma[i, j] <- ma_error[i, j] - 0.5 * ma_error[i - 1, j]
  }
}

# Create lists to store series
tma_series <- list()
sma_series <- list()

# Extracting series for tma and sma
for (k in 1:num_points) {
  tma_series[[k]] <- ma_error[k, ]
}

for (k in 1:num_realizations) {
  sma_series[[k]] <- mma[, k]
}

# Now tma_series and sma_series contain the required series

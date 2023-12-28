# Setting up the environment
n <- 1000  # Number of realizations
T <- 100   # Length of each realization
phi <- 0.5 # AR parameter

# Initialize a matrix to store the AR process realizations
mar <- matrix(0, nrow = T, ncol = n)

# Generate the AR process realizations
set.seed(123)  # For reproducibility
for (i in 1:n) {
  for (j in 2:T) {
    mar[j, i] <- phi * mar[j - 1, i] + rnorm(1)
  }
}

# Extracting series at each time point
tar <- list()
for (k in 1:T) {
  tar[[k]] <- mar[k, ]
}

# Extracting each realization of the AR process
sar <- list()
for (k in 1:n) {
  sar[[k]] <- mar[, k]
}

# Optionally, convert lists to time series objects or data frames for further analysis or plotting


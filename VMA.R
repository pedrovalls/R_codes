#
# Script to simulate MVA(1)
#


# Load necessary libraries
library(MASS)  # for mvrnorm and cov2cor
library(dplyr)  # for data manipulation

# Set seed for reproducibility
set.seed(123456789)

# Parameters
n <- 201
A <- matrix(c(0.55, 0.1, -0.1, 0.3), nrow=2, byrow=TRUE)  # A matrix
m <- c(0, 0)  # constant vector m

# Covariance matrix
cov_matrix <- matrix(c(0.3, 0.2, 0.2, 0.2), nrow=2, byrow=TRUE)
sigma <- cov2cor(cov_matrix)  # Convert to correlation matrix for mvrnorm

# Generate multivariate normal errors
e <- mvrnorm(n, mu=c(0,0), Sigma=sigma)

# Construct y_t
y <- matrix(0, nrow=n, ncol=2)
for (i in 2:n) {
  y[i,] <- m + e[i,] + A %*% e[i-1,]
}

# Convert to time-series data frame
df <- as.data.frame(y)
colnames(df) <- c("y1", "y2")

# Display the first few rows of the time-series data
head(df)

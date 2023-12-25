# Load necessary packages and library
#install.packages("tidyverse")
#install.packages("ggplot2") 
#install.packages("reshape2")

library(ggplot2)
library(reshape2) 



# 
# Stochastic process
#
#
# fixed the random seed for the stochastic process $y(t,w)$
#
set.seed(123456) 
#
# loop tyo simulate the stochastic process  y(t,w) standard normal 
#

# Define the dimensions
n_rows <- 100
n_cols <- 100

# Create a matrix of random numbers
# rnorm(n_rows * n_cols) generates n_rows*n_cols random numbers
# matrix() reshapes these numbers into a n_rows x n_cols matrix
y <- matrix(rnorm(n_rows * n_cols), nrow = n_rows, ncol = n_cols)



#
# random seed for the stochastic process z(t,w)
#
set.seed(654321) 

# Define the dimensions of the matrix
n_rows <- 100
n_cols <- 100

# Generate the matrix
# rnorm(n_rows * n_cols) creates a vector of n_rows*n_cols random numbers from a normal distribution
# 10 + sqrt(5) * ... scales these numbers according to your formula
z <- matrix(10 + sqrt(5) * rnorm(n_rows * n_cols), nrow = n_rows, ncol = n_cols)



#
# random seed for the stochastic process  u(t,w)
#

set.seed(123) 

#
# simulate the stochastic process  u(t,w) U[-$\phi$, $\phi$]
#

# Define the dimensions of the matrix
n_rows <- 100
n_cols <- 100

# Calculate the constant pi
pi <- acos(-1)

# Generate the u1 matrix with random numbers from uniform distribution
u1 <- matrix(runif(n_rows * n_cols), nrow = n_rows, ncol = n_cols)

# Transform the values in u1 to create the u matrix
# The transformation is -pi + pi * u1
u <- -pi + pi * u1



#
#  simulate the stochastic process  x(t,w) defined by  (1)
#

# Define the dimensions of the matrices
n_rows <- 100
n_cols <- 100

# Create a matrix of row indices
i_matrix <- matrix(rep(1:n_rows, each = n_cols), nrow = n_rows, ncol = n_cols)

# Calculate the x matrix
# The operation is vectorized over all elements
x <- y * cos(z + i_matrix * u)

# You can print the matrix or write it to a file
# Uncomment the line below to print the matrix
# print(x)


#
# Transpose the matrix x and assign it to xl
#
xl <- t(x)

#
# extract the first three trows of matrix xl 
#


first_three_columns <- xl[, 1:3]





# Assuming first_three_columns is already defined

# Create a data frame suitable for ggplot
df <- data.frame(Index = 1:nrow(first_three_columns),
                 x1 = first_three_columns[,1],
                 x2 = first_three_columns[,2],
                 x3 = first_three_columns[,3])

# Melt the data frame to long format for ggplot
long_df <- reshape2::melt(df, id.vars = 'Index')

# Plotting using ggplot
ggplot(long_df, aes(x = Index, y = value, color = variable)) +
  geom_line() +
  theme_minimal() +
  labs(title = "Line Graph of First Three Columns",
       x = "Index",
       y = "Value",
       color = "Column")



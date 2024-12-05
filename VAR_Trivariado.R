#
# Script to simulate VAR(1) trivariado 
#
##
# clean workspace
##
rm(list = ls()) 
##
# free up memory space
##
gc(reset = TRUE)


# Load package using a function load_package-----------------------------------------------------------------
load_package<-function(x){
  x<-as.character(match.call()[[2]])
  if (!require(x,character.only=TRUE)){
    install.packages(pkgs=x,repos="http://cran.r-project.org")
    require(x,character.only=TRUE)
  }
}

load_package(" MASS")
load_package("dplyr")
load_package("forecast")
load_package("xtable")
load_package("readxl")
load_package("stats")
load_package("ggplot2")


# Load necessary libraries
library(MASS)  # for mvrnorm and cov2cor
library(dplyr)  # for data manipulation
# Load necessary library
library(ggplot2)
# Melt the data frame to long format for easier plotting with ggplot2
library(tidyr)
library(forecast)

# Set seed for reproducibility
set.seed(123456789)

# Parameters
n <- 201
A <- matrix(c(1.0,1.0,0.9,0.0,0.5,0.5,0.0,0.0,0.3), nrow=3, byrow=TRUE)  # A matrix
A = t(A)
#A <- matrix(c(1.0,1.0,0.9,1.0,0.5,0.5,0.9,0.5,0.3), nrow=3, byrow=TRUE)  # A matrix
m <- c(0, 0,0)  # constant vector m

# Covariance matrix
#cov_matrix <- matrix(c(0.3, 0.2, 0.2, 0.2), nrow=2, byrow=TRUE)
cov_matrix <- matrix(c(1.0, 0.2, 0.2, 0.2,0.5,0.2,0.2,0.2,1.0), nrow=3, byrow=TRUE)
sigma <- cov2cor(cov_matrix)  # Convert to correlation matrix for mvrnorm

# Generate multivariate normal errors
e <- mvrnorm(n, mu=c(0,0,0), Sigma=sigma)

# Construct y_t
y <- matrix(0, nrow=n, ncol=3)
y[1,1] = e[1,1]
y[1,2] = e[1,2]
y[1,3] = e[1,3]
for (i in 2:n) {
  y[i,] <- m + e[i,] + A %*% y[i-1,]
}

# Convert to time-series data frame
df <- as.data.frame(y)
colnames(df) <- c("y1", "y2", "y3")

# Display the first few rows of the time-series data
head(df)



# Assuming 'df' is your data frame from the previous code
df$time <- 1:nrow(df)  # Create a time variable for plotting

# Plotting both components of y
p1 <- ggplot(df, aes(x = time)) + 
  geom_line(aes(y = y1), color = "blue") +
  ggtitle("Component y1 over Time") +
  xlab("Time") +
  ylab("y1") +
  theme_minimal()

p2 <- ggplot(df, aes(x = time)) + 
  geom_line(aes(y = y2), color = "red") +
  ggtitle("Component y2 over Time") +
  xlab("Time") +
  ylab("y2") +
  theme_minimal()

p3 <- ggplot(df, aes(x = time)) + 
  geom_line(aes(y = y3), color = "darkblue") +
  ggtitle("Component y3 over Time") +
  xlab("Time") +
  ylab("y3") +
  theme_minimal()
# Print the plots
print(p1)
print(p2)
print(p3)

# plotting in the same plot

df_long <- pivot_longer(df, cols = c(y1, y2,y3), names_to = "variable", values_to = "value")

# Plotting both components of y on the same graph
combined_plot <- ggplot(df_long, aes(x = time, y = value, color = variable)) + 
  geom_line() +
  ggtitle("Components y1, y2 and y3 over Time") +
  xlab("Time") +
  ylab("Value") +
  theme_minimal() +
  scale_color_manual(values = c("blue", "red", "darkblue"), labels = c("y1", "y2","y3"))

# Print the combined plot
print(combined_plot)





# Plotting both components of y, each on its own panel
side_by_side_plot <- ggplot(df_long, aes(x = time, y = value,color=variable)) +
  geom_line() +
  facet_wrap(~variable, scales = "free_y") +  # Create two panels, free_y allows independent y scales
  ggtitle("Comparison of Components y1, y2 and y3 over Time") +
  xlab("Time") +
  ylab("Value") +
  theme_minimal()+
  scale_color_manual(values = c("blue", "red","darkblue"), labels = c("y1", "y2","y3"))
# Print the side-by-side plot
print(side_by_side_plot)




##
# two cointegrated relationship
# y2 - y1 and y3 - y1
##

z1 <- y[,2]-y[,1]

z1_coint <- ggplot(df, aes(x = time)) + 
  geom_line(aes(y = z1), color = "red") +
  ggtitle("") +
  xlab("Time") +
  ylab("z1") +
  theme_minimal()
par(mfrow=c(1,1))
print(z1_coint)

z2 <- y[,3]-y[,1]

z2_coint <- ggplot(df, aes(x = time)) + 
  geom_line(aes(y = z2), color = "blue") +
  ggtitle("") +
  xlab("Time") +
  ylab("z2") +
  theme_minimal()
par(mfrow=c(1,1))
print(z2_coint)

par(mfrow=c(2,1))
plot(df$time,z1, xlab="Time", ylab="z1", type="l", col = "red" )
plot(df$time,z2, xlab="Time", ylab="z2", type= "l",col = "blue" )





##
# Acf and Pacf for y1
##
par(mfrow=c(2,1))
Acf(y[,1], lag.max = 24)
Pacf(y[,1], lag.max = 24)

##
# Acf and Pacf for y2
##
par(mfrow=c(2,1))
Acf(y[,2], lag.max = 24)
Pacf(y[,2], lag.max = 24)


##
# Acf and Pacf for y3
##
par(mfrow=c(2,1))
Acf(y[,3], lag.max = 24)
Pacf(y[,3], lag.max = 24)




##
# Ccf cross-correlation for y1 and y2
##
par(mfrow=c(1,1))
Ccf(y[,1],y[,2], lag.max=24, type = "correlation")


##
# Ccf cross-correlation for y1 and y3
##
par(mfrow=c(1,1))
Ccf(y[,1],y[,3], lag.max=24, type = "correlation")


##
# Ccf cross-correlation for y2 and y3
##
par(mfrow=c(1,1))
Ccf(y[,2],y[,3], lag.max=24, type = "correlation")




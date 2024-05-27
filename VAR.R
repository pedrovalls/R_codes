#
# Script to simulate VAR(1)
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
load_package("vars")
load_package("moments")


# Load necessary libraries
library(MASS)  # for mvrnorm and cov2cor
library(dplyr)  # for data manipulation
# Load necessary library
library(ggplot2)
# Melt the data frame to long format for easier plotting with ggplot2
library(tidyr)
library(forecast)
library(vars)
library(moments)

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
  y[i,] <- m + e[i,] + A %*% y[i-1,]
}

# Convert to time-series data frame
df <- as.data.frame(y)
colnames(df) <- c("y1", "y2")

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
# Print the plots
print(p1)
print(p2)

# plotting in the same plot

df_long <- pivot_longer(df, cols = c(y1, y2), names_to = "variable", values_to = "value")

# Plotting both components of y on the same graph
combined_plot <- ggplot(df_long, aes(x = time, y = value, color = variable)) + 
  geom_line() +
  ggtitle("Components y1 and y2 over Time") +
  xlab("Time") +
  ylab("Value") +
  theme_minimal() +
  scale_color_manual(values = c("blue", "red"), labels = c("y1", "y2"))

# Print the combined plot
print(combined_plot)





# Plotting both components of y, each on its own panel
side_by_side_plot <- ggplot(df_long, aes(x = time, y = value,color=variable)) +
  geom_line() +
  facet_wrap(~variable, scales = "free_y") +  # Create two panels, free_y allows independent y scales
  ggtitle("Comparison of Components y1 and y2 over Time") +
  xlab("Time") +
  ylab("Value") +
  theme_minimal()+
  scale_color_manual(values = c("blue", "red"), labels = c("y1", "y2"))
# Print the side-by-side plot
print(side_by_side_plot)




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
# Ccf cross-correlation for y1 and y2
##
par(mfrow=c(1,1))
Ccf(y[,1],y[,2], lag.max=24, type = "correlation")

##
# Similar model for univariate AR
##

y1 <- matrix(0, nrow=n, ncol=1)
y2 <- matrix(0, nrow=n, ncol=1)
for (i in 2:n) {
  y1[i] <- m[1] + e[i,1] + A[1,1] * y1[i-1]
  y2[i] <- m[2] + e[i,1] + A[2,2] * y2[i-1]
}


# Convert to time-series data frame
df_uni <- as.data.frame(cbind(y1,y2))
colnames(df_uni) <- c("y1", "y2")
df_uni$time <- 1:nrow(df)  # Create a time variable for plotting

# Display the first few rows of the time-series data
head(df_uni)

df_long_uni <- pivot_longer(df_uni, cols = c(y1, y2), names_to = "variable", values_to = "value")


# Plotting both components of y, each on its own panel
side_by_side_plot_uni <- ggplot(df_long_uni, aes(x = time, y = value,color=variable)) +
  geom_line() +
  facet_wrap(~variable, scales = "free_y") +  # Create two panels, free_y allows independent y scales
  ggtitle("Comparison of univariate components y1 and y2 over Time") +
  xlab("Time") +
  ylab("Value") +
  theme_minimal()+
  scale_color_manual(values = c("blue", "red"), labels = c("y1", "y2"))
# Print the side-by-side plot
print(side_by_side_plot_uni)



##
# Acf and Pacf for y1
##
par(mfrow=c(2,1))
Acf(y1, lag.max = 24, main ="Acf for y1")
Pacf(y1, lag.max = 24, main = "Pacf for y1")

##
# Acf and Pacf for y2
##
par(mfrow=c(2,1))
Acf(y2, lag.max = 24, main ="Acf for y2")
Pacf(y2, lag.max = 24, main = "Pacf for y2")


write.csv(y, file= "C:/Users/Pedro/Dropbox/EcoIII2021/Lecture7_var_vec/Script_R/y.csv" )

##
# Select the order of VAR
##
new_df <- as.data.frame(y)
colnames(new_df) <- c("y1", "y2") 
resultado_varselect <- VARselect(new_df, lag.max = 10, type = "both")
resultado_varselect
ordem_optima <- resultado_varselect$selection["SC(n)"]
ordem_optima

##
# Fixing the sample size
##
resultado_varselect_fix <- VARselect(new_df[11:201,], lag.max = 10, type = "both")
resultado_varselect_fix
ordem_optima_fix <- resultado_varselect_fix$selection["SC(n)"]
ordem_optima_fix

modelo_var1 <- VAR(new_df, p=1, type="both",ic = "BIC")
modelo_var1
summary(modelo_var1)

##
# VAR(1) without constant and linear trend
##

modelo_var2 <- VAR(new_df, p=1, type="none",ic = "BIC")
modelo_var2
summary(modelo_var2)

##
# Test for serial correlation
##
Serial2 <- serial.test(modelo_var2, lags.pt = 16, type = "PT.asymptotic")
Serial2
plot(Serial2, names="y1")
plot(Serial2, names = "y2")

##
# Test for normality
##
Normal2 <- normality.test(modelo_var2, multivariate.only = FALSE)
Normal2

##
# heteroscedasticity of ARCH type
##

Arch2 <- arch.test(modelo_var2, lags.single=16, multivariate.only = FALSE)
Arch2


##
# Show AR roots (eigenvalues of the companion matrix)
##
par(mfrow=c(1,1))
ar_roots <- roots(modelo_var2)
time = cbind(1,2)
plot(time,ar_roots, type = "b")

# Create points on the unit circle
theta <- seq(0, 2 * pi, length.out = 100)
x <- cos(theta)
y <- sin(theta)

# Points to be added
points_x <- c(0.573, 0.28)
points_y <- c(0, 0)

# Plot the unit circle
plot(x, y, type = 'l', asp = 1, main = "Inverse Roots of AR Characteristic Polynomial", xlab = "X", ylab = "Y")
abline(h = 0, col = "blue", lty = 2)  # x-axis
abline(v = 0, col = "blue", lty = 2)  # y-axis

# Add points to the plot
points(points_x, points_y, col = "red", pch = 19, cex = 1.5)


##
# Using ggplot
##


# Create a data frame with points on the unit circle
theta <- seq(0, 2 * pi, length.out = 100)
circle <- data.frame(
  x = cos(theta),
  y = sin(theta)
)

# Points to be added
points_df <- data.frame(
  x = c(0.573, 0.28),
  y = c(0, 0)
)

# Plot the unit circle and points using ggplot2
ggplot(circle, aes(x, y)) +
  geom_path() +
  geom_point(data = points_df, aes(x, y), color = "red", size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "blue") +  # x-axis
  geom_vline(xintercept = 0, linetype = "dashed", color = "blue") +  # y-axis
  coord_fixed() +  # Ensure the aspect ratio is 1:1
  labs(title = "Inverse Roots of AR Characteristic Polynomial",
       x = "X",
       y = "Y") +
  theme_minimal()

##
# Grangre Caulatity test
##
##
# H_{0} y1 do not Granger cause y2
#
causality(modelo_var2, cause=c("y1"))


# H_{0} y2 do not Granger cause y1
#
causality(modelo_var2, cause=c("y2"))



##
# IRF orthogonal 
##
irf_modelo_var2_ort <- irf(modelo_var2, n.ahead = 10, otho = TRUE)
plot(irf_modelo_var2_ort)

##
# IRF cumulative
##
irf_modelo_var2_cumul <- irf(modelo_var2, n.ahead = 10, ortho = FALSE, cumulative = TRUE)
plot(irf_modelo_var2_cumul)
##
# IRF ortho and cumulative equal to FALSE 
##

irf_modelo_var2 <- irf(modelo_var2, n.ahead = 10, ortho = FALSE, cumulative = FALSE)
plot(irf_modelo_var2)

#'
#' ## Variance Decomposition
#' 

fevd(modelo_var2, n.ahead=10)





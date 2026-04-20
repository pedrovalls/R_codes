load_package<-function(x){
  x<-as.character(match.call()[[2]])
  if (!require(x,character.only=TRUE)){
    install.packages(pkgs=x,repos="http://cran.r-project.org")
    require(x,character.only=TRUE)
  }
}

load_package('readxl')
load_package('ggplot2')
load_package('dplyr')



# ============================================
# Packages
# ============================================
library(readxl)
library(ggplot2)
library(dplyr)

# ============================================
# 0. Directory
# ============================================
setwd("C:/Users/Pedro/Dropbox/Time_Series_School_of_Methods/lecture1")

# Theoretical correlogram for AR(1) with phi = 0.5
phi     <- 0.5
max_lag <- 20

# Theoretical autocorrelations: rho_k = phi^k
lags <- 1:max_lag
rho  <- phi^lags

# Data frame for ggplot
df <- data.frame(lag = lags, rho = rho)

# Plot
g1 <- ggplot(df, aes(x = lag, y = rho)) +
  geom_hline(yintercept = 0, linewidth = 0.5, color = "black") +
  geom_segment(aes(x = lag, xend = lag, y = 0, yend = rho),
               color = "steelblue", linewidth = 1) +
  geom_point(color = "steelblue", size = 2.5) +
  scale_x_continuous(breaks = 1:max_lag) +
  scale_y_continuous(limits = c(-1, 1)) +
  labs(
    title    = expression("Theoretical ACF - AR(1) with " ~ phi == 0.5),
    x        = "Lag",
    y        = expression(rho[k])
  ) +
  theme_bw()
  
# Save Theoretical Correlogram for phi = 0.5  
ggsave("The_Correl_AR_P05.pdf", plot = g1, width = 8, height = 5)


# Theoretical correlogram for AR(1) with phi = -0.5
phi     <- - 0.5
max_lag <- 20

# Theoretical autocorrelations: rho_k = phi^k
lags <- 1:max_lag
rho  <- phi^lags

# Data frame for ggplot
df <- data.frame(lag = lags, rho = rho)

# Plot
g2 <- ggplot(df, aes(x = lag, y = rho)) +
  geom_hline(yintercept = 0, linewidth = 0.5, color = "black") +
  geom_segment(aes(x = lag, xend = lag, y = 0, yend = rho),
               color = "steelblue", linewidth = 1) +
  geom_point(color = "steelblue", size = 2.5) +
  scale_x_continuous(breaks = 1:max_lag) +
  scale_y_continuous(limits = c(-1, 1)) +
  labs(
    title    = expression("Theoretical ACF - AR(1) with " ~ phi == -0.5),
    x        = "Lag",
    y        = expression(rho[k])
  ) +
  theme_bw()
# Save Theoretical Correlogram for phi = -0.5  
ggsave("The_Correl_AR_N05.pdf", plot = g2, width = 8, height = 5)  
  
  
# Theoretical correlogram for AR(1) with phi = 0.99
phi     <-  0.99
max_lag <- 20

# Theoretical autocorrelations: rho_k = phi^k
lags <- 1:max_lag
rho  <- phi^lags

# Data frame for ggplot
df <- data.frame(lag = lags, rho = rho)

# Plot
g3 <- ggplot(df, aes(x = lag, y = rho)) +
  geom_hline(yintercept = 0, linewidth = 0.5, color = "black") +
  geom_segment(aes(x = lag, xend = lag, y = 0, yend = rho),
               color = "steelblue", linewidth = 1) +
  geom_point(color = "steelblue", size = 2.5) +
  scale_x_continuous(breaks = 1:max_lag) +
  scale_y_continuous(limits = c(-1, 1)) +
  labs(
    title    = expression("Theoretical ACF - AR(1) with " ~ phi == 0.99),
    x        = "Lag",
    y        = expression(rho[k])
  ) +
  theme_bw()  
  
# Save Theoretical Correlogram for phi = 0.99  
ggsave("The_Correl_AR_P099.pdf", plot = g3, width = 8, height = 5)  
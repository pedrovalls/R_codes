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
load_package("vars")
load_package("tseries")
load_package("tsDyn")

# Load necessary libraries
library(MASS)  # for mvrnorm and cov2cor
library(dplyr)  # for data manipulation
# Load necessary library
library(ggplot2)
# Melt the data frame to long format for easier plotting with ggplot2
library(tidyr)
library(forecast)
library(vars)
library(tseries)
library(tsDyn)


# Set seed for reproducibility
 set.seed(123456789)
#set.seed(12345)


# Parameters
n <- 201
#A <- matrix(c(1.0,0.0,0.0,0.0,1.0,1.0,0.0,0.0,0.5), nrow=3, byrow=TRUE)  # A matrix
A <- matrix(c(1.0,1.0,0.9,0.0,0.5,0.5,0.0,0.0,0.3), nrow=3, byrow=TRUE)  # A matrix
A = t(A)
A
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
#y[,2]=y[,1]+y[,2]
#y[,3] = y[,1]+y[,2]+y[,3]
#y[,1]=y[,1]+y[,3]
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

eq_y1_y3 <- lm(df$y1 ~ df$y3 -1 )
summary(eq_y1_y3)

eq_y2_y1 <- lm(df$y2 ~ df$y1 -1)
summary(eq_y2_y1)

eq_y3_y2 <- lm(df$y3 ~ df$y2 -1)
summary(eq_y3_y2)


##
# two cointegrated relationship
# y2 - y1 and y3 - y2
##

z21 <- y[,2]-1.96544*y[,1]

z21_coint <- ggplot(df, aes(x = time)) + 
  geom_line(aes(y = z21), color = "red") +
  ggtitle("") +
  xlab("Time") +
  ylab("z1") +
  theme_minimal()
par(mfrow=c(1,1))
print(z21_coint)

z32 <- y[,3]-1.35505*y[,2]

z32_coint <- ggplot(df, aes(x = time)) + 
  geom_line(aes(y = z32), color = "blue") +
  ggtitle("") +
  xlab("Time") +
  ylab("z32") +
  theme_minimal()
par(mfrow=c(1,1))
print(z32_coint)

par(mfrow=c(2,1))
plot(df$time,z21, xlab="Time", ylab="z1", type="l", col = "red" )
plot(df$time,z32, xlab="Time", ylab="z2", type= "l",col = "blue" )





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


write.csv(y, file= "C:/Users/Pedro/Dropbox/EcoIII2021/Lecture7_var_vec/Script_R/y_tri_2casoNew.csv" )
#write.csv(y, file= "C:/Users/pedro.valls/Dropbox/EcoIII2021/Lecture7_var_vec/Script_R/y_tri_2casoN.csv" )


EqCM_regressao <- lm(y[,2] ~ y[,3])
summary(EqCM_regressao)
EqCM_res=ts(EqCM_regressao$residuals)

#'
#' plot of EqCM_res
#'
par(mfrow=c(1,1))
plot(EqCM_res,type='l', col='steelblue', main = 'Equilibrium correction using beta from OLS')


#'
#' ACF and PACFfor EqCm_res
#'
par(mfrow=c(2,1))
Acf(EqCM_res, lag.max=12)
Pacf(EqCM_res, lag.max=12)


#'
#' test unit roots in EqCM_res
#'
EqCM_res_adf = adf.test(EqCM_res)
EqCM_res_adf

y1 <- y[2:201,1]
y2 <- y[2:201,2]
y3 <- y[2:201,3]
#'
#' lag one period transforme y2 and y3 
#'
y1_1=ts(c(NA,y1[1:length(y1)-1]))
y2_1=ts(c(NA,y2[1:length(y2)-1]))
y3_1=ts(c(NA,y3[1:length(y3)-1]))


#'
#' transforme y2 and y3 in first difference
#'
dy1 = y1-y1_1
dy2 = y2-y2_1
dy3 = y3-y3_1

#'
#' lag one period  dy2 and dy3 
#'
dy1_1=ts(c(NA,dy1[1:length(dy1)-1]))
dy2_1=ts(c(NA,dy2[1:length(dy2)-1]))
dy3_1=ts(c(NA,dy3[1:length(dy3)-1]))

#'
#' lag one period EqCM_res
#'
EqCM_res <- EqCM_res[2:201]
EqCM_res_1=ts(c(NA,EqCM_res[2:length(EqCM_res)-1]))

#'
#' ## Equation for $dsp500spot_{t} = \alpha_{1} + \beta_{11}*EqCMres_{t-1} + \beta_{12}* dsp500spot_{t-1} + \beta_{13}*dsp500fut_{t-1}$
#'
eq_dy2 <- lm(dy2 ~ EqCM_res_1 + dy2_1 +  dy3_1)
summary(eq_dy2)


#'
#' ## Equation for $dsp500fut_{t} = \alpha_{2} + \beta_{21}*EqCMres_{t-1} + \beta_{22}* dsp500spot_{t-1} + \beta_{23}*dsp500fut{t-1}$
#'
eq_dy3 <- lm(dy3 ~ EqCM_res_1 + dy2_1 +  dy3_1)
summary(eq_dy3)



#'
#' ## Define the multivariate time series $(spot_{t}, fut_{t})$
#'
sys_y2_y3=ts(cbind(y2,y3))
#sys_y1_y3=ts(cbind(y1,y3))



#'
#' ## VAR lag length selection using the package VARS
#'

VARselect(sys_y2_y3,lag.max=20, type=c("const"))
#VARselect(sys_y1_y3,lag.max=20, type=c("const")) 

#'
#' best model is VAR(2) using AIC, HQ and FPE and VAR(1) using SC
#'

#'
#' we are going to use VAR(1)
#'




VAR_y2_y3_1 <-VAR(sys_y2_y3 , p=1, type = "const")

VAR_y2_y3_1

#'
#' summary of equation spot
#'
summary(VAR_y2_y3_1, equation = "y2")

#'
#' plot fit and residuals and ACF and PACF for residuals for SPOT
#'#'
plot(VAR_y2_y3_1, names = "y2")

#'
#' summary of equation fut
#'
summary(VAR_y2_y3_1, equation = "y3")
#'
#' plot fit and residuals and ACF and PACF for residuals for FUT
#'
plot(VAR_y2_y3_1, names = "y3")


#'
#'  Portmanteau test for system 
#'
ser_y2_y3_1 <- serial.test(VAR_y2_y3_1, lags.pt = 1, type = "PT.asymptotic")
ser_y1_y3_1

#'
#' JB Test for system
#'
JB_y2_y3_1 <- normality.test(VAR_y2_y3_1)
JB_y1_y3_1$jb.mul

#'
#' ARCH test for system
#'
arch_y2_y3_1 <- arch.test(VAR_y2_y3_1, lags.multi = 1)
arch_y2_y3_1$arch.mul


#'
#' diagnostic for spot equation
#'
plot(arch_y2_y3_1,name= "y1")
plot(stability(VAR_y2_y3_1),nc=2)


plot(arch_y2_y3_1,name= "y3")



#'
#' Causality test
#' 

#'
#'Granger causality H0: FUT do not Granger-cause SPOT
#'

causality(VAR_y2_y3_1, cause=c("y2"))

#'
#'Granger causality H0: SPOT  do not Granger-cause FUT
#'

causality(VAR_y2_y3_1, cause=c("y3"))

#'
#' ## Johansen procedure for SPOT and FUT
#'
df_y2_y3 <- as.data.frame(cbind(y[,2],y[,3]))
df_y2_y3 <- df_y2_y3[2:201,]
vec_y2_y3_eigen <- (ca.jo(df_y2_y3, type=c("eigen"), ecdet = c("const"), K=2, spec = c("transitory")))
summary(vec_y2_y3_eigen)

vec_y2_y3_trace <- (ca.jo(df_y2_y3, type=c("trace"), ecdet = c("const"), K=2, spec = c("transitory")))
summary(vec_y2_y3_trace)

vec_z<- ca.jo(df_y2_y3, type=c("eigen"), ecdet = c("const"), K=2, spec = c("longrun"))
summary(vec_z)


#' ## Specifying a VEqCM
#' 
#' We will specify a **restricted** VEqCM with constant inside the EqCM
#' but no deterministic component in the short run dynamic. 
#' Since our VAR model is of order 20, the VEC model
#' will be of order 19, and the long term relation.
#' 
#' This long-term relation has already been estimated by our VAR(20) 
#' regression, thus, the VEqCM specification:
#' 
#' $$\Delta y_t = \alpha\beta^\prime y_{t-1}+\varepsilon_t \quad\text{VEqCM}$$
#' 
#' is equivalent to the VAR(1) specification
#' 
#' $$y_t = \Pi y_{t-1}+\varepsilon_t \quad\text{VAR}$$
#' 
#'  where $\Pi = \mathbf{I}+\alpha\beta^\prime$. Therefore, from our VAR(1) model
#'  already estimated we have all information we need to have the VEqCM 
#'  representation. 
#'  
#'  The estimated matrix $\Pi$ is given by:
#'  
Bcoef(VAR_y2_y3_1)


#'  We can use the function VECM from tsDyn package to estimate the VECM(19) with 
#'  one cointegrating vectors with constant in the EqCM (**LRinclude=c("const")**) 
#'  but not in the VAR (**include=c("none")**)
#'  
#'  


head(df_y2_y3)

colnames(df_y2_y3) <- c("y2", "y3")
vec <- VECM(df_y2_y3, lag = 1, r=1, estim = "ML", include = c("none"), LRinclude = c("const"))
summary(vec)


vec1 <- VECM(df_y2_y3, lag = 1, r=1, estim = "ML", include = c("none"))
summary(vec1)
toLatex(vec1)


#' where  EqCM = $\beta^{'}y_{t-1}$


#'
#' ## Test a restriction in $\alpha$ using alrtest $\alpha_{1}=0$ using the package urca weak exogeneity for spot in fut conditional model
#' 





DA_spot = matrix (c(0,1),c(2,1))
DA_spot
teste_alpha_spot <- alrtest(vec_spot_fut_eigen,A=DA_spot,r=1)
summary(teste_alpha_spot)



#'
#' ## Test a restriction in $\alpha$ using alrtest $\alpha_{2}=0$ using the package urca weak exogeneity for fut in spot conditional model 
#' 



DA_fut = matrix (c(1,0),c(2,1))
DA_fut
teste_alpha_fut <- alrtest(vec_spot_fut_eigen,A=DA_fut,r=1)
summary(teste_alpha_fut)


#'
#' ## Test a restriction in $\beta$ using blrtest $\beta_{1}=1$ and $\beta_{2}=-1$ 
#' using the package urca for an one-to-one in the long-run 
#' relationship between spot and fut 
#'  

HD_prop <- matrix(c(1,-1,0),c(3,1))
HD_prop
test_beta_prop <- blrtest(vec_spot_fut_eigen,H=HD_prop,r=1)
summary(test_beta_prop)

#'
#´## Test both restrictions $\beta_{1}=1$ and $\beta_{2}=-1$  and $\alpha_{2}=0$ 
#'
test_alpha_beta_fut <- ablrtest(vec_spot_fut_eigen,H=HD_prop,A=DA_fut,r=1)
summary(test_alpha_beta_fut)



#'
#´## Test both restrictions $\beta_{1}=1$ and $\beta_{2}=-1$  and $\alpha_{1}= 0$ 
#'
test_alpha_beta_spot <- ablrtest(vec_spot_fut_eigen,H=HD_prop,A=DA_spot,r=1)
summary(test_alpha_beta_spot)


#'
#' ## IRF
#' 


vec2var_vec_z <- vec2var(vec_z,r=1)


irf_vec_z_ort <- irf(vec2var_vec_z,n.ahead=10,ortho=TRUE)

plot(irf_vec_z_ort)

irf_vec_z_cumul <- irf(vec2var_vec_z,n.ahead=10,ortho=FALSE, cumulative = TRUE)

plot(irf_vec_z_cumul)

irf_vec_z <- irf(vec2var_vec_z,n.ahead=10,ortho=FALSE, cumulative = FALSE)

plot(irf_vec_z)


#'
#' ## Variance Decomposition
#' 

fevd(vec2var_vec_z, n.ahead=10)

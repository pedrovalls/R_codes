#' ---
#' title: "Gaussian White Noise"
#' author: "Pedro Valls"
#' date: "21/06/2021"
#' ---

# clear workspace
rm(list=ls())
# set local directory to where data set is 
setwd("C:/Users/Pedro/Dropbox/EcoIII2021/Livro/R_CODES/White_Noise")


# Load package using a function load_package-----------------------------------------------------------------
load_package<-function(x){
  x<-as.character(match.call()[[2]])
  if (!require(x,character.only=TRUE)){
    install.packages(pkgs=x,repos="http://cran.r-project.org")
    require(x,character.only=TRUE)
  }
}

load_package('tseries')
load_package('stats')
load_package('graphics')
load_package('ggplot2')


# Load Library
library(data.table)
library(ggplot2)
library(tseries)
library(xtable)
library(forecast)
library(FinTS)

#' 
#' # Fix the random seed
#' 
set.seed(123456)


#' 
#' # Generate a normal random variable with 1000 elements
#' 
gwn<-rnorm(1000)
#' 
#' # Transform to time series object 
#' 
gwn.ts=ts(gwn)


#' 
#' # plot the time series gwn 
#' 
plot(gwn.ts,col='blue', type='l', ylab='gwn', main = 'Gaussian White Noise')
#' 
#' # Define series lagged one period for the time series gwn call it gwn.ts_1
#' 
gwn.ts_1=lag(gwn.ts, k=-1)
#'
#' # Define series lagged two period for the time series gwn call it gwn.ts_2
#' 
gwn.ts_2=lag(gwn.ts, k=-2)
#' 
#' # Define series lagged twenty periods for the time series gwn call it gwn.ts_20
#' 
gwn.ts_20=lag(gwn.ts, k=-20)


#'
#' # Scatter plot between the series and its lags, idea of temporal dependence
#' 
par(mfrow=c(1,1))
plot(gwn.ts_1,gwn.ts, xlab=expression(textstyle(gwn[t-1])), ylab=expression(textstyle(gwn[t])), main = 'scatter plot gwn[t-1] X gwn[t]', col = 'blue') # Y and Y(-1)
plot(gwn.ts_2,gwn.ts, xlab=expression(textstyle(gwn[t-2])), ylab=expression(textstyle(gwn[t])), main = 'scatter plot gwn[t-2] X gwn[t]', col = 'red') # Y and Y(-2)
plot(gwn.ts_20,gwn.ts, xlab=expression(textstyle(gwn[t-20])), ylab=expression(textstyle(gwn[t])), main = 'scatter plot gwn[t-20] X gwn[t]', col = 'black') # Y and Y(-20)



#'
#' # Test normality using histogram
#' 
par(mfrow=c(1,1))


#'
#' # plot histogram
#' 
# Organize the data in a data frame structure
gwn.df <- data.frame(gwn=gwn.ts)
ggplot(gwn.df) + 
  geom_histogram(mapping = aes(x=gwn, y =..density..),
                 bins = 50, 
                 colour="black",
                 fill="blue", alpha=0.4) +
  labs(title = "Histogram", x = NULL, y = "Density")

#'
#' # Declare a data frame with points and pdf of a normal distribution with the same mean and standard deviation as the series gwn
#' 

# set the points of the normal distribution to be ploted
normal.points.to.plot <- seq(range(gwn.ts)[1], range(gwn.ts)[2], length.out=1000)

# Table with the normal distrinution points
Normal.df <- data.frame(point = normal.points.to.plot,
                        pdf = dnorm(normal.points.to.plot, mean(gwn.ts), sd(gwn.ts)))

# Clean up of unused variables
rm(list = c("normal.points.to.plot")) 

# plot the histogram of gwn, with the Kernel and a normal distribution
ggplot(gwn.df) + 
  geom_histogram(mapping = aes(x=gwn, y =..density.., fill = "gwn"),
                 bins = 50, 
                 colour="black", alpha=0.4) +
  geom_density(mapping = aes(x=gwn, fill="Kernel", colour=NULL), alpha = 0.3) + 
  geom_line(mapping = aes(x=point, y=pdf, colour="Normal"), data=Normal.df, size =1) +
  labs(title = "Histogram", x = NULL, y = "Density",
       colour = NULL, # legend name of colours
       fill = NULL # legend name of fill
  ) +
  scale_color_manual(breaks = c("Normal"),
                     values=c("red")) +
  scale_fill_manual(breaks = c("gwn", "Kernel"),
                    values=c("blue", "green")) + 
  theme_bw() +
  theme(legend.position="bottom")


#'
#'  # Testing normality using Jarque Bera Test
#'
jarque.bera.test(gwn.ts)


#'
#' # Test for constant mean and variance
#' 
mean_gwn<-mean(gwn.ts[1:100])
stdev_gwn<-sd(gwn.ts[1:100])
plot(gwn.ts, col='blue', type='l')
abline(h=mean_gwn+1.96*stdev_gwn, col='red', lty=2)
abline(h=mean_gwn-1.96*stdev_gwn, col='red', lty=2)


#'
#' # Normal Quantile-Quantile Plot 
#' 
qqnorm(gwn.ts, main = "Normal Q-Q Plot", xlab = "Theoretical Quantiles",
       ylab="Empirical Quantile", plot.it=TRUE, datax=FALSE)


p <- ggplot(gwn.df, aes(sample = gwn))
p + stat_qq() + stat_qq_line(color="red") +
  labs(
    title="Normal Q-Q Plot",
    x="Theoretical Quantiles",
    y="Empirical Quantile",
    colour = "Legenda"
  ) +
  theme_bw()


#'
#' # ACF for GWN
#' 
par(mfrow=c(1,1))
acf(gwn.ts, lag.max = 12)
pacf(gwn.ts, lag.max = 12)


#'
#' # Compute the $Q$ stats
#' 

#'
#' # use the 12 ACF and PACF to compute $Q$ stats
#' 
gwn.acf = acf(gwn.ts, lag.max = 12, plot = FALSE)
gwn.pacf = pacf(gwn.ts, lag.max = 12, plot = FALSE)

#'
#' # length of the series
#' 
T<-length(gwn.ts)

#'
#' # define the firts 5 Autocorrelation
#' 
a1<-gwn.acf$acf[2]
a2<-gwn.acf$acf[3]
a3<-gwn.acf$acf[4]
a4<-gwn.acf$acf[5]
a5<-gwn.acf$acf[6]
#'
#' # compute $Q$ and $Q^{*}$ Stats
#' 
q_gwn<-T*(a1*a1+a2*a2+a3*a3+a4*a4+a5*a5)
qs_gwn<-T*(T+2)*((a1*a1/(T-1)+a2*a2/(T-2)+a3*a3/(T-3)+a4*a4/(T-4)+a5*a5/(T-5)))

#'
#' # Compute the p-value for $Q$ stats
#' 
(pvalue_q_gwn<-pchisq(q_gwn,5, lower.tail = F))

#'
#' # Compute the p-value for $Q^{*}$ stats
#' 
(pvalue_qs_gwn<-pchisq(qs_gwn,5, lower.tail = F))

#'
#' # compute the theoretical population spectrum for a Gaussian white noise
#' 
#'
#' # Define the series $\frac{\gamma_{0}}{2 \pi}$ by multiplying a vector of ones by this quantity
#' 
n = 1000
one_vector = rep(1, n)
spec_gwn.ts=ts(one_vector*((sd(gwn.ts)^2/(2*pi))))
par(mfrow=c(1,1))
plot(spec_gwn.ts, ylab='population spectrum', xlab='', main = 'Population Spectrum for GWN', col='red')

#'
#' # compute the spectral density for the series gwn
#' 


spec.pgram(gwn.ts, kernel("daniell", c(30,50)), taper = 0.1,
           pad = 0, fast = TRUE, demean = FALSE, detrend = TRUE,
           plot = TRUE, na.action = na.fail)





# export gwn.ts to excel
write.csv2(gwn.ts, file='gwn1.csv')


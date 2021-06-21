#' 
#' ## "Generate Stochastic Processes and a realization of Time Series"
#' ## "Pedro Valls"
#' ## "June 20th, 2021"
#'



# Clear workspace

rm(list=ls())

#set local directory to where data set is
setwd("C:/Users/Pedro/Dropbox/EcoIII2021/Livro/R_CODES")

# Load Packages Functions
load_package<-function(x){
  x<-as.character(match.call()[[2]])
  if (!require(x,character.only=TRUE)){
    install.packages(pkgs=x,repos="http://cran.r-project.org")
    require(x,character.only=TRUE)
  }
}
load_package("readxl")
load_package("openxlsx")
load_package("xts")
load_package("fBasics")
load_package("qrmtools")
load_package("zoo")
load_package("graphics")
load_package("tseries")
load_package("stats")
load_package("moments")
load_package("MASS")
load_package("MultiRNG")
load_package("tidyverse")
load_package("cowplot")
load_package("urca")
load_package("vars")
load_package("tsDyn")
load_package("ggplot2")



# Load Library
library(data.table)
library(ggplot2)
library(tseries)
library(lmtest)
library(Metrics)
library(stargazer)
library(xtable)
library(forecast)
library(xtable)
library(FinTS)







#' ## Fix sample sixe
#' 


n <-100

#' # Initialize a matrix $y$ with dimension 100
#' 


y <- matrix(data=0, nrow=100, ncol=100)

#' # Initialize a matrix $z$ with dimension 100
#' 

z <- matrix(data=0, nrow=100, ncol=100)

#' # Initialize a matrix $u$ with dimension 100
#' 

u <- matrix(data=0, nrow=100, ncol=100)



#' # Initialize a matrix $x$ with dimension 100
#' 

x <- matrix(data=0, nrow=100, ncol=100)

#' # Initialize the vectors $x1$, $x2$, $x3$ e $x4$ with the realization of the time series $x$.
#' 
x1 <- matrix(data=0, nrow=100, ncol=1)
x2 <- matrix(data=0, nrow=100, ncol=1)
x3 <- matrix(data=0, nrow=100, ncol=1)
x4 <- matrix(data=0, nrow=100, ncol=1)

#' # Initialize  the $tsx$ vector that will contain the observable time series
#' 
tsx <- matrix(data=0, nrow=100, ncol=1)


#' # Fix the random seed for $y(t,w)$
#' 
set.seed(123456) 


#' # Generate a multivariade normal distribuition with mean = (0,...,0), $\Sigma$ = diag(1.0,...,1.0) 
#' # using mvrnorm para $y(t,w)$ with order 100
#' 
Sigmay <- diag(rep(1,100))
y <- mvrnorm(n, rep(0,100), Sigmay, tol = 1e-6, empirical = TRUE)

#' # Fix random seed for $z(t,w)$ 
#' 

set.seed(654321)

#' # Generate a multivariade normal distribuition with mean = (10,...,10), $\Sigma$ = diag(5.0,...,5.0) 
#' # using mvrnorm para $z(t,w)$ with order 100
#' 

Sigmaz <- diag(rep(5,100))
z <- mvrnorm(n, rep(10,100), Sigmaz, tol = 1e-6, empirical = TRUE)


#' # Fix random seed for $u(t,w)$ 
#' 

set.seed(123)

#' # Generate the process $u(t,w)$ using uniform between [$-\pi$, $\pi$] 
#' 



for(i in 1:n){
  for(j in 1:n){
    u1 <- runif(1, min=-pi,max=pi)
    u[i,j] =u1
  }
} 

#' # Generate $x(t,w)=y(t,w)*cos(z(t,w)+t*u(t,w))$
#' 

for(i  in  1:n){
  for(j in 1:n){
    x[i,j]=y[i,j]*cos(z[i,j]+i*u[i,j])
  }} 

#' # Select the first four columns are the first four realizations of the stochastic process $x(t,w)$ 
#' 

x1 <- ts(x[,1]) 
x2 <- ts(x[,2])
x3 <- ts(x[,3])
x4 <- ts(x[,4])

#' # Transform to data table
#' 

stoch.plot <- data.table('id'=1:100, 'x1'= x1, 'x2'=x2, 'x3'=x3, 'x4'=x4)

#' # Plot the four realizations of the time series
#' 

ggplot(stoch.plot)+
  geom_line(aes(x=id, y=x1,colour="x1"),linetype="solid") +
  geom_line(aes(x=id, y=x2,colour="x2"),linetype="solid") +
  geom_line(aes(x=id, y=x3,colour="x3"),linetype="solid") +
   geom_line(aes(x=id, y=x4,colour="x4"),linetype="solid") +
  scale_color_manual(values=c("x1"="#FF0000",
                                "x2"="blue",
                              "x3"="black",
                              "x4"="green"
                              ))+
  labs(
  title="Four realizations of the time series x(t,w)",
  x=NULL,
  y=NULL,
  colour = "Legenda"
  ) +
  theme_bw()
                              
#' # Generate the $tsx$ as $1^{rt}$ obs $x1$, $2^{nd}$ obs $x2$, $3^{rd}$ obs $x3$, $4^{th}$ obs $x4$,...
#' 

teste <- 0
teste1 <- 0
tsx <- 0



for(i  in  1:n){
  teste[i] = floor((i)/4)
  teste1[i]=i-(teste[i])*4
  
  if(teste1[i] == 1){
    tsx[i] <- x1[i] } 
  else if(teste1[i] == 2){
    tsx[i] <- x2[i] } 
  else if(teste1[i] == 3){
    tsx[i] <- x3[i]} 
  else if(teste1[i] == 0){
    tsx[i] <- x4[i]} 
  else {
   
  }
}


tsx = ts(tsx)


#' # Plot of the observable the time series
#' 

ggplot(stoch.plot)+
  geom_line(aes(x=id, y=tsx,colour="tsx"),linetype="solid") +
  scale_color_manual(values=c("tsx"="red"
                              ))+
  labs(
    title="Observable time series x(t,w)",
    x=NULL,
    y=NULL,
    colour = "Legenda"
  ) +
  theme_bw()
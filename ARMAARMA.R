#' 
#' # Simulac?o de AR(1), MA(1) e ARMA(1,1)
#' 
#' 
#' 
#'  # Pedro Valls
#'  
#'  
#'  # 20/07/2021
#'  

# clear workspace
  rm(list=ls())
# set local directory to where data set is 
# 
#setwd("C:/Users/pedrovallspereira/Dropbox/Econometria_Financeira_Morettin/Rcodes")
#
# Load package using a function load_package-----------------------------------------------------------------
#
load_package<-function(x){
  x<-as.character(match.call()[[2]])
  if (!require(x,character.only=TRUE)){
    install.packages(pkgs=x,repos="http://cran.r-project.org")
    require(x,character.only=TRUE)
  }
}

load_package('tseries')
load_package('ggplot2')
load_package('stats')
load_package('xtable')
load_package('strucchange') 
load_package('xts')
load_package('forecast')
load_package('cents')
library("strucchange")

#
# Fixa a semente aleat?ria
#
set.seed(12345)
#
# Gera uma N(0,1) com T=100
#
gwn<-rnorm(100)

# 
# Transforma em um objeto ts (time series)
#
epsilon=ts(gwn)


#
# Gera um AR(1) com psi=0.8 e constante igual a zero
#

ar1<- epsilon[1]
for(i in 2:100){
  ar1[i]= 0.8*ar1[i-1]+ epsilon[i] 
}
par(mfrow=c(1,1))
plot(ar1, type='l', xlabel = 't', ylabel='X(t)', col='blue', main = '')
ar1_ts = ts(ar1)

#
# F.A.C. e F.A.C.P. para ar1
#

par(mfrow=c(1,2))
acf(ar1_ts, lag.max=20,main='AR(1)')
pacf(ar1_ts, lag.max=20,main='AR(1)')

#
# Ajustar um AR aos dados por MV aproximada
#
ar1_fit <- arma(ar1_ts[2:100], order = c(1,0), include.intercept = FALSE)
summary(ar1_fit)

#
# Ajusta um AR aos dados por MV exata
#
ar1_fit_exato <- arima(ar1_ts, order=c(1,0,0), include.mean = FALSE, method= "ML", optim.method="BFGS")
ar1_fit_exato


# 
# Gera um MA(1) com o mesmo epsilon que o AR e theta =0,8
#

ma1<- epsilon[1]
for(i in 2:100){
  ma1[i]= epsilon[i]- 0.8*epsilon[i-1] 
}
par(mfrow=c(1,1))
plot(ma1, type='l', xlabel = 't', ylabel='X(t)', col='blue', main = '')
ma1_ts = ts(ma1)

#
# F.A.C. e F.A.C.P. para ar1
#

par(mfrow=c(1,2))
acf(ma1_ts, lag.max=20,main='MA(1)')
pacf(ma1_ts, lag.max=20,main='MA(1)')

#
# Ajustar um MA aos dados por MV aproximada
#
ma1_fit <- arma(ma1_ts, order = c(0,1), include.intercept = FALSE)
summary(ma1_fit)

#
# Ajusta um MA aos dados por MV exata
#
ma1_fit_exato <- arima(ma1_ts, order=c(0,0,1), include.mean = FALSE, method= "ML", optim.method="BFGS")
ma1_fit_exato



# 
# Gera um ARMA(1,1) com o mesmo epsilon psi=0,8 e theta = 0,3
#

arma11<- epsilon[1]
for(i in 2:100){
  arma11[i]=  0.8*arma11[i-1] + epsilon[i]- 0.3*epsilon[i-1] 
}
par(mfrow=c(1,1))
plot(arma11, type='l', xlabel = 't', ylabel='X(t)', col='blue', main = '')
arma11_ts = ts(arma11)

#
# F.A.C. e F.A.C.P. para ar1
#

par(mfrow=c(1,2))
acf(arma11_ts, lag.max=20,main='ARMA(1,1)')
pacf(arma11_ts, lag.max=20,main='ARMA(1,1)')

#
# Ajustar um ARMA aos dados por MV aproximada
#
arma11_fit <- arma(arma11_ts, order = c(1,1), include.intercept = FALSE)
summary(arma11_fit)

#
# Ajusta um ARMA aos dados por MV exata
#
arma11_fit_exato <- arima(arma11_ts, order=c(1,0,1), include.mean = FALSE, method= "ML", optim.method="BFGS")
arma11_fit_exato

#
# Exemplo 3.8 tabela a
#

n_ar1 = length(ar1)

tabAR1=matrix(0, nrow = 3, ncol =2)
#
# coef
#
tabAR1[1,1] = ar1_fit$coef
tabAR1[1,2]= ar1_fit_exato$coef
#
#erro padrao
#
tabAR1[2,1]=sqrt(ar1_fit$vcov)
tabAR1[2,2]=sqrt(ar1_fit_exato$var.coef)

#
# p-valor
#
t.value_ar1= tabAR1[1,1]/tabAR1[2,1]
t.value_ar1_exato = tabAR1[1,2]/tabAR1[2,2]
p.valor_ar1_exato=pt(t.value_ar1_exato,n_ar1-1, lower.tail=FALSE)
p.valor_ar1=pt(t.value_ar1,n_ar1-1, lower.tail=FALSE)
if(p.valor_ar1 < .00001) {p.valor_ar1=0.0000}

tex_tabAR1 = xtable(tabAR1,digits=4)
print(tex_tabAR1)




#
# Exemplo 3.8 tabela b
#

n_ma1 = length(ma1)

tabMA1=matrix(0, nrow = 3, ncol =2)
#
# coef
#
tabMA1[1,1] = (-1)*ma1_fit$coef
tabMA1[1,2]= (-1)*ma1_fit_exato$coef
#
#erro padrao
#
tabMA1[2,1]=sqrt(ma1_fit$vcov)
tabMA1[2,2]=sqrt(ma1_fit_exato$var.coef)

#
# p-valor
#
t.value_ma1= tabMA1[1,1]/tabMA1[2,1]
t.value_ma1_exato = tabMA1[1,2]/tabMA1[2,2]
p.valor_ma1_exato=pt(t.value_ma1_exato,n_ma1-1, lower.tail=FALSE)
p.valor_ma1=pt(t.value_ma1,n_ma1-1, lower.tail=FALSE)
if(p.valor_ma1 < .00001) {p.valor_ma1=0.0000}

tex_tabMA1 = xtable(tabMA1,digits=4)
print(tex_tabMA1)


#
# Exemplo 3.8 tabela c
#

n_arma11 = length(arma11)

tabARMA11=matrix(0, nrow = 2, ncol =6)
#
# coef ao AR
#
tabARMA11[1,1] = arma11_fit$coef[1]
tabARMA11[1,4]= arma11_fit_exato$coef[1]

#
# coef ao MA
#
tabARMA11[2,1] = (-1)*arma11_fit$coef[2]
tabARMA11[2,4]= (-1)*arma11_fit_exato$coef[2]
#
#erro padrao do AR
#
tabARMA11[1,2]=sqrt(arma11_fit$vcov[1,1])
tabARMA11[1,5]=sqrt(arma11_fit_exato$var.coef[1,1])

#
#erro padrao do MA
#
tabARMA11[2,2]=sqrt(arma11_fit$vcov[2,2])
tabARMA11[2,5]=sqrt(arma11_fit_exato$var.coef[2,2])



#
# p-valor para AR
#
t.value_ARMA11_AR= tabARMA11[1,1]/tabARMA11[1,2]
t.value_ARMA11_AR_exato = tabARMA11[1,4]/tabARMA11[1,5]
p.valor_ARMA11_AR_exato=pt(t.value_ARMA11_AR_exato,n_arma11-2, lower.tail=FALSE)
p.valor_ARMA11_AR=pt(t.value_ARMA11_AR,n_arma11-2, lower.tail=FALSE)
if(p.valor_ARMA11_AR < .00001) {p.valor_ARMA11_AR=0.0000}
if(p.valor_ARMA11_AR_exato < .00001) {p.valor_ARMA11_AR_exato=0.0000}
tabARMA11[1,3]=p.valor_ARMA11_AR
tabARMA11[1,6]=p.valor_ARMA11_AR_exato

#
# p-valor para MA
#
t.value_ARMA11_MA= tabARMA11[2,1]/tabARMA11[2,2]
t.value_ARMA11_MA_exato = tabARMA11[2,4]/tabARMA11[2,5]
p.valor_ARMA11_MA_exato=pt(t.value_ARMA11_MA_exato,n_arma11-2, lower.tail=FALSE)
p.valor_ARMA11_MA=pt(t.value_ARMA11_MA,n_arma11-2, lower.tail=FALSE)
if(p.valor_ARMA11_MA < .00001) {p.valor_ARMA11_MA=0.0000}
if(p.valor_ARMA11_MA_exato < .00001) {p.valor_ARMA11_MA_exato=0.0000}
tabARMA11[2,3]=p.valor_ARMA11_MA
tabARMA11[2,6]=p.valor_ARMA11_MA_exato



tex_tabARMA11 = xtable(tabARMA11,digits=4)
print(tex_tabARMA11)


#
# define a serie defasada de ar1 
#
ar1_ts_1=c(NA,ar1_ts[1:length(ar1_ts)-1])
#
# Estimacao do AR(1) por MQO
#
eq1 <- lm(ar1_ts ~ ar1_ts_1)
summary(eq1)

#
# Residuos Recursivos podem ser obtidos diretamente do pacote strucchange
#
par(mfrow = c(1,1))

residrls <- strucchange::recresid(eq1)
plot(residrls, type = "line", main = "Recursive Residuals")
ocus <- efp(eq1, type="Rec-CUSUM", data=ar1_ts)

bound.ocus <- boundary(ocus, alpha=0.05)
plot(ocus)
sctest(ocus)


## compute squared recursive residuals
w2 <- recresid(eq1, data = ar1_ts)^2
## compute CUSUM of squares process
sr <- ts(cumsum(c(0, w2))/sum(w2), end = end(ar1_ts), freq = 12)
## the border (r-k)/(T-k)
border <- ts(seq(0, 1, length = length(sr)),
             start = start(sr), freq = 12)

## nice plot
plot(sr, xaxs = "i", yaxs = "i", main = "CUSUM of Squares")
lines(border, col = grey(0.5))
lines(0.4 + border, col = grey(0.5))
lines(- 0.4 + border, col = grey(0.5))



#
# MQR forca bruta
#
beta_ar1=0
sd_beta_ar1=0
li_ar1=0
ls_ar1=0

for (i in  6:100){
  ar1_fit_rec <- arma(ar1_ts[2:i], order = c(1,0), include.intercept = FALSE)
  beta_ar1[i]=ar1_fit_rec$coef[1]
  sd_beta_ar1[i]=sqrt(ar1_fit_rec$vcov[1,1])
  li_ar1[i]=beta_ar1[i]-1.96*sd_beta_ar1[i]
  ls_ar1[i]=beta_ar1[i]+1.96*sd_beta_ar1[i]
  }

dat_ar <- matrix(0, nrow=95, ncol=4)
for (i in 1:95){
  dat_ar[i,1]=0
  dat_ar[i,2]=li_ar1[i+5]
  dat_ar[i,3]=beta_ar1[i+5]
  dat_ar[i,4]=ls_ar1[i+5]
  
}

par(mfrow = c(1,1))

matplot(dat_ar, type = "l",pch=1,col = 1:4, main="Estimativa Recursiva do Coef AR(1)") #plot
legend("bottom", legend = c("zero","li","coef","ls"), col=1:4, pch=15,trace=TRUE,horiz=TRUE) # optional legend


#
# MQR forca bruta para MA(1)
#
theta_const=0
theta_ma1=0
sd_theta_const=0
sd_theta_ma1=0
li_ma1=0
li_ma1_const=0
ls_ma1=0
ls_ma1_const=0

for (i in  6:100){
  ma1_fit_rec <- arima(ma1_ts[2:i], order=c(0,0,1), include.mean = FALSE, method= "ML", optim.method="BFGS")
#  theta_const[i]=ma1_fit_rec$coef[1]
  theta_ma1[i]=ma1_fit_rec$coef[1]
  sd_theta_ma1[i]=sqrt(ma1_fit_rec$var.coef[1,1])
 # sd_theta_const[i]=sqrt(ma1_fit_rec$var.coef[1,1])
  li_ma1[i]=theta_ma1[i]-1.96*sd_theta_ma1[i]
  ls_ma1[i]=theta_ma1[i]+1.96*sd_theta_ma1[i]
#  li_ma1_const[i]=theta_const[i]-1.96*sd_theta_const[i]
 # ls_ma1_const[i]=theta_const[i]+1.96*sd_theta_const[i]
  
}

dat_ma <- matrix(0, nrow=95, ncol=4)
for (i in 1:95){
  dat_ma[i,1]=0
  dat_ma[i,2]=li_ma1[i+5]
  dat_ma[i,3]=theta_ma1[i+5]
  dat_ma[i,4]=ls_ma1[i+5]
  
}
#dat_const <- matrix(0, nrow=95, ncol=4)
#for (i in 1:95){
  #dat_const[i,1]=0
  #dat_const[i,2]=li_ma1_const[i+5]
  #dat_const[i,3]=theta_const[i+5]
  #dat_const[i,4]=ls_ma1_const[i+5]
  
#}


par(mfrow = c(1,1))

matplot(dat_ma, type = "l",pch=1,col = 1:4, main="Estimativa Recursiva do Coef MA(1)") #plot
legend("bottom", legend = c("zero","li","coef","ls"), col=1:4, pch=15,trace=TRUE,horiz=TRUE) # optional legend

#matplot(dat_const, type = "l",pch=1,col = 1:4, main="Estimativa Recursiva do Coef da Constante") #plot
#legend("bottom", legend = c("zero","li","coef","ls"), col=1:3, pch=15,trace=TRUE,horiz=TRUE) # optional legend







#
# MQR forca bruta ARMA(1,1)
#
theta_arma11=0
beta_arma11=0
sd_theta_arma11=0
sd_beta_arma11=0
li_theta_arma11=0
li_beta_arma11=0
ls_theta_arma11=0
ls_beta_arma11=0


for (i in  6:100){
  
  arma11_fit_rec <- arima(arma11_ts[2:i], order=c(1,0,1), include.mean = FALSE, method= "ML", optim.method="BFGS")
  beta_arma11[i]=arma11_fit_rec$coef[1]
  theta_arma11[i]=arma11_fit_rec$coef[2]
  sd_beta_arma11[i]=sqrt(arma11_fit_rec$var.coef[1,1])
  sd_theta_arma11[i]=sqrt(arma11_fit_rec$var.coef[2,2])
  li_beta_arma11[i]=beta_arma11[i]-1.96*sd_beta_arma11[i]
  ls_beta_arma11[i]=beta_arma11[i]+1.96*sd_beta_arma11[i]
  li_theta_arma11[i]=theta_arma11[i]-1.96*sd_theta_arma11[i]
  ls_theta_arma11[i]=theta_arma11[i]+1.96*sd_theta_arma11[i]
  
}

dat_arma11_ar<- matrix(0, nrow=95, ncol=4)
for (i in 1:95){
  dat_arma11_ar[i,1]=0
  dat_arma11_ar[i,2]=li_beta_arma11[i+5]
  dat_arma11_ar[i,3]=beta_arma11[i+5]
  dat_arma11_ar[i,4]=ls_beta_arma11[i+5]
}

dat_arma11_ma <- matrix(0, nrow=95, ncol=4)
for (i in 1:95){
  dat_arma11_ma[i,1]=0
  dat_arma11_ma[i,2]=li_theta_arma11[i+5]
  dat_arma11_ma[i,3]=theta_arma11[i+5]
  dat_arma11_ma[i,4]=ls_theta_arma11[i+5]
  
}


par(mfrow = c(1,1))

matplot(dat_arma11_ar, type = "l",pch=1,col = 1:4, main="Estimativa Recursiva do Coef AR no ARMA(1,1)") #plot
legend("bottom", legend = c("zero","li","coef","ls"), col=1:4, pch=15,trace=TRUE,horiz=TRUE) # optional legend

matplot(dat_arma11_ma, type = "l",pch=1,col = 1:4, main="Estimativa Recursiva do Coef MA no ARMA(1,1)") #plot
legend("bottom", legend = c("zero","li","coef","ls"), col=1:4, pch=15,trace=TRUE,horiz=TRUE) # optional legend

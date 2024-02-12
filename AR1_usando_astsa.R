library(astsa)
library(forecast)
# Generate Data
set.seed(123456)
num = 1000
N = num+1
x = sarima.sim(n=N, ar=.75)
y = ts(x[-1])   

par(mfrow=c(1,2))
Acf(y)
Pacf(y)


# Initial Estimates 
u = ts.intersect(y, lag(y,-1), lag(y,-2)) 
varu = var(u)
coru = cor(u) 
phi = coru[1,3]/coru[1,2] 
q = (1-phi^2)*varu[1,2]/phi 
r = 0
(init.par = c(phi, sqrt(q), sqrt(r))) 

# Function to evaluate the likelihood 
Linn=function(para){
  phi = para[1]; sigw = para[2]; sigv = para[3] 
  Sigma0 = (sigw^2)/(1-phi^2); Sigma0[Sigma0<0]=0 
  kf = Kfilter(y,A=1,mu0=0,Sigma0,phi,sigw,sigv)
  return(kf$like)   
}


# Estimation
(est = optim(init.par, Linn, gr=NULL, method="BFGS", hessian=TRUE, control=list(trace=1,REPORT=1)))      
SE = sqrt(diag(solve(est$hessian)))
cbind(estimate=c(phi=est$par[1],sigw=est$par[2],sigv=est$par[3]), SE)
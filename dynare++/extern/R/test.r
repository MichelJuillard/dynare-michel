source("dynareR.r")

parameters <- c("beta","gamma","rho","alpha","delta")
varendo <- c("k","c","a")
varexo <- "eps"
parval <- c(.99,2,.9,.3,.025)
vcovmatrix <- matrix(1,1,1)
initval <- c(0.066, 0.43, 0.01)

modeleq <- c("(c/c(1))^gamma*beta*(alpha*exp(a(1))*k^(alpha-1)+1-delta)=1",
             "a=rho*a(-1)+eps",
             "k+c=exp(a)*k(-1)^alpha+(1-delta)*k(-1)")


dd <- calldynare(modeleq,varendo,varexo,parameters,2,parval,vcovmatrix,initval)

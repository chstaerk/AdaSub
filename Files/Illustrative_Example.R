
### Illustrative example of Section A2 of the supplement ############################

# Load main functions (Please adjust to the right directory!)
# setwd("C://...") 
# source("AdaSub_main_functions.R")

simdata.global.corr <- function (n,p,beta,sigma.normal,corr=0) {
  
  x=matrix(numeric(n*p),n,p) #initialize design matrix
  y=numeric(n) #initialize response variable
  
  mu=rep(0,p)
  Sigma=diag(1-corr,p)+matrix(corr,p,p)
  x = mvrnorm(n , mu, Sigma)
  
  linpred=x%*%beta
  y=rnorm(n,linpred,sigma.normal)
  
  return(list(x=x,y=y))
}

p=1000 # number of variables 
n=60 # sample size 
corr=0 # correlation between variables

S0=c(1,2,3,4,5) # true active set 
s0=length(S0)
beta1=numeric(p)
beta1[c(1,2,3,4,5)]=c(0.4,0.8,1.2,1.6,2) # true vector of regression coefficients

set.seed(2)
data = simdata.global.corr(n,p,beta=beta1,sigma.normal=1,corr= corr) # simulate data

Iterations=10000 # number of iterations in AdaSub 
K=n # learning rate in AdaSub 
q=10 #  
savings=1

# First run 
const=0.6 # constant gamma in EBIC
start.time <- Sys.time()
output06=AdaSub(data,Iter=Iterations,K=K,q=q,criterion="EBIC",savings=savings,const=const,U_C=40)
end.time <- Sys.time()
time06 <- end.time - start.time


# Second run 
const=1 # constant gamma in EBIC
start.time <- Sys.time()
output1=AdaSub(data,Iter=Iterations,K=K,q=q,criterion="EBIC",savings=savings,const=const,U_C=40)
end.time <- Sys.time()
time1 <- end.time - start.time

print(paste("Computation time of AdaSub for gamma=0.6 with", Iterations, "iterations: "))
print(time06)
print("#####################################")


print(paste("Computation time of AdaSub for gamma=1 with", Iterations, "iterations: "))
print(time1)
print("#####################################")


print("True model:")
print(S0)
print("#####################################")

print("Thresholded model of AdaSub for EBIC with gamma=0.6:")
print(which(output06$relfreq.final>0.9))
print("#####################################")

print("'Best' model found by AdaSub for EBIC with gamma=0.6:")
print(output06$best.S)
print("#####################################")

print("Thresholded model of AdaSub for EBIC with gamma=1:")
print(which(output1$relfreq.final>0.9))
print("#####################################")

print("'Best' model found by AdaSub for EBIC with gamma=1:")
print(output1$best.S)

#########Plots ##########################################################

font = 1
par(las=1,lwd=1,cex.main=font, cex.lab=font, cex.axis=font)
options(scipen=5)

dev.off()
selected06 = which(output06$relfreq.final>0.9)
selected1 = which(output1$relfreq.final>0.9)
win.graph()
par(mfrow=c(1,2))
plot(output06$values,pch=20,xlab="Iteration",ylab="EBIC value",cex.lab=1.5, main ="EBIC with gamma=0.6")
abline(h=EBIC(data,selected06,const=0.6),col="red")
plot(output1$values,pch=20,xlab="Iteration",ylab="EBIC value",cex.lab=1.5, main ="EBIC with gamma=1")
abline(h=EBIC(data,selected1,const=1),col="red")


#### Evolution of r_j^(t) (only every "savings" iteration)
win.graph()
par(mfrow=c(3,3))
zehner=0
for (i in 1:6)
  plot(output06$relfreq.hist[10*zehner+i,],pch=20,xlab="Iteration",cex.lab=1.5,cex.main=2,ylab=paste("Var",10*zehner+i),ylim=c(0,1))
plot(output06$relfreq.hist[519,],pch=20,xlab="Iteration",cex.lab=1.5,cex.main=2,ylab="Var 519",ylim=c(0,1))
plot(output06$relfreq.hist[731,],pch=20,xlab="Iteration",cex.lab=1.5,cex.main=2,ylab="Var 731",ylim=c(0,1))
plot(output06$relfreq.hist[950,],pch=20,xlab="Iteration",cex.lab=1.5,cex.main=2,ylab="Var 950",ylim=c(0,1))
mtext("Evolution of selection probabilities for EBIC with gamma=0.6",side=3,outer=TRUE,padj=3)

win.graph()
par(mfrow=c(3,3))
zehner=0
for (i in 1:6)
  plot(output1$relfreq.hist[10*zehner+i,],pch=20,xlab="Iteration",cex.lab=1.5,cex.main=2,ylab=paste("Var",10*zehner+i),ylim=c(0,1))
plot(output1$relfreq.hist[519,],pch=20,xlab="Iteration",cex.lab=1.5,cex.main=2,ylab="Var 519",ylim=c(0,1))
plot(output1$relfreq.hist[731,],pch=20,xlab="Iteration",cex.lab=1.5,cex.main=2,ylab="Var 731",ylim=c(0,1))
plot(output1$relfreq.hist[950,],pch=20,xlab="Iteration",cex.lab=1.5,cex.main=2,ylab="Var 950",ylim=c(0,1))
mtext("Evolution of selection probabilities for EBIC with gamma=1",side=3,outer=TRUE,padj=3)


### Plot of sampled sizes of V^(t) and f(V^(t))
win.graph()
par(mfrow=c(1,2))
grays=gray.colors(10)[c(1,4,8)]
plot(output06$V.size,pch=20,ylim=c(0,max(output06$V.size)),xlab="Iteration",ylab="Model size",cex.lab=1.5,col=grays[2], main ="EBIC with gamma=0.6")
points(output06$S.size,col='red',pch=4,cex=1.4)
plot(output1$V.size,pch=20,ylim=c(0,max(output1$V.size)),xlab="Iteration",ylab="Model size",cex.lab=1.5,col=grays[2], main ="EBIC with gamma=1")
points(output1$S.size,col='red',pch=4,cex=1.4)

# Plot of expected size of sampled V^(t)
win.graph()
par(mfrow=c(1,2))
plot(colSums(output06$relfreq.hist),pch=20, ylim=c(0,13), xlab="Iteration", ylab = "Expected Size of V^(t)", cex.lab=1, cex.main=1, main ="EBIC with gamma=0.6")
plot(colSums(output1$relfreq.hist),pch=20,  ylim=c(0,13), xlab="Iteration", ylab = "Expected Size of V^(t)", cex.lab=1, cex.main=1, main ="EBIC with gamma=1")














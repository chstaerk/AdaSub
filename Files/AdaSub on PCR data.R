
### PCR real data example (see Section A4.1 of the supplement) ############################

# Load main functions (Please adjust to the right directory!)
# setwd("C://...") 
# source("AdaSub_main_functions.R")


# leave-one-out cross-validation
LOOCV = function (data,indices){
  n=nrow(data$x)
  data.adj = list()
  data.adj$x = data$x[,indices]
  data.adj$x=cbind(c(rep(1,n)),data.adj$x)                 #include intercept!
  data.adj$y = data$y
  data.cur = list()
  y.hat = numeric(n)
  for (i in 1:n){
    data.cur$x = data.adj$x[-i,]
    data.cur$y = data.adj$y[-i]
    beta.hat.cur =  ginv(t(data.cur$x)%*%data.cur$x) %*% t(data.cur$x) %*% data.cur$y  # use pseudoinverse
    y.hat[i] = data.adj$x[i,] %*% beta.hat.cur 
  }
  mean_CV = mean((y.hat-data$y)^2)
  median_CV = median((y.hat-data$y)^2)
  return(list(mean_CV=mean_CV, median_CV = median_CV))
}

###############################################################
# Please download the PCR data from JRSSB Datasets Vol. 77(5), Song and Liang (2015)
# Link: https://rss.onlinelibrary.wiley.com/hub/journal/14679868/series-b-datasets/pre_2016
# Load data (Please adjust to the right directory!)
#X = t( read.table(".\\data\\PCR\\Xgenes.txt") )
#Y = scan(".\\data\\PCR\\Y3.txt")
#Xnames = scan(".\\data\\PCR\\gene_id.txt",what=character())
###############################################################

data=list()
data$x = X
colnames(data$x) = Xnames
data$y = Y

n = dim(data$x)[1]
p = dim(data$x)[2]



choices <-c( "Number of iterations as in paper (long computation time)", "Custom number of iterations")
mychoice <- menu( choices , graphics=TRUE, title="How many iterations of AdaSub do you want to run for each scenario?" )

if (mychoice==2) {
  Iterations = -1
  while(Iterations < 1 ){
    Iterations <- readline("Enter the number of iterations: ")
    Iterations <- ifelse(grepl("\\D",Iterations),-1,as.integer(Iterations))
    if(is.na(Iterations)){break}  # breaks when hit enter
  }
}


# Specify parameters for AdaSub algorithm

if (mychoice==1)  {
  Iterations= 500000
  Iter_vec = rep(1000,500)
} else {
  if (Iterations%%1000==0) Iter_vec = rep(1000,Iterations%/%1000) else
    Iter_vec = c(rep(1000,Iterations%/%1000),Iterations%%1000)
}
K=n
savings = 100
q=5

const = 0.6 # constant gamma in EBIC


set.seed(123)
start.time <- Sys.time()
#output06=AdaSub(data,Iter=Iterations,K=K,q=q,criterion="EBIC",const=const,savings=savings,plotting=FALSE)
wrapper_output06 = AdaSubWrapper(data,Iter_vec,K=K,q=q,criterion="EBIC",const=const,savings=savings)
output06 = wrapper_output06$output
end.time <- Sys.time()
print(paste("Computation time of AdaSub for gamma=0.6 with", Iterations, "iterations: "))
print(end.time - start.time)
print("#####################################")

model_6 = sort(order(output06$relfreq.final,decreasing=TRUE)[1:6])

if (mychoice==1) {
  Iterations= 50000
  Iter_vec = rep(1000,50)
}

const = 1 # constant gamma in EBIC

set.seed(123)
start.time <- Sys.time()
#output1=AdaSub(data,Iter=Iterations,K=K,q=q,criterion="EBIC",const=const,savings=savings,plotting=FALSE)
wrapper_output1 = AdaSubWrapper(data,Iter_vec,K=K,q=q,criterion="EBIC",const=const,savings=savings)
output1 = wrapper_output1$output
end.time <- Sys.time()
print(paste("Computation time of AdaSub for gamma=1 with", Iterations, "iterations: "))
print(end.time - start.time)
print("#####################################")



## Stability Selection with Lasso
#res = stabpath(data$y,data$x,mc.cores=2)
#error=1
#stability_nonzero = stabsel(res,error=error,type="pfer")$stable
#names(stability_nonzero)=NULL


split_merge_best = which(Xnames =="1429089_s_at" | Xnames =="1430779_at" | Xnames =="1432745_at" | Xnames =="1437871_at" | Xnames =="1440699_at" | Xnames =="1459563_x_at")


print("Thresholded model of AdaSub for gamma=0.6:")
cur.model =which(output06$relfreq.final>0.9)
print(Xnames[cur.model])
print(paste("LOOCV mean error: ", round(LOOCV(data, cur.model)[[1]],3), "    LOOCV median error:", round(LOOCV(data, cur.model)[[2]],3)))
print("#####################################")

print("'Best' model found by AdaSub for gamma=0.6:")
cur.model = wrapper_output06$best.S
print(Xnames[cur.model])
print(paste("LOOCV mean error: ", round(LOOCV(data, cur.model)[[1]],3), "    LOOCV median error:", round(LOOCV(data, cur.model)[[2]],3)))
print("#####################################")

print("Model including 6 variables with highest selection probability in AdaSub for gamma=0.6:")
cur.model = model_6
print(Xnames[cur.model])
print(paste("LOOCV mean error: ", round(LOOCV(data, cur.model)[[1]],3), "    LOOCV median error:", round(LOOCV(data, cur.model)[[2]],3)))
print("#####################################")

print("Thresholded model of AdaSub for gamma=1:")
cur.model = which(output1$relfreq.final>0.9)
print(Xnames[cur.model])
print(paste("LOOCV mean error: ", round(LOOCV(data, cur.model)[[1]],3), "    LOOCV median error:", round(LOOCV(data, cur.model)[[2]],3)))
print("#####################################")

print("'Best' model found by AdaSub for gamma=1:")
cur.model = wrapper_output1$best.S
print(Xnames[cur.model])
print(paste("LOOCV mean error: ", round(LOOCV(data, cur.model)[[1]],3), "    LOOCV median error:", round(LOOCV(data, cur.model)[[2]],3)))
print("#####################################")

print("SAM model with 6 variables:")
cur.model = split_merge_best
print(Xnames[cur.model])
print(paste("LOOCV mean error: ", round(LOOCV(data, cur.model)[[1]],3), "    LOOCV median error:", round(LOOCV(data, cur.model)[[2]],3)))
print("#####################################")

#print("Stability selection model:")
#cur.model = stability_nonzero
#print(Xnames[cur.model])
#print(paste("LOOCV mean error: ", round(LOOCV(data, cur.model)[[1]],3), "    LOOCV median error:", round(LOOCV(data, cur.model)[[2]],3)))




font = 1
par(las=1,lwd=1,cex.main=font, cex.lab=font, cex.axis=font)
options(scipen=5)

### Plots
dev.off()
selected06 = which(output06$relfreq.final>0.9)
selected1 = which(output1$relfreq.final>0.9)
win.graph(width = 12, height = 12)
par(mfrow=c(1,2))
plot(wrapper_output06$values,pch=20,xlab="Iteration",ylab="EBIC value",cex.lab=1, main ="EBIC with gamma=0.6")
abline(h=EBIC(data,selected06,const=0.6),col="red")
plot(wrapper_output1$values,pch=20,xlab="Iteration",ylab="EBIC value",cex.lab=1, main ="EBIC with gamma=1")
abline(h=EBIC(data,selected1,const=1),col="red")





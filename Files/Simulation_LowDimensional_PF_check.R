
### Reproduce results for low-dimensional simulation study (see Section 5.1 of paper) ############################

# Load main functions (Please adjust to the right directory!)
# setwd("C://...") 
# source("AdaSub_main_functions.R")

# Specify correlation structure between explanatory variables ("toeplitz", "global", "block")
mychoice <- menu( c("Toepliz c=0", "Global c=0.7", "Toepliz c=0.9", "Block c=0.5, b=10   "), graphics=TRUE, title="Choose correlation structure" )

if (mychoice==1){ 
  corr.type="global"
  corr=0 
} 

if (mychoice==2){
  corr.type="global"
  corr=0.7 
}

if (mychoice==3){
  corr.type="toeplitz"
  corr=0.9 
}

if (mychoice==4){
  corr.type="block"
  corr=0.5 
  blocks=10 # specify number of blocks in "block" correlation structure
}


choices <-c( "100 examples for each of 9 different sample sizes (as in paper)", "Custom number of examples and sample sizes")
mychoice <- menu( choices , graphics=TRUE, title="How many examples and different sample sizes should be considered?" )

if (mychoice==1) {
  N_Ex=100 
  Scenarios = 9
}
if (mychoice==2) {
  N_Ex = -1
  while(N_Ex < 1 ){
    N_Ex <- readline("Enter the number of examples per sample size: ")
    N_Ex <- ifelse(grepl("\\D",N_Ex),-1,as.integer(N_Ex))
    if(is.na(N_Ex)){break}  # breaks when hit enter
  }
  Scenarios = -1
  while(Scenarios < 1 ){
    Scenarios <- readline("Enter the number of scenarios: \n(no. of different sample sizes n starting from n=40 in steps of size 20) ")
    Scenarios <- ifelse(grepl("\\D",Scenarios),-1,as.integer(Scenarios))
    if(is.na(Scenarios)){break}  # breaks when hit enter
  }
}

p=30 # number of explanatory variables

n_test=100 # Sample size for independent test set (->prediction error)

Iterations=2000 # Number of iterations in AdaSub
q=5 # Initial expected search size in AdaSub 
#K=n later in loop (for different sample sizes n)

const=0 # constant gamma in EBIC 

#N_Ex=100 # number of replicates per scenario
#Scenarios = 9 # number of scenarios (different sample sizes n)
stepsize = 20 # stepsize for increasing sample size n
n.start = 40 # sample size for first scenario (note that number of variables is given by p=10n)


################################################################
## Execute all of the following code in oder to obtain results #
################################################################

AdaSub_false_pos_meanF = numeric(Scenarios)
BIC_Best_false_pos_meanF = numeric(Scenarios)
AdaSub_false_neg_meanF = numeric(Scenarios)
BIC_Best_false_neg_meanF = numeric(Scenarios)

Best_AdaSub_false_pos_meanF = numeric(Scenarios)
Best_AdaSub_false_neg_meanF = numeric(Scenarios)

AdaSub_estimation_error_mean = numeric(Scenarios)
BIC_Best_estimation_error_mean = numeric(Scenarios)
Best_AdaSub_estimation_error_mean = numeric(Scenarios)

AdaSub_test_prediction_error_mean = numeric(Scenarios)
BIC_Best_test_prediction_error_mean = numeric(Scenarios)
Best_AdaSub_test_prediction_error_mean = numeric(Scenarios)

EqualCount_Thres = numeric(Scenarios) # number of times thresholded AdaSub model equals best BIC model
EqualCount_Best =numeric(Scenarios)   # number of times best AdaSub model equals best BIC model

EqualCount_Thres0 = numeric(Scenarios) # number of times thresholded AdaSub model equals true model
EqualCount_Best0 =numeric(Scenarios)   # number of times best AdaSub model equals true model
EqualCount_BestBIC0 =numeric(Scenarios) # number of times best BIC model equals true model

AdaSub_Times = numeric(N_Ex*Scenarios)
FullEnum_Times = numeric(N_Ex*Scenarios)

## Additional info for AdaSUb: Check for PF and OIP 
AdaSub_add_info = list()

set.seed(22)
par(las=1,lwd=1,cex.main=1, cex.lab=1, cex.axis=1)

for (c in 1:Scenarios){
  
  n= n.start + stepsize*(c-1)
  #print(n)
  
  #########
  K = n
  #K = 100
  #########
  
  AdaSub_Thres9_List <- vector("list",N_Ex)
  AdaSub_Best_List <- vector("list",N_Ex)
  BIC_Best_List <- vector("list",N_Ex)
  
  AdaSub_false_pos = numeric(N_Ex)
  BIC_Best_false_pos = numeric(N_Ex)
  AdaSub_false_neg = numeric(N_Ex)
  BIC_Best_false_neg = numeric(N_Ex)
  
  Best_AdaSub_false_pos =numeric(N_Ex)
  Best_AdaSub_false_neg =numeric(N_Ex)
  
  AdaSub_estimation_error = numeric(N_Ex)
  BIC_Best_estimation_error = numeric(N_Ex)
  Best_AdaSub_estimation_error = numeric(N_Ex)
  
  AdaSub_test_prediction_error = numeric(N_Ex)
  BIC_Best_test_prediction_error = numeric(N_Ex)
  Best_AdaSub_test_prediction_error = numeric(N_Ex)
  
  AdaSub_add_info[[c]] = list()
  
  for (k in 1:N_Ex){
    s0 = sample(0:10,1)
    beta1=numeric(p)
    S0=sort (sample(1:p,s0))
    beta1[S0]= runif(s0,-2,2)
    
    if (corr.type=="global"){
      data      = simdata.global.corr(n,p,beta=beta1,sigma.normal=1,corr= corr) 
      data.test = simdata.global.corr(n_test,p,beta=beta1,sigma.normal=1,corr= corr)
    }
    if (corr.type=="block"){
      data      = simdata.block.corr(n,p,beta=beta1,sigma.normal=1,corr= corr, blocks=blocks) 
      data.test = simdata.block.corr(n_test,p,beta=beta1,sigma.normal=1,corr= corr, blocks=blocks) 
    }
    if (corr.type=="toeplitz"){
      data      = simdata.toeplitz.corr(n,p,beta=beta1,sigma.normal=1,corr= corr) 
      data.test = simdata.toeplitz.corr(n_test,p,beta=beta1,sigma.normal=1,corr= corr)
    }
    
    start.time <- Sys.time()
    output=AdaSub(data,Iter=Iterations,K=K,q=q,criterion="EBIC",const=const,U_C=40, plot.text=paste("\n n = ", n,", Example ", k))
    end.time <- Sys.time()
    AdaSub_Times[(c-1)*N_Ex+k]=difftime(end.time, start.time,units="secs")
    
    # Full Enumeration Result: 
    start.time <- Sys.time()
    regs=summary(regsubsets(as.matrix(data$x),data$y,intercept=TRUE,nvmax=30))
    end.time <- Sys.time()
    FullEnum_Times[(c-1)*N_Ex+k]=difftime(end.time, start.time,units="secs")
    
    modelmatrix=as.matrix(regs$which[,-1])
    vec=rep(0,nrow(modelmatrix)+1)
    vec[nrow(modelmatrix)+1] = EBIC(data,NULL,const)
    for (j in 1:nrow(modelmatrix)) {
      indices.help=as.vector(which(modelmatrix[j,]==TRUE))
      vec[j]=EBIC(data,indices.help,const)
    }
    mini = which.min(vec)
    if (mini<= nrow(modelmatrix)) { BestBIC=as.vector(which(modelmatrix[mini,]==TRUE)) } else { BestBIC = integer(0) }
    
    
    BIC_Best_List[[k]] <- BestBIC
    AdaSub_Thres9_List[[k]] <- which(output$relfreq.final>0.9)
    AdaSub_Best_List[[k]] <- output$best.S
    
    if (setequal(which(output$relfreq.final>0.9),BestBIC))  EqualCount_Thres[c] = EqualCount_Thres[c]+1
    if (setequal(output$best.S,BestBIC))  EqualCount_Best[c] = EqualCount_Best[c]+1
    
    if (setequal(which(output$relfreq.final>0.9),S0))  EqualCount_Thres0[c] = EqualCount_Thres0[c]+1
    if (setequal(output$best.S,S0))  EqualCount_Best0[c] = EqualCount_Best0[c]+1
    if (setequal(BestBIC,S0))  EqualCount_BestBIC0[c] = EqualCount_BestBIC0[c]+1
    
    AdaSub_false_pos[k] = sum(beta1[AdaSub_Thres9_List[[k]]]==0)
    BIC_Best_false_pos[k] = sum(beta1[BIC_Best_List[[k]]]==0)
    Best_AdaSub_false_pos[k] = sum(beta1[AdaSub_Best_List[[k]]]==0)
    
    AdaSub_false_neg[k] = s0-sum(beta1[AdaSub_Thres9_List[[k]]]!=0)
    BIC_Best_false_neg[k] = s0-sum(beta1[BIC_Best_List[[k]]]!=0)
    Best_AdaSub_false_neg[k] = s0-sum(beta1[AdaSub_Best_List[[k]]]!=0)
    
    AdaSub_beta.hat = beta.hat(data,which(output$relfreq.final>0.9))
    BIC_Best_beta.hat = beta.hat(data, BestBIC)
    Best_AdaSub_beta.hat = beta.hat(data,output$best.S)
    
    AdaSub_test_y.hat = y.hat(data.test,AdaSub_beta.hat)
    BIC_Best_test_y.hat = y.hat(data.test,BIC_Best_beta.hat)
    Best_AdaSub_y.hat = y.hat(data.test,Best_AdaSub_beta.hat)
    
    AdaSub_estimation_error[k] = sqnorm2 ( c(0,beta1) - AdaSub_beta.hat)
    BIC_Best_estimation_error[k] = sqnorm2 ( c(0,beta1) - BIC_Best_beta.hat)
    Best_AdaSub_estimation_error[k] = sqnorm2 ( c(0,beta1) - Best_AdaSub_beta.hat)
    
    AdaSub_test_prediction_error[k]   = sqrt( 1/n_test * sqnorm2 ( data.test$y - AdaSub_test_y.hat))
    BIC_Best_test_prediction_error[k] = sqrt( 1/n_test * sqnorm2 ( data.test$y -  BIC_Best_test_y.hat))
    Best_AdaSub_test_prediction_error[k]   = sqrt( 1/n_test * sqnorm2 ( data.test$y -  Best_AdaSub_y.hat))
    
    ## Check if PF, OIP and OIP' hold 
    AdaSub_add_info[[c]][[k]] <- check_OIP_PF(output = output, p = p, S_star = BestBIC)
    
    cat("n = ", n,", Example ", k, "\n")
    print("True model:")
    print(S0)
    print("Best BIC model:")
    print(BestBIC)
    print("AdaSub thresholded model:")
    print(which(output$relfreq.final>0.9))
    print("AdaSub 'best' model:")
    print(output$best.S)
    print("PF:")
    print(AdaSub_add_info[[c]][[k]]$PF)
    print("OIP:")
    print(AdaSub_add_info[[c]][[k]]$OIP)
    print("OIP_prime:")
    print(AdaSub_add_info[[c]][[k]]$OIP_prime)
    print("#####################################")
  }
  
  AdaSub_false_pos_meanF[c] = mean(AdaSub_false_pos)
  BIC_Best_false_pos_meanF[c] = mean(BIC_Best_false_pos)
  AdaSub_false_neg_meanF[c] = mean(AdaSub_false_neg)
  BIC_Best_false_neg_meanF[c] = mean(BIC_Best_false_neg)
  Best_AdaSub_false_pos_meanF[c] =  mean(Best_AdaSub_false_pos)
  Best_AdaSub_false_neg_meanF[c] = mean(Best_AdaSub_false_neg)
  
  AdaSub_estimation_error_mean[c] = mean(AdaSub_estimation_error)
  BIC_Best_estimation_error_mean[c] = mean(BIC_Best_estimation_error) 
  AdaSub_test_prediction_error_mean[c] = mean(AdaSub_test_prediction_error)
  BIC_Best_test_prediction_error_mean[c] = mean(BIC_Best_test_prediction_error)
  Best_AdaSub_estimation_error_mean[c] = mean(Best_AdaSub_estimation_error)
  Best_AdaSub_test_prediction_error_mean[c] = mean(Best_AdaSub_test_prediction_error)
}



### Compute additional information of AdaSub
PF_count = matrix(NA, Scenarios, N_Ex)
OIP_count = matrix(NA, Scenarios, N_Ex)
OIP_prime_count = matrix(NA, Scenarios, N_Ex)

PF_count_rel = matrix(NA, Scenarios, N_Ex)
OIP_count_rel = matrix(NA, Scenarios, N_Ex)
OIP_prime_count_rel = matrix(NA, Scenarios, N_Ex)


PF_holds_freq = rep(0, Scenarios)
OIP_holds_freq = rep(0, Scenarios)
OIP_prime_holds_freq = rep(0, Scenarios)
Thresholded_subseteq_best_freq = rep(0, Scenarios)
Thresholded_equals_best_freq = rep(0, Scenarios)
for (c in 1:Scenarios) {
  for (k in 1:N_Ex) {
    PF_count[c, k] = length(AdaSub_add_info[[c]][[k]]$PF) 
    OIP_count[c, k] = length(AdaSub_add_info[[c]][[k]]$OIP)  
    OIP_prime_count[c, k] = length(AdaSub_add_info[[c]][[k]]$OIP_prime) 
    if (length(AdaSub_add_info[[c]][[k]]$S_star)>0) PF_count_rel[c, k] = length(AdaSub_add_info[[c]][[k]]$PF) / length(AdaSub_add_info[[c]][[k]]$S_star)
    if (length(AdaSub_add_info[[c]][[k]]$S_star)>0) OIP_count_rel[c, k] = length(AdaSub_add_info[[c]][[k]]$OIP)  / length(AdaSub_add_info[[c]][[k]]$S_star)
    if (length(AdaSub_add_info[[c]][[k]]$S_star)>0) OIP_prime_count_rel[c, k] = length(AdaSub_add_info[[c]][[k]]$OIP_prime)  / length(AdaSub_add_info[[c]][[k]]$S_star)
    PF_holds_freq[c] = PF_holds_freq[c] + AdaSub_add_info[[c]][[k]]$PF_holds
    OIP_holds_freq[c] = OIP_holds_freq[c] + AdaSub_add_info[[c]][[k]]$OIP_holds
    OIP_prime_holds_freq[c] = OIP_prime_holds_freq[c] + AdaSub_add_info[[c]][[k]]$OIP_prime_holds
    Thresholded_subseteq_best_freq[c] = Thresholded_subseteq_best_freq[c] + AdaSub_add_info[[c]][[k]]$Thresholded_subseteq_best
    Thresholded_equals_best_freq[c] = Thresholded_equals_best_freq[c] + AdaSub_add_info[[c]][[k]]$Thresholded_equals_best
  }
}



##### Save ############################################################
Save1 = list ( AdaSub_Times = AdaSub_Times , FullEnum_Times = FullEnum_Times,
               AdaSub_false_pos_mean = AdaSub_false_pos_meanF, BIC_Best_false_pos_mean = BIC_Best_false_pos_meanF, Best_AdaSub_false_pos_mean = Best_AdaSub_false_pos_meanF,
               AdaSub_false_neg_mean = AdaSub_false_neg_meanF, BIC_Best_false_neg_mean = BIC_Best_false_neg_meanF, Best_AdaSub_false_neg_mean = Best_AdaSub_false_neg_meanF,
               EqualCount_Thres0 = EqualCount_Thres0, EqualCount_BestBIC0 = EqualCount_BestBIC0, EqualCount_Best0 = EqualCount_Best0,
               EqualCount_Thres = EqualCount_Thres, EqualCount_Best = EqualCount_Best, 
               AdaSub_estimation_error_mean = AdaSub_estimation_error_mean, BIC_Best_estimation_error_mean = BIC_Best_estimation_error_mean, Best_AdaSub_estimation_error_mean = Best_AdaSub_estimation_error_mean,
               AdaSub_test_prediction_error_mean = AdaSub_test_prediction_error_mean, BIC_Best_test_prediction_error_mean = BIC_Best_test_prediction_error_mean, Best_AdaSub_test_prediction_error_mean = Best_AdaSub_test_prediction_error_mean,
               p=p, n_test = n_test, Iterations = Iterations, q=q, corr = corr, const = const, N_Ex = N_Ex, Scenarios = Scenarios, stepsize = stepsize , n.start = n.start,
               AdaSub_add_info = AdaSub_add_info, 
               PF_holds_freq = PF_holds_freq, 
               OIP_holds_freq = OIP_holds_freq, 
               OIP_prime_holds_freq = OIP_prime_holds_freq,
               Thresholded_subseteq_best_freq = Thresholded_subseteq_best_freq,
               Thresholded_equals_best_freq = Thresholded_equals_best_freq,
               PF_count = PF_count,
               OIP_count = OIP_count,
               OIP_prime_count = OIP_prime_count,
               PF_count_rel = PF_count_rel,
               OIP_count_rel = OIP_count_rel,
               OIP_prime_count_rel = OIP_prime_count_rel)

setwd("C://Users//Staerk//scieboBonn//AdaSub paper for EJS//R simulation saves")
save(Save1, file="All_info_BIC_Global00_p30_simulations100.RData")

#########################################################################
#load("BIC_Toeplitz_09_p30_simulations100.RData")

All_info_EBIC06_Toeplitz09_p_increasing_400_2000_simulations500


###########################################################################################

color_AdaSub = "black"
color_AdaSubBest = "orange"
color_BIC = "blue"

legendnames = c("AdaSubThres", "AdaSubBest", "Best Subset BIC")

pchs = c(4,8,1)
cex.size = 1.5
ltys = c(1,2,3)


font = 2
par(las=1,lwd=2,cex.main=font, cex.lab=font, cex.axis=font)

dev.off()
win.graph(width = 14, height = 8, pointsize = 8)

############################################################################
#m <- matrix(c(1,2,3,4,5,3,6,7,3),nrow = 3,ncol = 3,byrow = FALSE)
#layout(mat = m, heights = c(0.45,0.45,0.1))
############################################################################

m <- matrix(c(1,2,3,4,1,5,6,4,1,7,8,4),nrow = 4,ncol = 3,byrow = FALSE)
layout(mat = m, heights = c(0.05,0.425,0.425,0.1))


par(mar=c(0, 0, 2, 0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="", cex.main=1.5,
     main=paste("Low-dimensional simulations for BIC with",N_Ex, "simulated datasets for", Scenarios, "scenarios"))
par(mar=c(5.1, 4.1, 4.1, 2.1))



font = 2
par(las=1,lwd=2,cex.main=font, cex.lab=font, cex.axis=font)



xaxis = n.start  + stepsize * ((1:Scenarios)-1)
#plot(xaxis,Save1$AdaSub_false_pos_mean,pch=pchs[1],cex=cex.size,ylim=c(0,max(Save1$AdaSub_false_pos_mean,Save1$Best_AdaSub_false_pos_mean,Save1$BIC_Best_false_pos_mean)), xlab="n",ylab="",main="Mean false positives",xaxt="n",col=color_AdaSub, type="b", lty = ltys[1])
plot(xaxis,Save1$AdaSub_false_pos_mean,pch=pchs[1],cex=cex.size,ylim=c(0,6), xlab="n",ylab="",main="Mean false positives",xaxt="n",col=color_AdaSub, type="b", lty = ltys[1])
points(xaxis,Save1$Best_AdaSub_false_pos_mean,pch=pchs[2],cex=cex.size,col=color_AdaSubBest, type="b", lty = ltys[2])
points(xaxis,Save1$BIC_Best_false_pos_mean,pch=pchs[3],cex=cex.size,col=color_BIC, type="b", lty = ltys[3])
axis(1,at=xaxis)

#plot(xaxis,Save1$AdaSub_false_neg_mean,pch=pchs[1],cex=cex.size,ylim=c(0,max(Save1$AdaSub_false_neg_mean,Save1$Best_AdaSub_false_neg_mean,Save1$BIC_Best_false_neg_mean)), xlab="n",ylab="",main="Mean false negatives",xaxt="n",col=color_AdaSub, type="b", lty = ltys[1])
plot(xaxis,Save1$AdaSub_false_neg_mean,pch=pchs[1],cex=cex.size,ylim=c(0,6), xlab="n",ylab="",main="Mean false negatives",xaxt="n",col=color_AdaSub, type="b", lty = ltys[1])
points(xaxis,Save1$Best_AdaSub_false_neg_mean,pch=pchs[2],cex=cex.size,col=color_AdaSubBest, type="b", lty = ltys[2])
points(xaxis,Save1$BIC_Best_false_neg_mean,pch=pchs[3],cex=cex.size,col=color_BIC, type="b", lty = ltys[3])
axis(1,at=xaxis)



#legend_order <- matrix(1:8,ncol=4,byrow = TRUE)
par(mar=c(0,0,0,0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
plot_colors <- c(color_AdaSub, color_AdaSubBest, color_BIC)
legend(x = "center",inset = 0,
       legend = legendnames , 
       col=plot_colors, cex=2, box.lty=1, ncol = 4,
       pch=pchs, pt.bg = 'white', lty=ltys) #,
par(mar=c(5.1, 4.1, 4.1, 2.1))



plot(xaxis,Save1$EqualCount_Thres0/N_Ex,pch=pchs[1],cex=cex.size,ylim=c(0,1), xlab="n",ylab="",main="Rel. freq. true model selected",xaxt="n",col=color_AdaSub, type="b", lty = ltys[1])
axis(1,at=xaxis)
points(xaxis,Save1$EqualCount_Best0/N_Ex,pch=pchs[2],cex=cex.size,col=color_AdaSubBest, type="b", lty = ltys[2])
points(xaxis,Save1$EqualCount_BestBIC0/N_Ex,pch=pchs[3],cex=cex.size,col=color_BIC, type="b", lty = ltys[3])

plot(xaxis,Save1$EqualCount_Thres/N_Ex,pch=pchs[1],cex=2.5,ylim=c(0,1), xlab="n",ylab="",main="Rel. freq. AdaSub yields best BIC",xaxt="n",col=color_AdaSub, type="b", lty = ltys[1])
axis(1,at=xaxis)
points(xaxis,Save1$EqualCount_Best/N_Ex,pch=pchs[2],ylim=c(0,1),cex=cex.size,col=color_AdaSubBest, type="b", lty = ltys[2])


plot(xaxis,Save1$AdaSub_estimation_error_mean,ylim=range(Save1$AdaSub_estimation_error_mean,Save1$BIC_Best_estimation_error_mean),col=color_AdaSub, type="b", lty = ltys[1],
     pch=pchs[1],cex=cex.size,xaxt="n",xlab="n",ylab="",main="Mean estimation error (MSE)")
points(xaxis,Save1$Best_AdaSub_estimation_error_mean,pch=pchs[2],cex=cex.size,col=color_AdaSubBest, type="b", lty = ltys[2])
points(xaxis,Save1$BIC_Best_estimation_error_mean,pch=pchs[3],cex=cex.size,col=color_BIC, type="b", lty = ltys[3])
axis(1,at=xaxis)

plot(xaxis,Save1$AdaSub_test_prediction_error_mean,ylim=range(Save1$AdaSub_test_prediction_error_mean,Save1$BIC_Best_test_prediction_error_mean), type="b", lty = ltys[1],
     pch=pchs[1],cex=cex.size,xaxt="n",xlab="n",ylab="",main="Mean prediction error (RMSE)")
points(xaxis,Save1$Best_AdaSub_test_prediction_error_mean,pch=pchs[2],cex=cex.size,col=color_AdaSubBest, type="b", lty = ltys[2])
points(xaxis,Save1$BIC_Best_test_prediction_error_mean,pch=pchs[3],cex=cex.size,col=color_BIC, type="b", lty = ltys[3])
axis(1,at=xaxis)

par(las=1,lwd=1,cex.main=1, cex.lab=1, cex.axis=1)



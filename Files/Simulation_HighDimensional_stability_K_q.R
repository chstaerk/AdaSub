
### Reproduce results for high-dimensional simulation study (see Section 5.2 of paper) ############################

# Load main functions (Please adjust to the right directory!)
# setwd("C://Users//Staerk//sciebo//AdaSub paper for S&C//R files//R files for submission to S&C") 
# source("AdaSub_main_functions.R")

# Specify correlation structure between explanatory variables ("toeplitz", "global", "block")
mychoice <- menu( c("Toepliz c=0", "Toepliz c=0.9"), graphics=TRUE, title="Choose correlation structure" )

if (mychoice==1){ 
  corr.type="global"
  corr=0 
} 

#if (mychoice==2){
#  corr.type="global"
#  corr=0.7 
#}

if (mychoice==2){
  corr.type="toeplitz"
  corr=0.9 
}

#if (mychoice==4){
#  corr.type="block"
#  corr=0.5 
#  blocks=10 # specify number of blocks in "block" correlation structure
#}


choices <-c( "20 examples for each of 9 different sample sizes (as in paper)", "Custom number of examples and sample sizes")
mychoice <- menu( choices , graphics=TRUE, title="How many examples and different sample sizes should be considered?" )

if (mychoice==1) {
  N_Ex=20 
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

choices <-c( "EBIC constant gamma = 0.6", "EBIC constant gamma = 1")
mychoice <- menu( choices , graphics=TRUE, title="Which EBIC constant (gamma) should be used?" )

if (mychoice==1) {
  const=0.6 
}
if (mychoice==2) {
  const = 1
}

n_test=100 # Sample size for independent test set (->prediction error)

Iterations=5000 # Number of iterations in AdaSu

################
# const = 1
###############

#N_Ex=100 # number of replicates per scenario
#Scenarios = 9 # number of scenarios (different sample sizes n)
stepsize = 20 # stepsize for increasing sample size n
n.start = 40 # sample size for first scenario (note that number of variables is given by p=10n)

#n = 100
#p = 1000

N_Rep = 10
#N_Ex = 2

p_max = 10 *(n.start + (Scenarios-1)*stepsize)

#AdaSubBest_results = matrix(0,nrow = N_Ex, ncol = p)
#AdaSub_results = matrix(0,nrow = N_Ex, ncol = p)


#############################################################################################
set.seed(123)
#############################################################################################

AdaSubBest_results = array(0, dim=c(Scenarios, N_Ex, N_Rep, p_max))
AdaSub_results =  array(0, dim=c(Scenarios, N_Ex, N_Rep, p_max))

AdaSubBest_false_pos = array(0, dim=c(Scenarios, N_Ex, N_Rep))
AdaSub_false_pos = array(0, dim=c(Scenarios, N_Ex, N_Rep))
AdaSubBest_false_neg = array(0, dim=c(Scenarios, N_Ex, N_Rep))
AdaSub_false_neg = array(0, dim=c(Scenarios, N_Ex, N_Rep))

comptimes = array(0,dim=c(Scenarios, N_Ex, N_Rep))
nb.iterations = array(0,dim=c(Scenarios, N_Ex, N_Rep))

q.array = array(0,dim=c(Scenarios, N_Ex, N_Rep))
K.array = array(0,dim=c(Scenarios, N_Ex, N_Rep))


par(las=1,lwd=1,cex.main=1, cex.lab=1, cex.axis=1)

for (c in 1:Scenarios){
  
  n= n.start + stepsize*(c-1)
  p = 10 * n # Growing p


  for (k in 1:N_Ex){
    
    
    # s0 = 5
    #S0 = 1:5
    #S0 = c(1,11,21,31,41)
    #beta1[S0]= c(0.4, 0.8, 1.2, 1.6, 2)
    
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
  

  values.matrix = matrix(0, nrow = N_Rep, ncol = Iterations) 
  
  for (l in 1:N_Rep){
    
    #K = n
    K = runif(1, min = 0.5*n, max = 2*n)
      # K = runif(1,min = 1, max = 500)
    #K = runif(1, min = 0.01*n, max = 10*n)
    #K = rexp(1, rate =1/n)
    #q = runif(1, min = 2, max = 10) # max = 20
    q = runif(1, min = 5, max = 15)
    
    q.array[c,k,l] = q 
    K.array[c,k,l] = K
    
    ## Adaptive Subspace Method
    start.time <- Sys.time()
    output=AdaSub(data,Iter=Iterations,K=K,q=q,criterion="EBIC",const=const,U_C=40,plot.text=paste("\n n =", n,", Example", k, ", Run", l))
    end.time <- Sys.time()
    time = as.double(end.time-start.time,units = "secs")
    comptimes[c,k,l] = time
    #nb.iterations[c,k,l] = which.min(output$values)
    
    model = which(output$relfreq.final>0.9)
    AdaSub_model = model
    
    model = output$best.S
    AdaSubBest_model = model
    
    values.matrix[l,] = output$values
    
    #######################################################################
    
    cat("n = ", n,", Example ", k, ", Run", l,  "\n")
    cat("q = ", q, "\n")
    cat("K = ", K, "\n")
    print("True model:")
    print(S0)
    print("AdaSub thresholded model:")
    print(AdaSub_model)
    #print(which(output$relfreq.final>0.9))
    print("AdaSub 'best' model:")
    print(AdaSubBest_model)
    print("#####################################")
    
    AdaSubBest_results[c, k, l, AdaSubBest_model] = 1
    AdaSub_results[c, k, l, AdaSub_model] = 1
    
    AdaSubBest_false_pos[c, k, l] = sum(beta1[AdaSubBest_model]==0)
    AdaSubBest_false_neg[c, k, l] = s0-sum(beta1[AdaSubBest_model]!=0)
    AdaSub_false_pos[c, k, l] = sum(beta1[AdaSub_model]==0)
    AdaSub_false_neg[c, k, l] = s0-sum(beta1[AdaSub_model]!=0)
    
  }
  min.value = min(values.matrix)
  for (l in 1:N_Rep){ 
    if (min(values.matrix[l,])==min.value){
      nb.iterations[c,k,l] = min(which(values.matrix[l,]==min.value)) 
     } else {
        nb.iterations[c,k,l] = Iterations + 1
       }
  }
}
}



#plot(c(q.array),c(comptimes))
#plot(c(K.array),c(comptimes))

##### Save ############################################################
Save1 = list (  n_test = n_test, Iterations = Iterations, corr = corr, const = const, corr.type = corr.type , 
                 N_Ex = N_Ex, N_Rep = N_Rep, Scenarios = Scenarios, stepsize = stepsize , n.start = n.start, 
                AdaSubBest_results = AdaSubBest_results , AdaSub_results = AdaSub_results ,
                q.array = q.array, K.array = K.array, nb.iterations = nb.iterations, comptimes = comptimes)

#setwd("C://Users//Staerk//sciebo//AdaSub paper for S&C//R simulation saves")
# save(Save1, file="Stability_EBIC1_Global00_p_increasing_600_2000_simulations20.RData")
#########################################################################
# load("Stability_EBIC06_Toeplitz09_p_increasing_600_1400_simulations100.RData")

# AdaSub_results = Save1$AdaSub_results
# AdaSubBest_results = Save1$AdaSubBest_results


####################################################################


Count_ident = function(Save1){
  AdaSub_ident = matrix(NA, nrow = Save1$Scenarios, ncol = Save1$N_Ex)
  AdaSubBest_ident = matrix(NA, nrow = Save1$Scenarios, ncol = Save1$N_Ex)
  AdaSub_var = matrix(NA, nrow = Save1$Scenarios, ncol = Save1$N_Ex)
  AdaSubBest_var = matrix(NA, nrow = Save1$Scenarios, ncol = Save1$N_Ex)
  for (c in 1:Save1$Scenarios){
    for(k in 1:Save1$N_Ex){
      AdaSub_ident_count = numeric(Save1$N_Rep)
      AdaSubBest_ident_count = numeric(Save1$N_Rep)
      AdaSub_model_size = numeric(Save1$N_Rep)
      AdaSubBest_model_size = numeric(Save1$N_Rep)
      for (l in 1:Save1$N_Rep){ 
        AdaSub_cur.sum = 0
        AdaSubBest_cur.sum = 0
        for (m in 1:Save1$N_Rep){
          if (all(Save1$AdaSub_results[c,k,l,] == Save1$AdaSub_results[c,k,m,] )) AdaSub_cur.sum = AdaSub_cur.sum + 1
          if (all(Save1$AdaSubBest_results[c,k,l,] == Save1$AdaSubBest_results[c,k,m,] )) AdaSubBest_cur.sum = AdaSubBest_cur.sum + 1
        }
        AdaSub_ident_count[l] = AdaSub_cur.sum
        AdaSubBest_ident_count[l] = AdaSubBest_cur.sum
        AdaSub_model_size[l] = sum(Save1$AdaSub_results[c,k,l,]) 
        AdaSubBest_model_size[l] = sum(Save1$AdaSubBest_results[c,k,l,]) 
      }
      AdaSub_ident[c,k] = max(AdaSub_ident_count)
      AdaSubBest_ident[c,k] = max(AdaSubBest_ident_count)
      AdaSub_var[c,k] = var(AdaSub_model_size) 
      AdaSubBest_var[c,k] = var(AdaSubBest_model_size)
    }
  }
  return(list(AdaSub_ident=AdaSub_ident, AdaSubBest_ident=AdaSubBest_ident, AdaSub_var = AdaSub_var, AdaSubBest_var = AdaSubBest_var ))
}



ident = Count_ident(Save1)
AdaSub_ident  = ident$AdaSub_ident
AdaSubBest_ident  = ident$AdaSubBest_ident
AdaSub_var  = ident$AdaSub_var
AdaSubBest_var  = ident$AdaSubBest_var

pchs = c(4,8,1,2)
cex.size = 1.4
ltys = c(1,4,3,2)
colors = c("black","gray40","orange","darkorange2")

dev.off()
win.graph(width = 12, height = 12)
###################################
#font = 1.5
font = 1
###################################
par(las=1,lwd=2,cex.main=font, cex.lab=font, cex.axis=font)

##################################################################
#m <- matrix(c(1,2,3),nrow = 3,ncol = 1,byrow = TRUE)
#layout(mat = m, heights = c(0.45,0.45,0.1))
##################################################################

m <- matrix(c(1,2,3,4),nrow = 4,ncol = 1,byrow = TRUE)
layout(mat = m, heights = c(0.05,0.425,0.425,0.1))

par(mar=c(0, 0, 2, 0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="", cex.main=1,
     main=paste("Sensitivity analysis regarding algorithmic stability of AdaSub \n for correlation c =", corr,"for EBIC with gamma =",const,"with",N_Ex, "simulated datasets for", Scenarios, "scenarios"))
par(mar=c(5.1, 4.1, 4.1, 2.1))

xaxis = Save1$n.start  + Save1$stepsize * ((1:Save1$Scenarios)-1)
plot(xaxis,rowMeans(AdaSub_ident)/Save1$N_Rep,ylim = c(0,1),pch=pchs[1],cex=cex.size, xlab="n",ylab="",main="Mean rel. freq. of agreement",xaxt="n",col=colors[1], type="b", lty = ltys[1])
points(xaxis,rowMeans(AdaSubBest_ident)/Save1$N_Rep,pch=pchs[3],cex=cex.size,col=colors[3], type="b", lty = ltys[3])
axis(1,at=xaxis)

xaxis = Save1$n.start  + Save1$stepsize * ((1:Save1$Scenarios)-1)
plot(xaxis,rowMeans(AdaSub_var),ylim = c(0,1),pch=pchs[1],cex=cex.size, xlab="n",ylab="",main="Mean variance of model size",xaxt="n",col=colors[1], type="b", lty = ltys[1])
points(xaxis,rowMeans(AdaSubBest_var),pch=pchs[3],cex=cex.size,col=colors[3], type="b", lty = ltys[3])
axis(1,at=xaxis)


legendnames = c("AdaSubThres     ", "AdaSubBest ")

par(mar=c(0,0,0,0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
plot_colors <- colors
legend(x = "center",inset = 0,
       legend = legendnames, 
       col=plot_colors[c(1,3)], cex=font, box.lty=1, ncol = 4,
       pch=pchs[c(1,3)], pt.bg = 'white', lty=ltys[c(1,3)]) #,
par(mar=c(5.1, 4.1, 4.1, 2.1))

par(las=1,lwd=1,cex.main=1, cex.lab=1, cex.axis=1)







































### Reproduce results for sensitivity analysis for effects of q and K (see Section 5.3, Figures 4 and 5 of paper) #######################

# Load main functions (Please adjust to the right directory!)
# setwd("C://Users//Staerk//sciebo//AdaSub paper for S&C//R files//R files for submission to S&C") 
# source("AdaSub_main_functions.R")

# Specify correlation structure between explanatory variables ("toeplitz", "global", "block")
mychoice <- menu( c("Toepliz c=0.9 (as in paper)","Toepliz c=0", "Global c=0.7", "Block c=0.5, b=10   "), graphics=TRUE, title="Choose correlation structure" )

if (mychoice==2){ 
  corr.type="global"
  corr=0 
} 

if (mychoice==3){
  corr.type="global"
  corr=0.7 
}

if (mychoice==1){
  corr.type="toeplitz"
  corr=0.9 
}

if (mychoice==4){
  corr.type="block"
  corr=0.5 
  blocks=10 # specify number of blocks in "block" correlation structure
}


choices <-c( "100 examples (as in paper)", "Custom number of examples")
mychoice <- menu( choices , graphics=TRUE, title="How many examples should be considered?" )

if (mychoice==1) {
  N_Ex=100 
}
if (mychoice==2) {
  N_Ex = -1
  while(N_Ex < 1 ){
    N_Ex <- readline("Enter the number of examples per sample size: ")
    N_Ex <- ifelse(grepl("\\D",N_Ex),-1,as.integer(N_Ex))
    if(is.na(N_Ex)){break}  # breaks when hit enter
  }
}

choices <-c( "EBIC constant gamma = 0.6 (as in paper)", "EBIC constant gamma = 1")
mychoice <- menu( choices , graphics=TRUE, title="Which EBIC constant (gamma) should be used?" )

if (mychoice==1) {
  const=0.6 
}
if (mychoice==2) {
  const = 1
}

n_test=100 # Sample size for independent test set (->prediction error)

Iterations=5000 # Number of iterations in AdaSub

################
# const = 0.6
###############

#N_Ex=100 # number of replicates per scenario
#Scenarios = 9 # number of scenarios (different sample sizes n)
stepsize = 100 # stepsize for increasing sample size n
n.start = 100 # sample size for first scenario (note that number of variables is given by p=10n)

Scenarios = 2

#n = 100
#p = 1000

N_Rep = 10
#N_Ex = 100

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

Ks = c(1,100,200,1000,2000)
qs = c(1,2,5,10,15)

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
    
    #K = runif(1, min = 0.5*n, max = 2*n)
    #q = runif(1, min = 5, max = 15)
    if (l <=5 ) {
      q = 10 
      K = Ks[l]
    }
    
    if (l >5 ) {
      q = qs[l-5] 
      K = n
    }
    
    q.array[c,k,l] = q 
    K.array[c,k,l] = K
    
    ## Adaptive Subspace Method
    start.time <- Sys.time()
    output=AdaSub(data,Iter=Iterations,K=K,q=q,criterion="EBIC",const=const,U_C=40,plot.text=paste("\n n =", n,", Example", k, ", Run", l, ", q =", q,", K =", K))
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


##### Save ############################################################
Save1 = list (  n_test = n_test, Iterations = Iterations, corr = corr, const = const, corr.type = corr.type , 
                 N_Ex = N_Ex, N_Rep = N_Rep, Scenarios = Scenarios, stepsize = stepsize , n.start = n.start, 
                AdaSubBest_results = AdaSubBest_results , AdaSub_results = AdaSub_results ,
                q.array = q.array, K.array = K.array, nb.iterations = nb.iterations, comptimes = comptimes)

# setwd("C://Users//Staerk//sciebo//AdaSub paper for S&C//R simulation saves")
# save(Save1, file="Choice_Kq_EBIC06_Toeplitz09_p_1000_2000_simulations100.RData")
#########################################################################
# load("Choice_Kq_EBIC06_Toeplitz09_p_1000_2000_simulations100.RData")
# nb.iterations = Save1$nb.iterations
# comptimes = Save1$comptimes

Nbs_q_100 = numeric(5)
Nbs_q_200 = numeric(5)
Nbs_K_100 = numeric(5)
Nbs_K_200 = numeric(5)
for (i in 1:5){
  Nbs_q_100[i] = sum(nb.iterations[1,,5+i]>5000)
  Nbs_q_200[i] = sum(nb.iterations[2,,5+i]>5000)
  Nbs_K_100[i] = sum(nb.iterations[1,,i]>5000)
  Nbs_K_200[i] = sum(nb.iterations[2,,i]>5000)
}

nb.iterations[nb.iterations==5001] = 5000

font = 1
par(las=1,lwd=1,cex.main=font, cex.lab=font, cex.axis=font)

cexax = 1
cexlab = 1.2

dev.off()
win.graph(width = 12, height = 12)


m <- matrix(c(1,2,3,1,4,5),nrow = 3,ncol = 2,byrow = FALSE)
layout(mat = m, heights = c(0.05,0.475,0.475))

par(mar=c(0, 0, 2, 0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="", cex.main=1.2,
     main=paste("Sensitivity analysis for effect of q for EBIC with gamma =",const,"with",N_Ex, "simulated datasets for", Scenarios, "scenarios"))
par(mar=c(5.1, 4.1, 4.1, 2.1))

font = 1.5

####################################
#par(mfrow=c(2,2))
####################################

par(las=1,lwd=2,cex.main=font, cex.lab=font, cex.axis=font)

boxplot(nb.iterations[1,,6],nb.iterations[1,,7],nb.iterations[1,,8],nb.iterations[1,,9],nb.iterations[1,,10],
        names=c("1", "2", "5", "10", "15"), xlab = "q", main="No. iterations after 'best' model found (n = 100, K = n)"
        , ylim=c(0,5400),cex.main=1.3,cex.axis=cexax, cex.lab = cexlab)

text( c(1:5) , 5300 , paste("f = ",Nbs_q_100,sep="")  )

boxplot(comptimes[1,,6],comptimes[1,,7],comptimes[1,,8],comptimes[1,,9],comptimes[1,,10],
        names=c("1", "2", "5", "10", "15"), xlab = "q", main="Comp. time (n = 100, K = n)"
        ,ylab="Time (s)",cex.main=1.3,cex.axis=cexax, cex.lab = cexlab)

boxplot(nb.iterations[2,,6],nb.iterations[2,,7],nb.iterations[2,,8],nb.iterations[2,,9],nb.iterations[2,,10],
        names=c("1", "2", "5", "10", "15"), xlab = "q", main="No. iterations after 'best' model found (n = 200, K = n)"
        , ylim=c(0,5400),cex.main=1.3,cex.axis=cexax, cex.lab = cexlab)

text( c(1:5) , 5300 , paste("f = ",Nbs_q_200,sep="")  )

boxplot(comptimes[2,,6],comptimes[2,,7],comptimes[2,,8],comptimes[2,,9],comptimes[2,,10],
        names=c("1", "2", "5", "10", "15"), xlab = "q", main="Comp. time (n = 200, K = n)"
        ,ylab="Time (s)",cex.main=1.3,cex.axis=cexax, cex.lab = cexlab)

######################################################################################################

win.graph(width = 12, height = 12)

m <- matrix(c(1,2,3,1,4,5),nrow = 3,ncol = 2,byrow = FALSE)
layout(mat = m, heights = c(0.05,0.475,0.475))

par(mar=c(0, 0, 2, 0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="", cex.main=1.2,
     main=paste("Sensitivity analysis for effect of K for EBIC with gamma =",const,"with",N_Ex, "simulated datasets for", Scenarios, "scenarios"))
par(mar=c(5.1, 4.1, 4.1, 2.1))

###################################################
#par(mfrow=c(2,2))
###################################################

boxplot(nb.iterations[1,,1],nb.iterations[1,,2],nb.iterations[1,,3],nb.iterations[1,,4],nb.iterations[1,,5],
        names=c("1", "100", "200", "1000", "2000"), xlab = "K", main="No. iterations after 'best' model found (n = 100, q = 10)"
        , ylim=c(0,5400),cex.main=1.3,cex.axis=cexax, cex.lab = cexlab)

text( c(1:5) , 5300 , paste("f = ",Nbs_K_100,sep="")  )

boxplot(comptimes[1,,1],comptimes[1,,2],comptimes[1,,3],comptimes[1,,4],comptimes[1,,5],
        names=c("1", "100", "200", "1000", "2000"), xlab = "K", main="Comp. time (n = 100, q = 10)"
        ,ylab="Time (s)",cex.main=1.3,cex.axis=cexax, cex.lab = cexlab)


boxplot(nb.iterations[2,,1],nb.iterations[2,,2],nb.iterations[2,,3],nb.iterations[2,,4],nb.iterations[2,,5],
        names=c("1", "100", "200", "1000", "2000"), xlab = "K", main="No. iterations after 'best' model found (n = 200, q = 10)"
        , ylim=c(0,5400),cex.main=1.3,cex.axis=cexax, cex.lab = cexlab)

text( c(1:5) , 5300 , paste("f = ",Nbs_K_200,sep="")  )

boxplot(comptimes[2,,1],comptimes[2,,2],comptimes[2,,3],comptimes[2,,4],comptimes[2,,5],
        names=c("1", "100", "200", "1000", "2000"), xlab = "K", main="Comp. time (n = 200, q = 10)"
        ,ylab="Time (s)",cex.main=1.3,cex.axis=cexax, cex.lab = cexlab)

par(las=1,lwd=1,cex.main=1, cex.lab=1, cex.axis=1)




####### Alternative plots #######################################################################


#cexax = 1

#par(mfrow=c(2,2))
#boxplot(nb.iterations[1,,6],nb.iterations[1,,7],nb.iterations[1,,8],nb.iterations[1,,9],nb.iterations[1,,10],
#        names=c("q = 1", "q = 2", "q = 5", "q = 10", "q = 15"), main="Nb. iterations after 'best' model found (n = 100, K = n)"
#        , ylim=c(0,5400),cex.main=1.3,cex.axis=cexax)

#text( c(1:5) , 5300 , paste("f = ",Nbs_q_100,sep="")  )

#boxplot(comptimes[1,,6],comptimes[1,,7],comptimes[1,,8],comptimes[1,,9],comptimes[1,,10],
#        names=c("q = 1", "q = 2", "q = 5", "q = 10", "q = 15"), main="Comp. time (n = 100, K = n)"
#        ,ylab="Time (s)",cex.main=1.3,cex.axis=cexax)

#boxplot(nb.iterations[2,,6],nb.iterations[2,,7],nb.iterations[2,,8],nb.iterations[2,,9],nb.iterations[2,,10],
#        names=c("q = 1", "q = 2", "q = 5", "q = 10", "q = 15"), main="Nb. iterations after 'best' model found (n = 200, K = n)"
#        , ylim=c(0,5400),cex.main=1.3,cex.axis=cexax)

#text( c(1:5) , 5300 , paste("f = ",Nbs_q_200,sep="")  )

#boxplot(comptimes[2,,6],comptimes[2,,7],comptimes[2,,8],comptimes[2,,9],comptimes[2,,10],
#        names=c("q = 1", "q = 2", "q = 5", "q = 10", "q = 15"), main="Comp. time (n = 200, K = n)"
#        ,ylab="Time (s)",cex.main=1.3,cex.axis=cexax)

######################################################################################################

#par(mfrow=c(2,2))
#boxplot(nb.iterations[1,,1],nb.iterations[1,,2],nb.iterations[1,,3],nb.iterations[1,,4],nb.iterations[1,,5],
#        names=c("K = 1", "K = 100", "K = 200", "K = 1000", "K = 2000"), main="Nb. iterations after 'best' model found (n = 100, q = 10)"
#        , ylim=c(0,5400),cex.main=1.3,cex.axis=cexax)

#text( c(1:5) , 5300 , paste("f = ",Nbs_K_100,sep="")  )

#boxplot(comptimes[1,,1],comptimes[1,,2],comptimes[1,,3],comptimes[1,,4],comptimes[1,,5],
#        names=c("K = 1", "K = 100", "K = 200", "K = 1000", "K = 2000"), main="Comp. time (n = 100, q = 10)"
#        ,ylab="Time (s)",cex.main=1.3,cex.axis=cexax)


#boxplot(nb.iterations[2,,1],nb.iterations[2,,2],nb.iterations[2,,3],nb.iterations[2,,4],nb.iterations[2,,5],
#        names=c("K = 1", "K = 100", "K = 200", "K = 1000", "K = 2000"), main="Nb. iterations after 'best' model found (n = 200, q = 10)"
#        , ylim=c(0,5400),cex.main=1.3,cex.axis=cexax)

#text( c(1:5) , 5300 , paste("f = ",Nbs_K_200,sep="")  )

#boxplot(comptimes[2,,1],comptimes[2,,2],comptimes[2,,3],comptimes[2,,4],comptimes[2,,5],
#        names=c("K = 1", "K = 100", "K = 200", "K = 1000", "K = 2000"), main="Comp. time (n = 200, q = 10)"
#        ,ylab="Time (s)",cex.main=1.3,cex.axis=cexax)




















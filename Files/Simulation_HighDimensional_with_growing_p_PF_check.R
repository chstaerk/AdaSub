
### Reproduce results for high-dimensional simulation study (see Section 5.2 of paper) ############################

# Load main functions (Please adjust to the right directory!)
# setwd("C://Users//Staerk//ScieboBonn/AdaSub paper for EJS//R files//Files") 
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


choices <-c( "500 examples for each of 9 different sample sizes (as in paper)", "Custom number of examples and sample sizes")
mychoice <- menu( choices , graphics=TRUE, title="How many examples and different sample sizes should be considered?" )

if (mychoice==1) {
  N_Ex=500 
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

Iterations=5000 # Number of iterations in AdaSub
#q=5 # Initial expected search size in AdaSub
q = 10 # Initial expected search size in AdaSub
#K=n later in loop (for different sample sizes n)

##################################################
# const=1 # constant gamma in EBIC 

#N_Ex=500 # number of replicates per scenario
#Scenarios = 9 # number of scenarios (different sample sizes n)
######################################################

stepsize = 20 # stepsize for increasing sample size n
n.start = 40 # sample size for first scenario (note that number of variables is given by p=10n)

error=1 # for Stability Selection: expected number of type I errors 
Stab_Iterations=100 # Number of subsamples considered in Stability Selection

PC = FALSE # run PC-simple algorithm?

plotting = TRUE # set to FALSE to speed up computations 
savings = 1


Forward_summary = method_init(Scenarios,n.start=n.start,stepsize=stepsize)
AdaSub_summary = method_init(Scenarios,n.start=n.start,stepsize=stepsize)
AdaSubBest_summary = method_init(Scenarios,n.start=n.start,stepsize=stepsize)
Lasso_summary = method_init(Scenarios,n.start=n.start,stepsize=stepsize)
Stability_summary = method_init(Scenarios,n.start=n.start,stepsize=stepsize)
AdaLasso_summary = method_init(Scenarios,n.start=n.start,stepsize=stepsize)
SCAD_summary = method_init(Scenarios,n.start=n.start,stepsize=stepsize)
Tilting_summary = method_init(Scenarios,n.start=n.start,stepsize=stepsize)
PC_summary = method_init(Scenarios,n.start=n.start,stepsize=stepsize)

### Lists which include all results for all scenarios (first entry, matrix of results for first scenario)
Forward_everything = list()
AdaSub_everything = list()
AdaSubBest_everything = list()
Lasso_everything = list()
Stability_everything = list()
AdaLasso_everything = list()
SCAD_everything = list()
Tilting_everything = list()
PC_everything = list()

## Additional info for AdaSUb: Check for PF and OIP 
AdaSub_add_info = list()

#############################################################################################
#set.seed(123)
#set.seed(22)
# if (corr.type=="toeplitz") set.seed(22) # Technical remark: This seed was used for Toeplitz correlation structures  
# [the seed of 123 shows also comparable (and favourable) results for AdaSub in this setting]
#############################################################################################

par(las=1,lwd=1,cex.main=1, cex.lab=1, cex.axis=1)

for (c in 1:Scenarios){
  
  n= n.start + stepsize*(c-1)
  p = 10 * n # Growing p
  #p = 1/4*n^2
  #########
  K= n
  #########
  Forward_scenario = method_init_scenario(N_Ex)
  AdaSub_scenario = method_init_scenario(N_Ex)
  AdaSubBest_scenario = method_init_scenario(N_Ex)
  Lasso_scenario = method_init_scenario(N_Ex)
  Stability_scenario = method_init_scenario(N_Ex)
  AdaLasso_scenario = method_init_scenario(N_Ex)
  SCAD_scenario = method_init_scenario(N_Ex)
  Tilting_scenario = method_init_scenario(N_Ex)
  PC_scenario = method_init_scenario(N_Ex)
  
  AdaSub_add_info[[c]] = list()
  
  
  for (k in 1:N_Ex){
    
    set.seed(123*k*c)
    
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
    
    ## Adaptive Subspace Method
    start.time <- Sys.time()
    output=AdaSub(data,Iter=Iterations,K=K,q=q,criterion="EBIC",const=const,U_C=40, plot.text=paste("\n n =", n,", Example", k), plotting = plotting, savings = savings)
    end.time <- Sys.time()
    time = as.double(end.time-start.time,units = "secs")
    
    ## Check if PF, OIP and OIP' hold 
    AdaSub_add_info[[c]][[k]] <- check_OIP_PF(output = output, p = p, S_star = NA)
    
    #check_AdaSub$PF
    #check_AdaSub$OIP
    #check_AdaSub$OIP_prime
    #check_AdaSub$S_star
    #check_AdaSub$counter_OIP
    #check_AdaSub$counter_OIP_prime
    
    model = which(output$relfreq.final>0.9)
    AdaSub_model = model
    
    AdaSub_scenario[,k] = method_results_scenario(model=AdaSub_model,beta1=beta1,s0=s0,data.test=data.test,n_test=n_test,time=time,data=data)
    if (setequal(AdaSub_model,S0))  AdaSub_summary["equal_count_0",c] = AdaSub_summary["equal_count_0",c] + 1
    
    model = output$best.S
    AdaSubBest_model = model
    
    AdaSubBest_scenario[,k] = method_results_scenario(model=AdaSubBest_model,beta1=beta1,s0=s0,data.test=data.test,n_test=n_test,time=time,data=data)
    if (setequal(AdaSubBest_model,S0))  AdaSubBest_summary["equal_count_0",c] = AdaSubBest_summary["equal_count_0",c] + 1
    
    
    # Forward Stepwise Selection
    start.time <- Sys.time()
      regs=summary(regsubsets(as.matrix(data$x),data$y,intercept=TRUE,nvmax=30,method="forward")) #,nvmax=p)) # Maximal number of Variables in model is 30 

      modelmatrix=as.matrix(regs$which[,-1])
      vec=rep(0,nrow(modelmatrix)+1)
      vec[1] = EBIC(data,NULL,const)
      for (j in 1:nrow(modelmatrix)) {
        indices.help=as.vector(which(modelmatrix[j,]==TRUE))
        vec[j+1]=EBIC(data,indices.help,const)
      }
      # Only continue if EBIC value is improved (i.e. decreased)
      differences = vec[-length(vec)]-vec[-1]
      mini = min(which(differences<0))
      #mini = which.min(vec)
      if (mini> 1) { Forward_model=as.vector(which(modelmatrix[mini-1,]==TRUE)) } else { Forward_model = integer(0) }
    end.time <- Sys.time()
    time = as.double(end.time-start.time,units = "secs")
      
      
    Forward_scenario[,k] = method_results_scenario(model=Forward_model,beta1=beta1,s0=s0,data.test=data.test,n_test=n_test,time=time,data=data)
    if (setequal(Forward_model,S0))  Forward_summary["equal_count_0",c] = Forward_summary["equal_count_0",c] + 1
      
      
    ## Stability Selection with Lasso
    start.time <- Sys.time()
    #res = stabpath(data$y,data$x,steps=Stab_Iterations, size = 0.5)  
    #stability_nonzero = stabsel(res,error=error,type="pfer",pi_thr=0.6)$stable # stabsel(res,error=1/p,type="pcer",pi_thr=0.6)$stable
    res <- stabsel(x = data$x, y = data$y,fitfun = glmnet.lasso, cutoff = 0.6 , PFER = error, sampling.type  =  "SS") # Sampling type is complementary pairs
    stability_nonzero = res$selected
    end.time <- Sys.time()
    time = as.double(end.time-start.time,units = "secs")
    names(stability_nonzero)=NULL
    
    Stability_scenario[,k] = method_results_scenario(model=stability_nonzero,beta1=beta1,s0=s0,data.test=data.test,n_test=n_test,time=time,data=data)
    if (setequal(stability_nonzero,S0))  Stability_summary["equal_count_0",c] = Stability_summary["equal_count_0",c] + 1
    
    
    
    ## Lasso with EBIC_gamma (use least squares estimator with l_0 penalty, not the fitted Lasso coefficients)
    start.time <- Sys.time()
    fit.lasso = glmnet(data$x,data$y,intercept=TRUE,family="gaussian",standardize=FALSE,alpha=1,dfmax=30)
    nonzeros.path = predict(fit.lasso,type="nonzero")
    length.path = length(nonzeros.path)
    vec = numeric(length.path)
    for (j in 1:length.path){
      vec[j] = EBIC(data,nonzeros.path[[j]],const) 
    }
    mini = which.min(vec)
    Lasso_nonzero = nonzeros.path[[mini]]
    if (is.null(Lasso_nonzero)) Lasso_nonzero = integer(0)
    end.time <- Sys.time()
    time = as.double(end.time-start.time,units = "secs")
    
    Lasso_scenario[,k] = method_results_scenario(model=Lasso_nonzero,beta1=beta1,s0=s0,data.test=data.test,n_test=n_test,time=time,data=data)
    if (setequal(Lasso_nonzero,S0))  Lasso_summary["equal_count_0",c] = Lasso_summary["equal_count_0",c] + 1
  
    #set.seed(seedsafe)
    
    ## Adaptive Lasso with Cross-Validation in first step and EBIC in second step
    start.time <- Sys.time()
    cv.fit1=cv.glmnet(data$x,data$y,intercept=TRUE,family="gaussian",standardize=FALSE,alpha=1)
    betas.cv.fit1=as.vector(predict(cv.fit1,type="coefficients",s="lambda.1se"))[-1]
    new_penalty=rep(0,p)
    new_penalty=1/pmax(abs(betas.cv.fit1),.Machine$double.eps)
    AdaLasso.fit2= glmnet(data$x,data$y,penalty.factor=new_penalty,family="gaussian",intercept=TRUE,standardize=FALSE,alpha=1,dfmax=n-3)
    AdaLasso.path <- predict(AdaLasso.fit2,type="coefficients")
    nonzeros.AdaLasso.path = predict(AdaLasso.fit2,type="nonzero")
    length.path = dim(AdaLasso.path)[2]
    vec = numeric(length.path)
    for (j in 1:length.path){
      vec[j] = EBIC_beta(data,AdaLasso.path[,j],const) 
    }
    mini = which.min(vec)
    AdaLasso_nonzero = nonzeros.AdaLasso.path[[mini]]
    if (is.null(AdaLasso_nonzero)) AdaLasso_nonzero = integer(0)
    AdaLasso_beta = AdaLasso.path[,mini]
    end.time <- Sys.time()
    time = as.double(end.time-start.time,units = "secs")
    
    AdaLasso_scenario[,k] = method_results_scenario(model=AdaLasso_nonzero,beta1=beta1,s0=s0,data.test=data.test,n_test=n_test,time=time,beta=AdaLasso_beta,data=data)
    if (setequal(AdaLasso_nonzero,S0))  AdaLasso_summary["equal_count_0",c] = AdaLasso_summary["equal_count_0",c] + 1
    
    
    ## SCAD
    start.time <- Sys.time()
    scad.fit <- ncvreg(data$x, data$y, family="gaussian", penalty="SCAD", alpha=1 )  # penalty="MCP"
    beta.path <- predict(scad.fit,type="coefficients")
    length.path = dim(beta.path)[2]
    vec = numeric(length.path)
    for (j in 1:length.path){
      vec[j] = EBIC_beta(data,beta.path[,j],const) 
    }
    mini = which.min(vec)
    SCAD_nonzero = which(beta.path[-1,mini]!= 0)
    names(SCAD_nonzero) =NULL
    if (is.null(SCAD_nonzero)) SCAD_nonzero = integer(0)
    SCAD_beta = beta.path[,mini] 
    end.time <- Sys.time()
    time = as.double(end.time-start.time,units = "secs")
    
    SCAD_scenario[,k] = method_results_scenario(model=SCAD_nonzero,beta1=beta1,s0=s0,data.test=data.test,n_test=n_test,time=time,data=data) #,beta=SCAD_beta
    if (setequal(SCAD_nonzero,S0))  SCAD_summary["equal_count_0",c] = SCAD_summary["equal_count_0",c] + 1
    
    
    if (p<=1200){
    ## Tilting 
    start.time <- Sys.time()
    tilt<-tilting(data$x, data$y, op=2, bic.gamma=1, max.count = 10) #, max.count = 10) # ,max.size = 2) # bic.gamma=1 #,bic.gamma=const) # Tilting with the EBIC constant 0.6 heavily overfits
    Tilting_nonzero = sort(tilt$active.hat)
    end.time <- Sys.time()
    time = as.double(end.time-start.time,units = "secs")
    
    Tilting_scenario[,k] = method_results_scenario(model=Tilting_nonzero,beta1=beta1,s0=s0,data.test=data.test,n_test=n_test,time=time,data=data)
    if (setequal(Tilting_nonzero,S0))  Tilting_summary["equal_count_0",c] = Tilting_summary["equal_count_0",c] + 1
    }
    
    if (PC==T){
    ## PC-Algo
    start.time <- Sys.time()
    PC_nonzero <- which(pcSelect(data$y, data$x, alpha=0.05)$G)
    end.time <- Sys.time()
    time = as.double(end.time-start.time,units = "secs")
    
    PC_scenario[,k] = method_results_scenario(model=PC_nonzero,beta1=beta1,s0=s0,data.test=data.test,n_test=n_test,time=time,data=data)
    if (setequal(PC_nonzero,S0))  PC_summary["equal_count_0",c] = PC_summary["equal_count_0",c] + 1
    }
    
    #######################################################################
    
    cat("n = ", n,", Example ", k, "\n")
    print("True model:")
    print(S0)
    print("Stepwise model:")
    print(Forward_model)
    print("AdaSub thresholded model:")
    print(AdaSub_model)
    #print(which(output$relfreq.final>0.9))
    print("AdaSub 'best' model:")
    print(AdaSubBest_model)
    print("Stability Sel. model:")
    print(stability_nonzero)
    print("Lasso model:")
    print(Lasso_nonzero)
    print("AdaLasso model:")
    print(AdaLasso_nonzero)
    print("SCAD model:")
    print(SCAD_nonzero)
    if (p<=1200){
    print("Tilting model:")
    print(Tilting_nonzero)
    }
    if (PC==T){
    print("PC-simple model:")
    print(PC_nonzero)
    }
    print("#####################################")
  }
  
  Forward_summary[1:5,c] = 1/N_Ex*rowSums(Forward_scenario)
  AdaSub_summary[1:5,c] = 1/N_Ex*rowSums(AdaSub_scenario)
  AdaSubBest_summary[1:5,c] = 1/N_Ex*rowSums(AdaSubBest_scenario)
  Lasso_summary[1:5,c] = 1/N_Ex*rowSums(Lasso_scenario)
  Stability_summary[1:5,c] = 1/N_Ex*rowSums(Stability_scenario)
  AdaLasso_summary[1:5,c] = 1/N_Ex*rowSums(AdaLasso_scenario)
  SCAD_summary[1:5,c] = 1/N_Ex*rowSums(SCAD_scenario)
  Tilting_summary[1:5,c] = 1/N_Ex*rowSums(Tilting_scenario)
  PC_summary[1:5,c] = 1/N_Ex*rowSums(PC_scenario)
  
  
  Forward_everything[[c]] = Forward_scenario
  AdaSub_everything[[c]] = AdaSub_scenario
  AdaSubBest_everything[[c]] = AdaSubBest_scenario
  Lasso_everything[[c]] = Lasso_scenario
  Stability_everything[[c]] = Stability_scenario
  AdaLasso_everything[[c]] = AdaLasso_scenario
  SCAD_everything[[c]] = SCAD_scenario
  Tilting_everything[[c]] = Tilting_scenario
  PC_everything[[c]] = PC_scenario
  
}

#####
indi = n.start  + stepsize * ((1:Scenarios)-1)
Tilting_indi = indi[indi<=120]
Tilting_summary = Tilting_summary[,as.character(Tilting_indi),drop=F]
#####


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
#PF_holds_freq
#OIP_holds_freq
#OIP_prime_holds_freq
#Thresholded_subseteq_best_freq
#Thresholded_equals_best_freq


##### Save ############################################################
Save1 = list (   Forward_summary = Forward_summary,
                 AdaSub_summary = AdaSub_summary,
                 AdaSubBest_summary = AdaSubBest_summary,
                 Stability_summary = Stability_summary,
                 Lasso_summary = Lasso_summary,
                 AdaLasso_summary = AdaLasso_summary,
                 SCAD_summary = SCAD_summary,
                 Tilting_summary = Tilting_summary,
                 PC_summary = PC_summary,
                 Forward_everything = Forward_everything,
                 AdaSub_everything = AdaSub_everything,
                 AdaSubBest_everything = AdaSubBest_everything,
                 Lasso_everything = Lasso_everything,
                 Stability_everything = Stability_everything,
                 AdaLasso_everything = AdaLasso_everything,
                 SCAD_everything = SCAD_everything,
                 Tilting_everything = Tilting_everything,
                 PC_everything = PC_everything,
                 p=p, n_test = n_test, Iterations = Iterations, q=q, corr = corr, const = const, 
                 N_Ex = N_Ex, Scenarios = Scenarios, stepsize = stepsize , n.start = n.start,
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
                 OIP_prime_count_rel = OIP_prime_count_rel
                 )


#setwd("C://Users//Staerk//scieboBonn//AdaSub paper for EJS//R simulation saves")
#save(Save1, file="All_info_EBIC06_Global00_p_increasing_400_2000_simulations500.RData")

#########################################################################

N_Ex = Save1$N_Ex

ploteq =FALSE

PC = FALSE

###########################################################################################

# Color Legend: 
# BLACK: AdaSub Thresholded
# ORANGE: AdaSub Best
# RED: Stability Selection
# BLUE: Lasso
# GREEN: Stepwise
# GRAY: Adaptive Lasso 
# PURPLE: SCAD
# PINK: Tilting
# YELLOW: PC-Simple


color_AdaSub = "black"
color_AdaSubBest = "orange"
color_Stability = "red"
color_Forward = "green"
color_Lasso = "blue"
color_AdaLasso = "gray"
color_SCAD = "purple"
color_Tilting = "pink"
color_PC = "yellow"

legendnames = c("AdaSubThres", "AdaSubBest", "StabSel", "Forward", "Lasso", "AdaLasso", "SCAD", "Tilting")

pchs = c(4,8,1,2,3,5,6,0)
cex.size = 1.5
ltys = c(1,2,3,4,5,6,7,8)

dev.off()
win.graph(width = 14, height = 8, pointsize = 8)

font = 2
par(las=1,lwd=2,cex.main=font, cex.lab=font, cex.axis=font)
par(mfcol=c(2,3))

xaxis = Save1$n.start  + Save1$stepsize * ((1:Save1$Scenarios)-1)
xaxis.Tilting = Tilting_indi

############################################################################
#m <- matrix(c(1,2,3,4,5,3,6,7,3),nrow = 3,ncol = 3,byrow = FALSE)
#layout(mat = m, heights = c(0.45,0.45,0.1))
############################################################################

m <- matrix(c(1,2,3,4,1,5,6,4,1,7,8,4),nrow = 4,ncol = 3,byrow = FALSE)
layout(mat = m, heights = c(0.05,0.425,0.425,0.1))


par(mar=c(0, 0, 2, 0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="", cex.main=1.5,
     main=paste("High-dimensional simulations (p=10n) for EBIC with gamma =",const,"with",N_Ex, "simulated datasets for", Scenarios, "scenarios"))
par(mar=c(5.1, 4.1, 4.1, 2.1))


values_false_pos = range(c(0,5,
                           Save1$Forward_summary["false_pos_mean",],
                           Save1$Lasso_summary["false_pos_mean",],
                           Save1$AdaLasso_summary["false_pos_mean",],
                           Save1$SCAD_summary["false_pos_mean",],
                           Save1$Stability_summary["false_pos_mean",], 
                           Save1$AdaSub_summary["false_pos_mean",],
                           Save1$Tilting_summary["false_pos_mean",],
                           Save1$PC_summary["false_pos_mean",]))


plot(xaxis,Save1$AdaSub_summary["false_pos_mean",],pch=pchs[1],cex=cex.size,ylim=c(0,5), xlab="n",ylab="",main="Mean false positives",xaxt="n",col=color_AdaSub, type="b", lty = ltys[1])
points(xaxis,Save1$AdaSubBest_summary["false_pos_mean",],pch=pchs[2],cex=cex.size,col=color_AdaSubBest, type="b", lty = ltys[2])
points(xaxis,Save1$Stability_summary["false_pos_mean",],pch=pchs[3],cex=cex.size,col=color_Stability, type="b", lty = ltys[3])
points(xaxis,Save1$Forward_summary["false_pos_mean",],pch=pchs[4],cex=cex.size,col=color_Forward, type="b", lty = ltys[4])
points(xaxis,Save1$Lasso_summary["false_pos_mean",],pch=pchs[5],cex=cex.size,col=color_Lasso, type="b", lty = ltys[5])
points(xaxis,Save1$AdaLasso_summary["false_pos_mean",],pch=pchs[6],cex=cex.size,col=color_AdaLasso, type="b", lty = ltys[6])
points(xaxis,Save1$SCAD_summary["false_pos_mean",],pch=pchs[7],cex=cex.size,col=color_SCAD, type="b", lty = ltys[7])
points(xaxis.Tilting,Save1$Tilting_summary["false_pos_mean",],pch=pchs[8],cex=cex.size,col=color_Tilting, type="b", lty = ltys[8])
if (PC){
points(xaxis,Save1$PC_summary["false_pos_mean",],pch=20,cex=2.5,col=color_PC, type="b", lty = ltys[1])
}
axis(1,at=xaxis)
points(xaxis,Save1$AdaSubBest_summary["false_pos_mean",],pch=pchs[2],cex=cex.size,col=color_AdaSubBest, type="b", lty = ltys[2])
points(xaxis,Save1$AdaSub_summary["false_pos_mean",],pch=pchs[1],cex=cex.size,col=color_AdaSub, type="b", lty = ltys[1])



plot(xaxis,Save1$AdaSub_summary["false_neg_mean",],pch=pchs[1],cex=cex.size,ylim=c(0,5), xlab="n",ylab="",main="Mean false negatives",xaxt="n",col=color_AdaSub, type="b", lty = ltys[1])
points(xaxis,Save1$AdaSubBest_summary["false_neg_mean",],pch=pchs[2],cex=cex.size,col=color_AdaSubBest, type="b", lty = ltys[2])
points(xaxis,Save1$Stability_summary["false_neg_mean",],pch=pchs[3],cex=cex.size,col=color_Stability, type="b", lty = ltys[3])
points(xaxis,Save1$Forward_summary["false_neg_mean",],pch=pchs[4],cex=cex.size,col=color_Forward, type="b", lty = ltys[4])
points(xaxis,Save1$Lasso_summary["false_neg_mean",],pch=pchs[5],cex=cex.size,col=color_Lasso, type="b", lty = ltys[5])
points(xaxis,Save1$AdaLasso_summary["false_neg_mean",],pch=pchs[6],cex=cex.size,col=color_AdaLasso, type="b", lty = ltys[6])
points(xaxis,Save1$SCAD_summary["false_neg_mean",],pch=pchs[7],cex=cex.size,col=color_SCAD, type="b", lty = ltys[7])
points(xaxis.Tilting,Save1$Tilting_summary["false_neg_mean",],pch=pchs[8],cex=cex.size,col=color_Tilting, type="b", lty = ltys[8])
if (PC){
points(xaxis,Save1$PC_summary["false_neg_mean",],pch=pchs[1],cex=cex.size,col=color_PC)
}
axis(1,at=xaxis)
points(xaxis,Save1$AdaSubBest_summary["false_neg_mean",],pch=pchs[2],cex=cex.size,col=color_AdaSubBest, type="b", lty = ltys[2])
points(xaxis,Save1$AdaSub_summary["false_neg_mean",],pch=pchs[1],cex=cex.size,col=color_AdaSub, type="b", lty = ltys[1])




legend_order <- matrix(1:8,ncol=4,byrow = TRUE)
par(mar=c(0,0,0,0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
plot_colors <- c(color_AdaSub, color_AdaSubBest, color_Stability, color_Forward, color_Lasso, color_AdaLasso, color_SCAD, color_Tilting)
legend(x = "center",inset = 0,
       legend = legendnames[legend_order], 
       col=plot_colors[legend_order], cex=2, box.lty=1, ncol = 4,
       pch=pchs[legend_order], pt.bg = 'white', lty=ltys[legend_order]) #,
par(mar=c(5.1, 4.1, 4.1, 2.1))




plot(xaxis,Save1$AdaSub_summary["equal_count_0",]/N_Ex,pch=pchs[1],cex=cex.size,ylim=c(0,0.5), xlab="n",ylab="",main="Rel. freq. true model selected",xaxt="n",col=color_AdaSub, type="b", lty = ltys[1])
axis(1,at=xaxis)
points(xaxis,Save1$AdaSubBest_summary["equal_count_0",]/N_Ex,pch=pchs[2],cex=cex.size,col=color_AdaSubBest, type="b", lty = ltys[2])
points(xaxis,Save1$Stability_summary["equal_count_0",]/N_Ex,pch=pchs[3],cex=cex.size,col=color_Stability, type="b", lty = ltys[3])
points(xaxis,Save1$Forward_summary["equal_count_0",]/N_Ex,pch=pchs[4],cex=cex.size,col=color_Forward, type="b", lty = ltys[4])
points(xaxis,Save1$Lasso_summary["equal_count_0",]/N_Ex,pch=pchs[5],cex=cex.size,col=color_Lasso, type="b", lty = ltys[5])
points(xaxis,Save1$AdaLasso_summary["equal_count_0",]/N_Ex,pch=pchs[6],cex=cex.size,col=color_AdaLasso, type="b", lty = ltys[6])
points(xaxis,Save1$SCAD_summary["equal_count_0",]/N_Ex,pch=pchs[7],cex=cex.size,col=color_SCAD, type="b", lty = ltys[7])
points(xaxis.Tilting,Save1$Tilting_summary["equal_count_0",]/N_Ex,pch=pchs[8],cex=cex.size,col=color_Tilting, type="b", lty = ltys[8])
if (PC){
points(xaxis,Save1$PC_summary["equal_count_0",]/N_Ex,pch=20,cex=2.5,col=color_PC)
}
points(xaxis,Save1$AdaSubBest_summary["equal_count_0",]/N_Ex,pch=pchs[2],cex=cex.size,col=color_AdaSubBest, type="b", lty = ltys[2])
points(xaxis,Save1$AdaSub_summary["equal_count_0",]/N_Ex,pch=pchs[1],cex=cex.size,col=color_AdaSub, type="b", lty = ltys[1])



#### Plot computational times ##############
values_times = range(
  Save1$Forward_summary["times",],
  Save1$Lasso_summary["times",],
  Save1$AdaLasso_summary["times",],
  Save1$SCAD_summary["times",],
  Save1$Stability_summary["times",], 
  Save1$AdaSub_summary["times",],
  Save1$Tilting_summary["times",],
  Save1$PC_summary["times",]
)

plot(xaxis,Save1$AdaSub_summary["times",],ylim=values_times,
     pch=pchs[1],cex=cex.size,xaxt="n",xlab="n",ylab="",main="Mean comp. time (in s)",col=color_AdaSub, type="b", lty = ltys[1])
lines(xaxis,Save1$AdaSub_summary["times",],col=color_AdaSubBest, lty = ltys[2])
points(xaxis,Save1$Stability_summary["times",],pch=pchs[3],cex=cex.size,col=color_Stability, type="b", lty = ltys[3])
points(xaxis,Save1$Forward_summary["times",],pch=pchs[4],cex=cex.size,col=color_Forward, type="b", lty = ltys[4])
points(xaxis,Save1$Lasso_summary["times",],pch=pchs[5],cex=cex.size,col=color_Lasso, type="b", lty = ltys[5])
points(xaxis,Save1$AdaLasso_summary["times",],pch=pchs[6],cex=cex.size,col=color_AdaLasso, type="b", lty = ltys[6])
points(xaxis,Save1$SCAD_summary["times",],pch=pchs[7],cex=cex.size,col=color_SCAD, type="b", lty = ltys[7])
points(xaxis.Tilting,Save1$Tilting_summary["times",],pch=pchs[8],cex=cex.size,col=color_Tilting, type="b", lty = ltys[8])
if (PC){
points(xaxis,Save1$PC_summary["times",],pch=20,cex=cex.size,col=color_PC)
}
axis(1,at=xaxis)



values_estimation = range(
  Save1$Forward_summary["estimation_error_mean",],
  Save1$Lasso_summary["estimation_error_mean",],
  Save1$AdaLasso_summary["estimation_error_mean",],
  Save1$SCAD_summary["estimation_error_mean",],
  Save1$Stability_summary["estimation_error_mean",], 
  Save1$AdaSub_summary["estimation_error_mean",],
  Save1$AdaSubBest_summary["estimation_error_mean",],
  Save1$Tilting_summary["estimation_error_mean",],
  Save1$PC_summary["estimation_error_mean",]
)

plot(xaxis,Save1$AdaSub_summary["estimation_error_mean",],ylim=values_estimation,
     pch=pchs[1],cex=cex.size,xaxt="n",xlab="n",ylab="",main="Mean estimation error (MSE)",col=color_AdaSub, type="b", lty = ltys[1])
points(xaxis,Save1$AdaSubBest_summary["estimation_error_mean",],pch=pchs[2],cex=cex.size,col=color_AdaSubBest, type="b", lty = ltys[2])
points(xaxis,Save1$Stability_summary["estimation_error_mean",],pch=pchs[3],cex=cex.size,col=color_Stability, type="b", lty = ltys[3])
points(xaxis,Save1$Forward_summary["estimation_error_mean",],pch=pchs[4],cex=cex.size,col=color_Forward, type="b", lty = ltys[4])
points(xaxis,Save1$Lasso_summary["estimation_error_mean",],pch=pchs[5],cex=cex.size,col=color_Lasso, type="b", lty = ltys[5])
points(xaxis,Save1$AdaLasso_summary["estimation_error_mean",],pch=pchs[6],cex=cex.size,col=color_AdaLasso, type="b", lty = ltys[6])
points(xaxis,Save1$SCAD_summary["estimation_error_mean",],pch=pchs[7],cex=cex.size,col=color_SCAD, type="b", lty = ltys[7])
points(xaxis.Tilting,Save1$Tilting_summary["estimation_error_mean",],pch=pchs[8],cex=cex.size,col=color_Tilting, type="b", lty = ltys[8])
if (PC){
points(xaxis,Save1$PC_summary["estimation_error_mean",],pch=20,cex=cex.size,col=color_PC)
}
axis(1,at=xaxis)
points(xaxis,Save1$AdaSubBest_summary["estimation_error_mean",],pch=pchs[2],cex=cex.size,col=color_AdaSubBest, type="b", lty = ltys[2])
points(xaxis,Save1$AdaSub_summary["estimation_error_mean",],pch=pchs[1],cex=cex.size,col=color_AdaSub, type="b", lty = ltys[1])





values_prediction = range( Save1$Forward_summary["prediction_error_mean",],
                           Save1$Lasso_summary["prediction_error_mean",],
                           Save1$AdaLasso_summary["prediction_error_mean",],
                           Save1$SCAD_summary["prediction_error_mean",],
                           Save1$Stability_summary["prediction_error_mean",], 
                           Save1$AdaSub_summary["prediction_error_mean",],
                           Save1$AdaSubBest_summary["prediction_error_mean",],
                           Save1$Tilting_summary["prediction_error_mean",] #,
                        #   Save1$PC_summary["prediction_error_mean",]
)

plot(xaxis,Save1$AdaSub_summary["prediction_error_mean",],ylim=values_prediction,
     pch=pchs[1],cex=cex.size,xaxt="n",xlab="n",ylab="",main="Mean prediction error (RMSE)",col=color_AdaSub, type="b", lty = ltys[1])
points(xaxis,Save1$AdaSubBest_summary["prediction_error_mean",],pch=pchs[2],cex=cex.size,col=color_AdaSubBest, type="b", lty = ltys[2])
points(xaxis,Save1$Stability_summary["prediction_error_mean",],pch=pchs[3],cex=cex.size,col=color_Stability, type="b", lty = ltys[3])
points(xaxis,Save1$Forward_summary["prediction_error_mean",],pch=pchs[4],cex=cex.size,col=color_Forward, type="b", lty = ltys[4])
points(xaxis,Save1$Lasso_summary["prediction_error_mean",],pch=pchs[5],cex=cex.size,col=color_Lasso, type="b", lty = ltys[5])
points(xaxis,Save1$AdaLasso_summary["prediction_error_mean",],pch=pchs[6],cex=cex.size,col=color_AdaLasso, type="b", lty = ltys[6])
points(xaxis,Save1$SCAD_summary["prediction_error_mean",],pch=pchs[7],cex=cex.size,col=color_SCAD, type="b", lty = ltys[7])
points(xaxis.Tilting,Save1$Tilting_summary["prediction_error_mean",],pch=pchs[8],cex=cex.size,col=color_Tilting, type="b", lty = ltys[8])
if (PC){
points(xaxis,Save1$PC_summary["prediction_error_mean",],pch=20,cex=cex.size,col=color_PC)
}
axis(1,at=xaxis)
points(xaxis,Save1$AdaSubBest_summary["prediction_error_mean",],pch=pchs[2],cex=cex.size,col=color_AdaSubBest, type="b", lty = ltys[2])
points(xaxis,Save1$AdaSub_summary["prediction_error_mean",],pch=pchs[1],cex=cex.size,col=color_AdaSub, type="b", lty = ltys[1])


par(las=1,lwd=1,cex.main=1, cex.lab=1, cex.axis=1)



##############################################












### Reproduce results for high-dimensional simulation study (see Section 5.2 of paper) ############################

# Load main functions (Please adjust to the right directory!)
# setwd("C://Users//Staerk//scieboBonn//AdaSub paper for EJS//R files//files") 
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


choices <-c( "500 examples for each of 9 different sample sizes (as in paper)", "Custom number of examples per sample size")
mychoice <- menu( choices , graphics=TRUE, title="How many examples per sample size should be considered?" )

if (mychoice==1) {
  N_Ex=500 #100 
  Scenarios = 9
}
if (mychoice==2) {
  N_Ex = -1
  while(N_Ex < 1 ){
    N_Ex <- readline("Enter the number of examples per sample size: ")
    N_Ex <- ifelse(grepl("\\D",N_Ex),-1,as.integer(N_Ex))
    if(is.na(N_Ex)){break}  # breaks when hit enter
  }
  Scenarios = 9
}


n_test=100 # Sample size for independent test set (->prediction error)

#N_Ex=500 # number of replicates per scenario
#Scenarios = 9 # number of scenarios (different sample sizes n)
stepsize = 20 # stepsize for increasing sample size n
n.start = 40 # sample size for first scenario (note that number of variables is given by p=10n)

LassoCV_summary = method_init(Scenarios,n.start=n.start,stepsize=stepsize)
AdaLassoCV_summary = method_init(Scenarios,n.start=n.start,stepsize=stepsize)
SCAD_CV_summary = method_init(Scenarios,n.start=n.start,stepsize=stepsize)

### Lists which include all results for all scenarios (first entry, matrix of results for first scenario)
LassoCV_everything = list()
AdaLassoCV_everything = list()
SCAD_CV_everything = list()


for (c in 1:Scenarios){
  
  n= n.start + stepsize*(c-1)
  p = 10 * n # Growing p
  #p = 1/4*n^2
  
  LassoCV_scenario = method_init_scenario(N_Ex)
  AdaLassoCV_scenario = method_init_scenario(N_Ex)
  SCAD_CV_scenario = method_init_scenario(N_Ex)
  
  
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
    
    ## Lasso with cross-validation (CV)
    #seedsafe <- .Random.seed
    start.time <- Sys.time()
    cv.fit=cv.glmnet(data$x,data$y,intercept=TRUE,family="gaussian",standardize=FALSE,alpha=1)
    betas.cv.fit=as.vector(predict(cv.fit,type="coefficients", s="lambda.1se"))
    LassoCV_nonzero = which(betas.cv.fit[-1]!=0)
    end.time <- Sys.time()
    time = as.double(end.time-start.time,units = "secs")
  
    LassoCV_scenario[,k] = method_results_scenario(model=LassoCV_nonzero,beta1=beta1,s0=s0,data.test=data.test,n_test=n_test,time=time,beta=betas.cv.fit,data=data)
    if (setequal(LassoCV_nonzero,S0))  LassoCV_summary["equal_count_0",c] = LassoCV_summary["equal_count_0",c] + 1
    
    #set.seed(seedsafe)
    
    ## Adaptive Lasso with CV in first step and CV in second step
    start.time <- Sys.time()
    cv.fit1=cv.glmnet(data$x,data$y,intercept=TRUE,family="gaussian",standardize=FALSE,alpha=1)
    betas.cv.fit1=as.vector(predict(cv.fit1,type="coefficients",s="lambda.1se"))[-1]
    new_penalty=rep(0,p)
    new_penalty=1/pmax(abs(betas.cv.fit1),.Machine$double.eps)
    cv.fit2=cv.glmnet(data$x,data$y,penalty.factor=new_penalty,intercept=TRUE,family="gaussian",standardize=FALSE,alpha=1)
    AdaLassoCV_beta=as.vector(predict(cv.fit2,type="coefficients",s="lambda.1se"))
    AdaLassoCV_nonzero = which(AdaLassoCV_beta[-1]!=0)
    if (is.null(AdaLassoCV_nonzero)) AdaLassoCV_nonzero = integer(0)
    end.time <- Sys.time()
    time = as.double(end.time-start.time,units = "secs")
    
    AdaLassoCV_scenario[,k] = method_results_scenario(model=AdaLassoCV_nonzero,beta1=beta1,s0=s0,data.test=data.test,n_test=n_test,time=time,beta=AdaLassoCV_beta,data=data)
    if (setequal(AdaLassoCV_nonzero,S0))  AdaLassoCV_summary["equal_count_0",c] = AdaLassoCV_summary["equal_count_0",c] + 1
    
    
    ## SCAD wich cross-validation (CV)
    start.time <- Sys.time()
    cvfit <- cv.ncvreg(data$x, data$y, family="gaussian", penalty="SCAD", alpha=1 ) 
    fit <- cvfit$fit
    min1se <- min(which(cvfit$cve <= cvfit$cve[cvfit$min] + cvfit$cvse[cvfit$min]))
    SCAD_CV_beta <- fit$beta[,min1se] 
    SCAD_CV_nonzero = which(SCAD_CV_beta[-1]!= 0)
    names(SCAD_CV_nonzero) =NULL
    if (is.null(SCAD_CV_nonzero)) SCAD_CV_nonzero = integer(0)
    end.time <- Sys.time()
    time = as.double(end.time-start.time,units = "secs")
    
    SCAD_CV_scenario[,k] = method_results_scenario(model=SCAD_CV_nonzero,beta1=beta1,s0=s0,data.test=data.test,n_test=n_test,time=time,beta=SCAD_CV_beta,data=data) 
    if (setequal(SCAD_CV_nonzero,S0))  SCAD_CV_summary["equal_count_0",c] = SCAD_CV_summary["equal_count_0",c] + 1
    
    #######################################################################
    
    cat("n = ", n,", Example ", k, "\n")
    print("True model:")
    print(S0)
    print("Lasso CV model:")
    print(LassoCV_nonzero)
    print("AdaLasso CV model:")
    print(AdaLassoCV_nonzero)
    print("SCAD CV model:")
    print(SCAD_CV_nonzero)
    print("#####################################")
  }

  LassoCV_summary[1:5,c] = 1/N_Ex*rowSums(LassoCV_scenario)
  AdaLassoCV_summary[1:5,c] = 1/N_Ex*rowSums(AdaLassoCV_scenario)
  SCAD_CV_summary[1:5,c] = 1/N_Ex*rowSums(SCAD_CV_scenario)
  
  LassoCV_everything[[c]] = LassoCV_scenario
  AdaLassoCV_everything[[c]] = AdaLassoCV_scenario
  SCAD_CV_everything[[c]] = SCAD_CV_scenario
  
}


##### Save ############################################################
Save_CV = list ( LassoCV_summary = LassoCV_summary,
                 AdaLassoCV_summary = AdaLassoCV_summary,
                 SCAD_CV_summary = SCAD_CV_summary,
                 LassoCV_everything = LassoCV_everything,
                 AdaLassoCV_everything = AdaLassoCV_everything,
                 SCAD_CV_everything = SCAD_CV_everything,
                 p=p, n_test = n_test, corr = corr, 
                 N_Ex = N_Ex, Scenarios = Scenarios, stepsize = stepsize , n.start = n.start)

#setwd("C://Users//Staerk//scieboBonn//AdaSub paper for EJS//R simulation saves")
#save(Save_CV, file="CV_Blocks10_05_p_increasing_400_2000_simulations500.RData")
#########################################################################

N_Ex = Save1$N_Ex

ploteq =FALSE

PC = FALSE

#####
if (corr.type =="global" & corr==0) {
  load("./Files/R_simulation_saves/EBIC06_Global00_p_increasing_400_2000_simulations500.RData")
  Save06 = Save1
  load("./Files/R_simulation_saves/EBIC1_Global00_p_increasing_400_2000_simulations500.RData")
  #load("./Files/R_simulation_saves/CV_Global00_p_increasing_400_2000_simulations500.RData")
}
if (corr.type =="global" & corr==0.7) {
  load("./Files/R_simulation_saves/EBIC06_Global07_p_increasing_400_2000_simulations500.RData")
  Save06 = Save1
  load("./Files/R_simulation_saves/EBIC1_Global07_p_increasing_400_2000_simulations500.RData")
  #load("./Files/R_simulation_saves/CV_Global07_p_increasing_400_2000_simulations500.RData")
}
if (corr.type =="toeplitz" & corr==0.9) {
  load("./Files/R_simulation_saves/EBIC06_Toeplitz09_p_increasing_400_2000_simulations500.RData")
  Save06 = Save1
  load("./Files/R_simulation_saves/EBIC1_Toeplitz09_p_increasing_400_2000_simulations500.RData")
  #load("./Files/R_simulation_saves/CV_Toeplitz09_p_increasing_400_2000_simulations500.RData")
}
if (corr.type =="block" & corr==0.5) {
  load("./Files/R_simulation_saves/EBIC06_Blocks10_05_p_increasing_400_2000_simulations500.RData")
  Save06 = Save1
  load("./Files/R_simulation_saves/EBIC1_Blocks10_05_p_increasing_400_2000_simulations500.RData")
  #load("./Files/R_simulation_saves/CV_Blocks10_05_p_increasing_400_2000_simulations500.RData")
}
#####

###########################################################################################
### Plots of CV-results (in comparison to AdaSub) ########################################
###########################################################################################

# Color Legend: 
# BLACK: AdaSub Thresholded
# ORANGE: AdaSub Best
# BLUE: Lasso CV
# GRAY: Adaptive Lasso CV
# PURPLE: SCAD CV

colors = c("black","gray40","orange","darkorange2")


color_AdaSub0.6 = "black"
color_AdaSubBest0.6 = "orange"
color_AdaSub1 = "gray40"
color_AdaSubBest1 = "darkorange2"
color_Lasso = "blue"
color_AdaLasso = "gray"
color_SCAD = "purple"

legendnames = c(expression(paste("AdaSubThres ", gamma," = 0.6  ")), 
                expression(paste("AdaSubBest ", gamma," = 0.6  ")),
                expression(paste("AdaSubThres ", gamma," = 1  ")),
                expression(paste("AdaSubBest ", gamma," = 1  ")),
                "Lasso CV", "AdaLasso CV", "SCAD CV")

pchs = c(4,1,8,2,3,5,6)
cex.size = 1.5
ltys = c(1,3,2,4,5,6,7)

#dev.off()
win.graph(width = 14, height = 8, pointsize = 8)

font = 2
par(las=1,lwd=2,cex.main=font, cex.lab=font, cex.axis=font)
par(mfcol=c(2,3))

xaxis = Save_CV$n.start  + Save_CV$stepsize * ((1:Save_CV$Scenarios)-1)

############################################################################
m <- matrix(c(1,2,3,4,5,3,6,7,3),nrow = 3,ncol = 3,byrow = FALSE)
layout(mat = m, heights = c(0.45,0.45,0.1))
############################################################################

#m <- matrix(c(1,2,3,4,1,5,6,4,1,7,8,4),nrow = 4,ncol = 3,byrow = FALSE)
#layout(mat = m, heights = c(0.05,0.425,0.425,0.1))


#par(mar=c(0, 0, 2, 0))
#plot(1, type = "n", axes=FALSE, xlab="", ylab="", cex.main=1.5,
#     main=paste("High-dimensional simulations (p=10n) for cross-validated models with",N_Ex, "simulated datasets for", Scenarios, "scenarios"))
par(mar=c(5.1, 4.1, 4.1, 2.1))

values_false_pos = range(c(0,5,
  Save06$AdaSub_summary["false_pos_mean",],
  Save1$AdaSub_summary["false_pos_mean",],
  Save_CV$LassoCV_summary["false_pos_mean",],
  Save_CV$AdaLassoCV_summary["false_pos_mean",],
  Save_CV$SCAD_CV_summary["false_pos_mean",])
)


plot(xaxis,Save06$AdaSub_summary["false_pos_mean",],pch=pchs[1],cex=cex.size,ylim=values_false_pos, xlab="n",ylab="",main="Mean false positives",xaxt="n",col=color_AdaSub0.6, type="b", lty = ltys[1])
points(xaxis,Save_CV$LassoCV_summary["false_pos_mean",],pch=pchs[5],cex=cex.size,col=color_Lasso, type="b", lty = ltys[5])
points(xaxis,Save_CV$SCAD_CV_summary["false_pos_mean",],pch=pchs[7],cex=cex.size,col=color_SCAD, type="b", lty = ltys[7])
points(xaxis,Save_CV$AdaLassoCV_summary["false_pos_mean",],pch=pchs[6],cex=cex.size,col=color_AdaLasso, type="b", lty = ltys[6])
points(xaxis,Save06$AdaSubBest_summary["false_pos_mean",],pch=pchs[2],cex=cex.size,col=color_AdaSubBest0.6, type="b", lty = ltys[2])
points(xaxis,Save1$AdaSub_summary["false_pos_mean",],pch=pchs[3],cex=cex.size,col=color_AdaSub1, type="b", lty = ltys[3])
points(xaxis,Save1$AdaSubBest_summary["false_pos_mean",],pch=pchs[4],cex=cex.size,col=color_AdaSubBest1, type="b", lty = ltys[4])
axis(1,at=xaxis)
#points(xaxis,Save1$AdaSubBest_summary["false_pos_mean",],pch=pchs[2],cex=cex.size,col=color_AdaSubBest, type="b", lty = ltys[2])
#points(xaxis,Save1$AdaSub_summary["false_pos_mean",],pch=pchs[1],cex=cex.size,col=color_AdaSub, type="b", lty = ltys[1])

plot(xaxis,Save06$AdaSub_summary["false_neg_mean",],pch=pchs[1],cex=cex.size,ylim=c(0,5), xlab="n",ylab="",main="Mean false negatives",xaxt="n",col=color_AdaSub0.6, type="b", lty = ltys[1])
points(xaxis,Save_CV$LassoCV_summary["false_neg_mean",],pch=pchs[5],cex=cex.size,col=color_Lasso, type="b", lty = ltys[5])
points(xaxis,Save_CV$SCAD_CV_summary["false_neg_mean",],pch=pchs[7],cex=cex.size,col=color_SCAD, type="b", lty = ltys[7])
points(xaxis,Save_CV$AdaLassoCV_summary["false_neg_mean",],pch=pchs[6],cex=cex.size,col=color_AdaLasso, type="b", lty = ltys[6])
points(xaxis,Save06$AdaSubBest_summary["false_neg_mean",],pch=pchs[2],cex=cex.size,col=color_AdaSubBest0.6, type="b", lty = ltys[2])
points(xaxis,Save1$AdaSub_summary["false_neg_mean",],pch=pchs[3],cex=cex.size,col=color_AdaSub1, type="b", lty = ltys[3])
points(xaxis,Save1$AdaSubBest_summary["false_neg_mean",],pch=pchs[4],cex=cex.size,col=color_AdaSubBest1, type="b", lty = ltys[4])
axis(1,at=xaxis)
#points(xaxis,Save1$AdaSubBest_summary["false_neg_mean",],pch=pchs[2],cex=cex.size,col=color_AdaSubBest, type="b", lty = ltys[2])
#points(xaxis,Save1$AdaSub_summary["false_neg_mean",],pch=pchs[1],cex=cex.size,col=color_AdaSub, type="b", lty = ltys[1])


legend_order <- matrix(1:8,ncol=4,byrow = TRUE)
par(mar=c(0,0,0,0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
plot_colors <- c(color_AdaSub0.6, color_AdaSubBest0.6, color_AdaSub1, color_AdaSubBest1, color_Lasso, color_AdaLasso, color_SCAD)
legend(x = "center",inset = 0,
       legend = legendnames[legend_order], 
       col=plot_colors[legend_order], cex=2, box.lty=1, ncol = 4,
       pch=pchs[legend_order], pt.bg = 'white', lty=ltys[legend_order]) #,
par(mar=c(5.1, 4.1, 4.1, 2.1))



plot(xaxis,Save06$AdaSub_summary["equal_count_0",]/Save06$N_Ex,pch=pchs[1],cex=cex.size,ylim=c(0,0.5), xlab="n",ylab="",main="Rel. freq. true model selected",xaxt="n",col=color_AdaSub0.6, type="b", lty = ltys[1])
points(xaxis,Save_CV$LassoCV_summary["equal_count_0",]/Save_CV$N_Ex,pch=pchs[5],cex=cex.size,col=color_Lasso, type="b", lty = ltys[5])
points(xaxis,Save_CV$SCAD_CV_summary["equal_count_0",]/Save_CV$N_Ex,pch=pchs[7],cex=cex.size,col=color_SCAD, type="b", lty = ltys[7])
points(xaxis,Save_CV$AdaLassoCV_summary["equal_count_0",]/Save_CV$N_Ex,pch=pchs[6],cex=cex.size,col=color_AdaLasso, type="b", lty = ltys[6])
points(xaxis,Save06$AdaSubBest_summary["equal_count_0",]/Save06$N_Ex,pch=pchs[2],cex=cex.size,col=color_AdaSubBest0.6, type="b", lty = ltys[2])
points(xaxis,Save1$AdaSub_summary["equal_count_0",]/Save1$N_Ex,pch=pchs[3],cex=cex.size,col=color_AdaSub1, type="b", lty = ltys[3])
points(xaxis,Save1$AdaSubBest_summary["equal_count_0",]/Save1$N_Ex,pch=pchs[4],cex=cex.size,col=color_AdaSubBest1, type="b", lty = ltys[4])
axis(1,at=xaxis)

#points(xaxis,Save1$AdaSubBest_summary["equal_count_0",]/N_Ex,pch=pchs[2],cex=cex.size,col=color_AdaSubBest, type="b", lty = ltys[2])
#points(xaxis,Save1$AdaSub_summary["equal_count_0",]/N_Ex,pch=pchs[1],cex=cex.size,col=color_AdaSub, type="b", lty = ltys[1])



#### Plot computational times ##############
values_times = range(
  Save06$AdaSub_summary["times",],
  Save1$AdaSub_summary["times",],
  Save_CV$LassoCV_summary["times",],
  Save_CV$AdaLassoCV_summary["times",],
  Save_CV$SCAD_CV_summary["times",]
)

plot(xaxis,Save06$AdaSub_summary["times",],pch=pchs[1],cex=cex.size,ylim=values_times, xlab="n",ylab="",main="Mean comp. time (in s)",xaxt="n",col=color_AdaSub0.6, type="b", lty = ltys[1])
lines(xaxis,Save06$AdaSub_summary["times",],col=color_AdaSubBest0.6, lty = ltys[2])
points(xaxis,Save1$AdaSub_summary["times",],pch=pchs[3],cex=cex.size,col=color_AdaSub1, type="b", lty = ltys[3])
lines(xaxis,Save1$AdaSub_summary["times",],col=color_AdaSubBest1, lty = ltys[4])
points(xaxis,Save_CV$LassoCV_summary["times",],pch=pchs[5],cex=cex.size,col=color_Lasso, type="b", lty = ltys[5])
points(xaxis,Save_CV$SCAD_CV_summary["times",],pch=pchs[7],cex=cex.size,col=color_SCAD, type="b", lty = ltys[7])
points(xaxis,Save_CV$AdaLassoCV_summary["times",],pch=pchs[6],cex=cex.size,col=color_AdaLasso, type="b", lty = ltys[6])
axis(1,at=xaxis)



values_estimation = range(
  Save_CV$LassoCV_summary["estimation_error_mean",],
  Save_CV$AdaLassoCV_summary["estimation_error_mean",],
  Save_CV$SCAD_CV_summary["estimation_error_mean",],
  Save1$AdaSub_summary["estimation_error_mean",],
  Save1$AdaSubBest_summary["estimation_error_mean",],
  Save06$AdaSub_summary["estimation_error_mean",],
  Save06$AdaSubBest_summary["estimation_error_mean",]
)

plot(xaxis,Save06$AdaSub_summary["estimation_error_mean",],pch=pchs[1],cex=cex.size,ylim=values_estimation, xlab="n",ylab="",main="Mean estimation error (MSE)",xaxt="n",col=color_AdaSub0.6, type="b", lty = ltys[1])
points(xaxis,Save_CV$LassoCV_summary["estimation_error_mean",],pch=pchs[5],cex=cex.size,col=color_Lasso, type="b", lty = ltys[5])
points(xaxis,Save_CV$SCAD_CV_summary["estimation_error_mean",],pch=pchs[7],cex=cex.size,col=color_SCAD, type="b", lty = ltys[7])
points(xaxis,Save_CV$AdaLassoCV_summary["estimation_error_mean",],pch=pchs[6],cex=cex.size,col=color_AdaLasso, type="b", lty = ltys[6])
points(xaxis,Save06$AdaSubBest_summary["estimation_error_mean",],pch=pchs[2],cex=cex.size,col=color_AdaSubBest0.6, type="b", lty = ltys[2])
points(xaxis,Save1$AdaSub_summary["estimation_error_mean",],pch=pchs[3],cex=cex.size,col=color_AdaSub1, type="b", lty = ltys[3])
points(xaxis,Save1$AdaSubBest_summary["estimation_error_mean",],pch=pchs[4],cex=cex.size,col=color_AdaSubBest1, type="b", lty = ltys[4])
axis(1,at=xaxis)




values_prediction = range(Save_CV$LassoCV_summary["prediction_error_mean",],
                          Save_CV$AdaLassoCV_summary["prediction_error_mean",],
                          Save_CV$SCAD_CV_summary["prediction_error_mean",],
                          Save1$AdaSub_summary["prediction_error_mean",],
                          Save1$AdaSubBest_summary["prediction_error_mean",],
                          Save06$AdaSub_summary["prediction_error_mean",],
                          Save06$AdaSubBest_summary["prediction_error_mean",]
)

plot(xaxis,Save06$AdaSub_summary["prediction_error_mean",],pch=pchs[1],cex=cex.size,ylim=values_prediction, xlab="n",ylab="",main="Mean prediction error (RMSE)",xaxt="n",col=color_AdaSub0.6, type="b", lty = ltys[1])
points(xaxis,Save_CV$LassoCV_summary["prediction_error_mean",],pch=pchs[5],cex=cex.size,col=color_Lasso, type="b", lty = ltys[5])
points(xaxis,Save_CV$SCAD_CV_summary["prediction_error_mean",],pch=pchs[7],cex=cex.size,col=color_SCAD, type="b", lty = ltys[7])
points(xaxis,Save_CV$AdaLassoCV_summary["prediction_error_mean",],pch=pchs[6],cex=cex.size,col=color_AdaLasso, type="b", lty = ltys[6])
points(xaxis,Save06$AdaSubBest_summary["prediction_error_mean",],pch=pchs[2],cex=cex.size,col=color_AdaSubBest0.6, type="b", lty = ltys[2])
points(xaxis,Save1$AdaSub_summary["prediction_error_mean",],pch=pchs[3],cex=cex.size,col=color_AdaSub1, type="b", lty = ltys[3])
points(xaxis,Save1$AdaSubBest_summary["prediction_error_mean",],pch=pchs[4],cex=cex.size,col=color_AdaSubBest1, type="b", lty = ltys[4])
axis(1,at=xaxis)

par(las=1,lwd=1,cex.main=1, cex.lab=1, cex.axis=1)









########################################################################
############################ Load Packages #############################
########################################################################

#library("leaps") 
#library("MASS")


########################################################################
############################ Load Functions ############################
########################################################################

# Compute extended BIC for model specified by "indices" with constant "const" (BIC: const=0)
# Number of parameters is corrected for estimating intercept and variance (only constant summand, not important for minimization)
# "data" should be a list with data$x as the design matrix and data$y as the response vector
EBIC <- function(data,indices,const) {  
  n=nrow(data$x)
  p=ncol(data$x)
  x.cur=data$x[,indices]
  x.cur=cbind(c(rep(1,n)),x.cur)                 #include intercept!
  lm.out = .lm.fit(x.cur, data$y)
  deviance = n*(1+log(2*pi)+ log(sum(lm.out$residuals^2) /n))
  EBIC = deviance + log(n)*(length(indices)+2) + 2*(length(indices)+2)*const*log(p)
  return(EBIC)
}


# Adaptive Subspace Method

# Input: 
# data = (data$x,data$y) ("data" is list with data$x as the design matrix and data$y as the response vector)
# Iter: number of Iterations
# K: learning rate 
# q: initial search size
#
# criterion "EBIC" (with constant "const", BIC: const=0) 
# only every "savings" iterations is saved (less memory usage)
# Upper computational bound U_C (maximal size of subsets for low-dimensional problems)

# Output: 
# relfreq.hist: relative frequency history (matrix of dimension p X floor(Iter/savings)) 
# relfreq.final: vector of final selection probabilities (after the last iteration)
# V.size, S.size: sizes of sampled models V and selected models S=f(V^(t)) in each iteration (vectors of dimension Iter)
# values: values of the Information Crition of the "best" models f(V^(t)) along iterations (vector of dimension Iter)
# best.S: "Best" model (according to the criterion) that has been found during the whole search 
# best.models: List of "best" submodels (of length Iter)

AdaSub <-function (data,Iter,K=100,q=10,criterion="EBIC",const=0,savings=10, U_C = 40, CV.null=NULL, CS.null=NULL, plot.text=NULL, plotting=TRUE) { 
  
  p=ncol(data$x)
  n=nrow(data$x)
  
  CV=matrix(numeric(Iter/savings*p),p,floor(Iter/savings))
  CS=matrix(numeric(Iter/savings*p),p,floor(Iter/savings))
  
  CV.cur = numeric(p)  
  CS.cur = numeric(p) 
  
  if (is.null(CV.null)) CV.null = numeric(p)
  if (is.null(CS.null)) CS.null = numeric(p)
  
  #relfreq=rep(q/p,p)
  relfreq = (q+CS.null+CS.cur)/(p+CV.null+CV.cur)
  
  V.size<-numeric(Iter)
  S.size<-numeric(Iter)
  
  values <- rep(NA,Iter)
  
  criterion.EBIC <-  criterion=="EBIC"
  
  best.models = vector("list", Iter)
  
  S = integer(0)

  for (t in 1:Iter) {
    
    b = rbinom(p,1,relfreq)
    V = which(b==1)
    
    if (length(V)>U_C){ V=sample(V,U_C) } 
    
    if (length(V)>1){
        regs=summary(regsubsets(as.matrix(data$x[,V]),data$y,intercept=TRUE,nvmax=min(c(length(V),n-3))))
        modelmatrix=as.matrix(regs$which[,-1])
        vec=rep(0,nrow(modelmatrix)+1)
        vec[nrow(modelmatrix)+1] = EBIC(data,NULL,const)
        for (j in 1:nrow(modelmatrix)) {
          indices.help=V[as.vector(which(modelmatrix[j,]==TRUE))]
          vec[j]=EBIC(data,indices.help,const)
        }
        mini = which.min(vec)
        if (mini<= nrow(modelmatrix)) { S=V[as.vector(which(modelmatrix[which.min(vec),]==TRUE))] } else { S = integer(0) }
        values[t] = min(vec)
        CV.cur[V] = CV.cur[V] +1
        CS.cur[S] = CS.cur[S] +1
        V.size[t]=length(V)
        S.size[t]=length(S)
        relfreq[V] = (q + CS.null[V] + K*CS.cur[V])/(p +CV.null[V] + K*CV.cur[V])
    }
    
    if (length(V)==1) { 
      model_null = EBIC(data,NULL,const)
      model_one = EBIC(data,V,const)
      if (model_null<model_one) { S=integer(0) } else { S=V }
      values[t] = min(model_null,model_one)
      CV.cur[V] = CV.cur[V] +1
      CS.cur[S] = CS.cur[S] +1
      V.size[t]=length(V)
      S.size[t]=length(S)
      relfreq[V] = (q + CS.null[V] + K*CS.cur[V])/(p + CV.null[V] + K*CV.cur[V])
    }
    
    if (length(V)==0) { values[t] = EBIC(data,NULL,const) }
    
    if (t %% savings == 0) {
      CV[,t/savings] = CV.cur
      CS[,t/savings] = CS.cur
    }
    
    if (values[t]==min(values[1:t])) best.S = S
    
    best.models[[t]] = S
    
    if (t %% 1000 == 0 && plotting) { 
      par(mfrow=c(1,1))
      plot(values,pch=20,main=paste("Dynamic plot of EBIC values for gamma =", const, "(after", t ,"iterations) \n Red line: Thresholded model (0.9)"),
           #sub="Red line: Thresholded model (0.9)",
           #font.sub=3,col.sub = "red",
           sub = plot.text ,
           xlab="Iteration",ylab="EBIC value of selected model")
      abline(h=EBIC(data,which(relfreq>0.9),const),col="red")
    }
  }
  
  #importance = CS.cur/CV.cur
  relfreq.hist=(q + CS.null + K*CS)/(p + CV.null + K*CV) 
  relfreq.final=relfreq.hist[,floor(Iter/savings)]
  return(list(relfreq.hist=relfreq.hist, relfreq.final=relfreq.final,S.size=S.size,V.size=V.size,values=values,best.S=best.S,best.models=best.models, CV.cur=CV.cur, CS.cur=CS.cur,
              CS=CS, CV=CV))
}   

# Function for checking whether (empirical) PF, OIP and OIP' hold

check_OIP_PF <- function(output, p,  S_star=NA, rho = 0.9) {
  
  if (all(is.na(S_star))) S_star = output$best.S
  
  sampled <- matrix(NA, nrow = p, ncol = dim(output$CV)[2])
  rownames(sampled) <- 1:p
  selected <- matrix(NA, nrow = p, ncol = dim(output$CV)[2])
  rownames(selected) <- 1:p
  
  for (j in 1:p) {
    sampled[as.character(j),] <-  output$CV[j,] - c(0,output$CV[j,-dim(output$CV)[2]]) # vector of indicators whether variable X_j included in V^(t) along iterations t
    selected[as.character(j),] <-  output$CS[j,] - c(0,output$CS[j,-dim(output$CS)[2]]) # vector of indicators whether variable X_j selected in S^(t) along iterations t
  }
  
  PF <- c()
  for (j in S_star) {
    if (all(selected[as.character(j),][sampled[as.character(j),]==1]==1)) PF = c(PF,j) # check whether PF holds
  }
  
  
  OIP_cur <- PF
  OIP_prime_cur <- PF
  N_cur <- c()
  counter_OIP <- 0
  counter_OIP_prime <- 0
  
  repeat {
    OIP_old <- OIP_cur
    OIP_prime_old <- OIP_prime_cur
    
    all_included_sampled <- colSums(sampled[as.character(OIP_cur),,drop=F]) == length(OIP_cur)
    
    for (j in setdiff(1:p, union(N_cur,S_star))) {
      if (all(selected[as.character(j),][sampled[as.character(j),]==1 &  all_included_sampled]==0)) N_cur = c(N_cur,j) 
    }
    
    all_excluded_noise <- colSums(sampled[as.character(N_cur),,drop=F]) == 0
    
    for (j in setdiff(S_star,OIP_cur)) {
      if (all(selected[as.character(j),][sampled[as.character(j),]==1 
                                         &  all_included_sampled 
                                         & all_excluded_noise]==1)) OIP_prime_cur = c(OIP_prime_cur,j) 
      if (all(selected[as.character(j),][sampled[as.character(j),]==1 
                                         &  all_included_sampled]==1)) OIP_cur = c(OIP_cur,j) 
    }
    if (setequal(OIP_old,OIP_cur) & setequal(OIP_prime_old, OIP_prime_cur)) { 
      OIP_prime = unique(sort(OIP_prime_cur)) 
      OIP = sort(OIP_cur) 
      break  
    }
    if (!setequal(OIP_old,OIP_cur)) counter_OIP <- counter_OIP + 1
    if (!setequal(OIP_prime_old,OIP_prime_cur)) counter_OIP_prime <- counter_OIP_prime + 1
  }
  
  PF_holds = setequal(PF, S_star)
  OIP_holds = setequal(OIP, S_star)
  OIP_prime_holds = setequal(OIP_prime, S_star)
  Thresholded_equals_best = setequal(which(output$relfreq.final > rho), S_star)
  Thresholded_subseteq_best = all(is.element(which(output$relfreq.final > rho), S_star))
  
           
  return(list(PF = PF, OIP = OIP, OIP_prime = OIP_prime, S_star = S_star, 
              counter_OIP = counter_OIP, counter_OIP_prime = counter_OIP_prime,
              PF_holds = PF_holds, OIP_holds = OIP_holds, OIP_prime_holds = OIP_prime_holds, 
              Thresholded_equals_best = Thresholded_equals_best, Thresholded_subseteq_best = Thresholded_subseteq_best))
}




AdaSubWrapper <-function (data,Iter_vec,K=100,q=10,criterion="EBIC",const=0,savings=10, U_C = 40) { 
  CV.null = numeric(p)
  CS.null = numeric(p)
  values = rep(NA,sum(Iter_vec))
  best.models = vector("list", sum(Iter_vec))
  Iter.cur = 0
 for (Iter in Iter_vec){
   output.cur = AdaSub(data = data, Iter=Iter, K=K, q=q, criterion=criterion, const=const, savings=savings, U_C = 40, CV.null=CV.null, CS.null=CS.null) 
   CV.null = CV.null + K * output.cur$CV.cur
   CS.null = CS.null + K * output.cur$CS.cur
   values[(Iter.cur+1):(Iter.cur+Iter)] = output.cur$values
   for (i in 1:Iter) best.models[[Iter.cur+i]] = output.cur$best.models[[i]]
   par(mfrow=c(1,1))
   plot((Iter.cur+1):(Iter.cur+Iter),values[(Iter.cur+1):(Iter.cur+Iter)],pch=20,
        main=paste("Dynamic plot of EBIC values for gamma =", const, "(from iteration", Iter.cur + 1 ,"to",  Iter.cur + Iter ,")\n Red line: Thresholded model (0.9)"),xlab="Iteration",ylab="EBIC value of selected model")
   abline(h=EBIC(data,which((q+K*CS.null)/(p+K*CV.null)>0.9),const),col="red")
   Iter.cur = Iter.cur + Iter
   # print(paste("Current number of iterations: ", Iter_cur))
 }
  best.S = best.models[[which.min(values)]]
 return(list(output=output.cur, values=values, best.models= best.models, best.S =best.S))
}

## Simulate data from linear model with different correlation settings (see paper for description) ###############

###### Simulate data with global correlation structure

simdata.global.corr <- function (n,p,beta,sigma.normal,corr=0) {
  
  mu=rep(0,p)
  Sigma=diag(1-corr,p)+matrix(corr,p,p)
  if (p<=140) x = mvrnorm(n , mu, Sigma) else x = rmvn(n,mu, Sigma)
  linpred=x%*%beta
  y=rnorm(n,linpred,sigma.normal)
  
  return(list(x=x,y=y))
}

###### Simulate data with block correlation structure (same correlation "corr" in blocks, no correlation between blocks) )

simdata.block.corr <- function (n,p,blocks,beta,sigma.normal,corr=0) { # "blocks" = no. of blocks
  
  mu=rep(0,p)
  Sigma=matrix(0,p,p)
  for (k in 1:p){
    for (m in 1:p){
      if ((k-m) %% blocks == 0 ) Sigma[k,m]= corr 
    }
  }
  
  Sigma=diag(1-corr,p)+Sigma
  if (p<=140) x = mvrnorm(n , mu, Sigma) else x = rmvn(n, mu, Sigma)
  
  linpred=x%*%beta
  y=rnorm(n,linpred,sigma.normal)
  
  return(list(x=x,y=y))
}

###### Simulate data with Toeplitz correlation structure

simdata.toeplitz.corr <- function (n,p,beta,sigma.normal,corr=0) {
  
  mu=rep(0,p)
  help=numeric(p)
  for (k in 1:p) help[k]=corr^(k-1)
  
  Sigma=toeplitz(help)
  #x = mvrnorm(n , mu, Sigma)
  if (p<=140) x = mvrnorm(n , mu, Sigma) else x = rmvn(n, mu, Sigma)
  
  linpred=x%*%beta
  y=rnorm(n,linpred,sigma.normal)
  
  return(list(x=x,y=y))
}



########################################################################
############# Functions for simulation study ###########################
########################################################################


method_init <- function (Scenarios,n.start,stepsize){ # Initialize performance measures for a certain method
  
  matr = matrix(0,nrow=7,ncol=Scenarios)
  rownames(matr) = c("false_pos_mean",
                     "false_neg_mean",
                     "estimation_error_mean", 
                     "prediction_error_mean",
                     "times",
                     "equal_count_0",
                     "equal_count_best"
  )
  matr[c("equal_count_0", "equal_count_best"),] = 0
  colnames(matr) =  n.start  + stepsize * ((1:Scenarios)-1)
  return(matr)
}

method_init_scenario <- function (N_Ex){ # Initialize performance measures for a certain method for one Scenario
  
  matr = matrix(0,nrow=5,ncol=N_Ex)
  rownames(matr) = c("false_pos",
                     "false_neg",
                     "estimation_error", 
                     "prediction_error",
                     "times")
  return(matr)
}

# beta=beta : Provide fitted beta (see e.g. Adaptive Lasso), else use (approximate) ML estimator 
method_results_scenario <- function (model,beta1,s0,data.test,n_test,time,beta=NULL,data){ # Compute performance measures for a certain method for one Scenario
  
  false_pos = sum(beta1[model]==0)
  false_neg = s0-sum(beta1[model]!=0)
  
  if (is.null(beta)){
    beta.hat =  beta.hat(data,model)
    #estimation_error = sqnorm2 ( c(0,beta1) - beta.hat)
    #test_y.hat = beta.hat_glm(data.test,model,family=family)$y.hat 
  } else {
    beta.hat =  beta
  }
  estimation_error = sqnorm2 ( c(0,beta1) - beta.hat)

  x.cur=cbind(c(rep(1,n_test)),data.test$x)  
  linpred = x.cur %*% beta.hat
  test_y.hat = as.vector(linpred)
  
  ### Changed to RMSE ###########
  prediction_error   = sqrt( 1/n_test * sqnorm2 ( data.test$y - test_y.hat))
  
  vect = c(false_pos,false_neg,estimation_error,prediction_error,time)
  return(vect)
}




########################################################################
############# Additional functions #####################################
########################################################################


beta.hat <- function(data,indices) {
  #x.matrix = cbind(c(rep(1,n)),data$x)  
  n=nrow(data$x)
  x.cur=data$x[,indices]
  x.cur=cbind(c(rep(1,n)),x.cur)                 #include intercept!
  hat.beta.cur = solve(t(x.cur)%*%x.cur) %*% t(x.cur) %*% data$y 
  hat.beta = numeric(ncol(data$x)+1)
  hat.beta[c(1,indices+1)] = hat.beta.cur
  #hat.y = x.matrix %*% hat.beta
  return(hat.beta)
}

y.hat <- function(data,beta.hat) {
  n=nrow(data$x)
  x.matrix = cbind(c(rep(1,n)),data$x) 
  return (x.matrix %*% beta.hat)
}

sqnorm2 <- function(x) return(sum(x^2)) # squared 2-norm

# Compute extended BIC for given estimate beta of length p+1 
EBIC_beta <- function(data, beta,const) {  
  n=nrow(data$x)
  p=ncol(data$x)
  df=sum(beta!=0) 
  x.cur=data$x
  x.cur=cbind(c(rep(1,n)),x.cur)                 #include intercept!
  RSS = sum( (data$y - x.cur %*% beta)^2 )
  BIC=n*(1+log(2*pi)+log(RSS/n))+log(n)*df #compute BIC, nb. of parameters=df+1 (+sigma^2)
  EBIC=BIC+2*df*const*log(p)
  return(EBIC)
}






 
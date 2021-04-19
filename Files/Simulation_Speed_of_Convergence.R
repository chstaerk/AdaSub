
AdaSub_evolution <-function (Iter, p = 1000, s = 2, K = 100, q = 10, rho = 0.9, assumption = "PF", best_only=FALSE) { 
  
  CV.cur = numeric(p)  
  CS.cur = numeric(p) 
  
  CV.null = numeric(p)
  CS.null = numeric(p)
  
  #relfreq=rep(q/p,p)
  relfreq = (q+CS.null+CS.cur)/(p+CV.null+CV.cur)
  
  V.size<-numeric(Iter)
  S.size<-numeric(Iter)
  
  S = integer(0)
  
  iteration_best_model_found <- NA
  iteration_thresholded_model <- NA
  
  for (t in 1:Iter) {
    
    b = rbinom(p,1,relfreq)
    V = which(b==1)
    
    if (length(V)>1){
      if (assumption=="PF") S = intersect(V, 1:s)
      if (assumption=="OIP") {
        included = is.element(1:s,sort(V)) 
        S=c()
        for (j in sort(intersect(V, 1:s))) {
          if (all(included[1:j])) S = c(S,j) 
        }
      }
      CV.cur[V] = CV.cur[V] +1
      CS.cur[S] = CS.cur[S] +1
      V.size[t]=length(V)
      S.size[t]=length(S)
      relfreq[V] = (q + CS.null[V] + K*CS.cur[V])/(p +CV.null[V] + K*CV.cur[V])
      
      if (is.na(iteration_thresholded_model)) {
        if (setequal(which(relfreq>rho),1:s))
          iteration_thresholded_model <- t
      }
      
      if (is.na(iteration_best_model_found)) {
        if (setequal(S,1:s)) 
          iteration_best_model_found <- t
      }
    }
    
    if (!is.na(iteration_thresholded_model) & !is.na(iteration_best_model_found)) break()
    if (best_only) {
      if (!is.na(iteration_best_model_found)) break()
    }
    
  }
  relfreq.final=relfreq
  return(list(relfreq.final=relfreq.final,S.size=S.size,V.size=V.size, CV.cur=CV.cur, CS.cur=CS.cur,
              iteration_best_model_found = iteration_best_model_found, 
              iteration_thresholded_model = iteration_thresholded_model))
}    


#### Basic setting 

q = 10
K = 200
p = 2000 
rho = 0.9
s = 3 
Iter = 50000



choices <-c( "500 simulations per setting (as in paper)", "Custom number of simulations per setting")
mychoice <- menu( choices , graphics=TRUE, title="How many simulations per setting should be considered?" )

if (mychoice==1) {
  nsim=500 
}
if (mychoice==2) {
  nsim = -1
  while(nsim < 1 ){
    nsim <- readline("Enter the number of simulations per setting: ")
    nsim <- ifelse(grepl("\\D",nsim),-1,as.integer(nsim))
    if(is.na(nsim)){break}  # breaks when hit enter
  }
}

#nsim <- 500 

#### Different settings

p_values <- c(500,1000,2000,4000,8000,16000) 
s_values <- c(1,2,3,4,5,7,10,20,40)
q_values <- c(2.5,5,10,20,40)
K_values <- c(10,50,100,200,1000,2000,5000,10000)


p_settings <- length(p_values)
s_settings <- length(s_values)
q_settings <- length(q_values)
K_settings <- length(K_values)


iteration_best_model_found_mat_p_PF <- matrix(NA, nrow = nsim, ncol = p_settings)
iteration_thresholded_model_mat_p_PF <- matrix(NA, nrow = nsim, ncol = p_settings)
iteration_best_model_found_mat_p_OIP <- matrix(NA, nrow = nsim, ncol = p_settings)
iteration_thresholded_model_mat_p_OIP <- matrix(NA, nrow = nsim, ncol = p_settings)
colnames(iteration_best_model_found_mat_p_PF) <-  p_values
colnames(iteration_thresholded_model_mat_p_PF) <-  p_values
colnames(iteration_best_model_found_mat_p_OIP) <-  p_values
colnames(iteration_thresholded_model_mat_p_OIP) <-  p_values


iteration_best_model_found_mat_s_PF <- matrix(NA, nrow = nsim, ncol = s_settings)
iteration_thresholded_model_mat_s_PF <- matrix(NA, nrow = nsim, ncol = s_settings)
iteration_best_model_found_mat_s_OIP <- matrix(NA, nrow = nsim, ncol = s_settings)
iteration_thresholded_model_mat_s_OIP <- matrix(NA, nrow = nsim, ncol = s_settings)
colnames(iteration_best_model_found_mat_s_PF) <-  s_values
colnames(iteration_thresholded_model_mat_s_PF) <-  s_values
colnames(iteration_best_model_found_mat_s_OIP) <-  s_values
colnames(iteration_thresholded_model_mat_s_OIP) <-  s_values

iteration_best_model_found_mat_q_PF <- matrix(NA, nrow = nsim, ncol = q_settings)
iteration_thresholded_model_mat_q_PF <- matrix(NA, nrow = nsim, ncol = q_settings)
iteration_best_model_found_mat_q_OIP <- matrix(NA, nrow = nsim, ncol = q_settings)
iteration_thresholded_model_mat_q_OIP <- matrix(NA, nrow = nsim, ncol = q_settings)
colnames(iteration_best_model_found_mat_q_PF) <-  q_values
colnames(iteration_thresholded_model_mat_q_PF) <-  q_values
colnames(iteration_best_model_found_mat_q_OIP) <-  q_values
colnames(iteration_thresholded_model_mat_q_OIP) <-  q_values

iteration_best_model_found_mat_K_PF <- matrix(NA, nrow = nsim, ncol = K_settings)
iteration_thresholded_model_mat_K_PF <- matrix(NA, nrow = nsim, ncol = K_settings)
iteration_best_model_found_mat_K_OIP <- matrix(NA, nrow = nsim, ncol = K_settings)
iteration_thresholded_model_mat_K_OIP <- matrix(NA, nrow = nsim, ncol = K_settings)
colnames(iteration_best_model_found_mat_K_PF) <-  K_values
colnames(iteration_thresholded_model_mat_K_PF) <-  K_values
colnames(iteration_best_model_found_mat_K_OIP) <-  K_values
colnames(iteration_thresholded_model_mat_K_OIP) <-  K_values


#### Vary p 
start_time <- Sys.time()
set.seed(1234)
for (k in 1:p_settings) {
  for (i in 1:nsim) {
    output <- AdaSub_evolution(Iter = Iter, p = p_values[k], s = s, K = K, q = q, rho = rho, assumption = "PF")
    iteration_best_model_found_mat_p_PF[i,k] <- output$iteration_best_model_found
    iteration_thresholded_model_mat_p_PF[i,k] <- output$iteration_thresholded_model
    
    output <- AdaSub_evolution(Iter = Iter, p = p_values[k], s = s, K = K, q = q, rho = rho, assumption = "OIP")
    iteration_best_model_found_mat_p_OIP[i,k] <- output$iteration_best_model_found
    iteration_thresholded_model_mat_p_OIP[i,k] <- output$iteration_thresholded_model
  }
}
end_time <- Sys.time()
time_p <- end_time - start_time

#### Vary s
start_time <- Sys.time()
set.seed(1234)
for (k in 1:s_settings) {
  for (i in 1:nsim) {
    output <- AdaSub_evolution(Iter = Iter, p = p, s = s_values[k], K = K, q = q, rho = rho, assumption = "PF")
    iteration_best_model_found_mat_s_PF[i,k] <- output$iteration_best_model_found
    iteration_thresholded_model_mat_s_PF[i,k] <- output$iteration_thresholded_model
    if (k<=5) {
      output <- AdaSub_evolution(Iter = Iter, p = p, s = s_values[k], K = K, q = q, rho = rho, assumption = "OIP")
      iteration_best_model_found_mat_s_OIP[i,k] <- output$iteration_best_model_found
      iteration_thresholded_model_mat_s_OIP[i,k] <- output$iteration_thresholded_model
    }
    if (k==6) {
      output <- AdaSub_evolution(Iter = Iter, p = p, s = s_values[k], K = K, q = q, rho = rho, assumption = "OIP",  best_only=TRUE)
      iteration_best_model_found_mat_s_OIP[i,k] <- output$iteration_best_model_found
      iteration_thresholded_model_mat_s_OIP[i,k] <- output$iteration_thresholded_model
    }
  }
}
end_time <- Sys.time()
time_s <- end_time - start_time

#### Vary q
start_time <- Sys.time()
set.seed(1234)
for (k in 1:q_settings) {
  for (i in 1:nsim) {
    output <- AdaSub_evolution(Iter = Iter, p = p, s = s, K = K, q = q_values[k], rho = rho, assumption = "PF")
    iteration_best_model_found_mat_q_PF[i,k] <- output$iteration_best_model_found
    iteration_thresholded_model_mat_q_PF[i,k] <- output$iteration_thresholded_model
    
    output <- AdaSub_evolution(Iter = Iter, p = p, s = s, K = K, q = q_values[k], rho = rho, assumption = "OIP")
    iteration_best_model_found_mat_q_OIP[i,k] <- output$iteration_best_model_found
    iteration_thresholded_model_mat_q_OIP[i,k] <- output$iteration_thresholded_model
  }
}
end_time <- Sys.time()
time_q <- end_time - start_time

#### Vary K
start_time <- Sys.time()
set.seed(1234)
for (k in 1:K_settings) {
  for (i in 1:nsim) {
    output <- AdaSub_evolution(Iter = Iter, p = p, s = s, K = K_values[k], q = q, rho = rho, assumption = "PF")
    iteration_best_model_found_mat_K_PF[i,k] <- output$iteration_best_model_found
    iteration_thresholded_model_mat_K_PF[i,k] <- output$iteration_thresholded_model
    
    output <- AdaSub_evolution(Iter = Iter, p = p, s = s, K = K_values[k], q = q, rho = rho, assumption = "OIP")
    iteration_best_model_found_mat_K_OIP[i,k] <- output$iteration_best_model_found
    iteration_thresholded_model_mat_K_OIP[i,k] <- output$iteration_thresholded_model
  }
}
end_time <- Sys.time()
time_K <- end_time - start_time

time_p
time_q
time_s
time_K

##
Save_iter <- list(iteration_best_model_found_mat_p_PF = iteration_best_model_found_mat_p_PF,
                  iteration_thresholded_model_mat_p_PF = iteration_thresholded_model_mat_p_PF,
                  iteration_best_model_found_mat_p_OIP = iteration_best_model_found_mat_p_OIP,
                  iteration_thresholded_model_mat_p_OIP = iteration_thresholded_model_mat_p_OIP,
                  
                  iteration_best_model_found_mat_s_PF = iteration_best_model_found_mat_s_PF,
                  iteration_thresholded_model_mat_s_PF = iteration_thresholded_model_mat_s_PF,
                  iteration_best_model_found_mat_s_OIP = iteration_best_model_found_mat_s_OIP,
                  iteration_thresholded_model_mat_s_OIP = iteration_thresholded_model_mat_s_OIP,
                  
                  iteration_best_model_found_mat_q_PF = iteration_best_model_found_mat_q_PF,
                  iteration_thresholded_model_mat_q_PF = iteration_thresholded_model_mat_q_PF,
                  iteration_best_model_found_mat_q_OIP = iteration_best_model_found_mat_q_OIP,
                  iteration_thresholded_model_mat_q_OIP = iteration_thresholded_model_mat_q_OIP,
                  
                  iteration_best_model_found_mat_K_PF = iteration_best_model_found_mat_K_PF,
                  iteration_thresholded_model_mat_K_PF = iteration_thresholded_model_mat_K_PF,
                  iteration_best_model_found_mat_K_OIP = iteration_best_model_found_mat_K_OIP,
                  iteration_thresholded_model_mat_K_OIP = iteration_thresholded_model_mat_K_OIP,
                  
                  q = q,
                  K = K,
                  p = p,
                  rho = rho,
                  s = s,
                  Iter = Iter,
                  nsim = nsim,
                  
                  p_values = p_values,
                  s_values = s_values, 
                  q_values = q_values, 
                  K_values = K_values,
                  
                  time_p = time_p,
                  time_s = time_s,
                  time_q = time_q,
                  time_K = time_K
)

#setwd("C:/Users/Staerk/ScieboBonn/AdaSub paper for EJS/R simulation saves")
#save( Save_iter, file = "Save_iter_200.RData" )
#save( Save_iter, file = "Save_iter_500_s7.RData" )
#load("Save_iter_500_s7.RData")

#dim(Save_iter$iteration_best_model_found_mat_p_PF)

### Plots

#dev.off()
win.graph(width = 14, height = 8, pointsize = 8)

font = 1.5
par(las=1,lwd=2,cex.main=font, cex.lab=font, cex.axis=font, cex = 2)


############################################################################
m <- matrix(c(1,2,3,4,5,6,7,8,9,10,11,11,12,12),nrow = 7,ncol = 2,byrow = TRUE)
layout(mat = m, heights = c( 0.235, 0.235, 0.235, 0.235 ,0.025, 0.05, 0.025))
par(las=1, cex.axis = 1.1)
############################################################################

#par(mfrow=c(4,2), las=1)
ylim1 <- c(0,11000)
ylim2a <- c(0,3700)
ylim2b <- c(0,3700)
ylim3 <- c(0,3500)
ylim4 <- c(0,5000)

plot(Save_iter$p_values, colMeans(Save_iter$iteration_best_model_found_mat_p_PF), 
     ylim = ylim1, 
     type ="b", xlab ="Number of covariates p", ylab ="", main= "Best model: Mean iterations for varying p", pch=19, xaxt="n")
axis(1, at = Save_iter$p_values, las=1)
lines(Save_iter$p_values, colMeans(Save_iter$iteration_best_model_found_mat_p_OIP), col ="red", type ="b", lty = 2)
plot(Save_iter$p_values, colMeans(Save_iter$iteration_thresholded_model_mat_p_PF), 
     ylim =  ylim1
     , type ="b", xlab ="Number of covariates p", ylab ="", main= "Thresholded model: Mean iterations for varying p", pch=19, xaxt="n")
axis(1, at = Save_iter$p_values, las=1)
lines(Save_iter$p_values, colMeans(Save_iter$iteration_thresholded_model_mat_p_OIP), col ="red", type ="b", lty = 2)

plot(Save_iter$s_values, colMeans(Save_iter$iteration_best_model_found_mat_s_PF), 
     ylim =  ylim2a, 
     type ="b", xlab ="Number of variables in S*", ylab ="", main= "Best model: Mean iterations for varying s*", pch=19, xaxt="n")
axis(1, at = Save_iter$s_values, las=1)
lines(Save_iter$s_values, colMeans(Save_iter$iteration_best_model_found_mat_s_OIP), col ="red", type ="b", lty = 2)
plot(Save_iter$s_values, colMeans(Save_iter$iteration_thresholded_model_mat_s_PF), 
     ylim =  ylim2b
     , type ="b", xlab ="Number of variables in S*", ylab ="", main= "Thresholded model: Mean iterations for varying s*", pch=19, xaxt="n")
lines(Save_iter$s_values[1:6], colMeans(Save_iter$iteration_thresholded_model_mat_s_OIP)[1:6], col ="red", type ="b", lty = 2)
axis(1, at = Save_iter$s_values, las=1)

plot(Save_iter$q_values, colMeans(Save_iter$iteration_best_model_found_mat_q_PF), 
     ylim =  ylim3, 
     type ="b", xlab ="Initial expected search size q", ylab ="", main= "Best model: Mean iterations for varying q", pch=19, xaxt="n")
axis(1, at = Save_iter$q_values, las=1)
lines(Save_iter$q_values, colMeans(Save_iter$iteration_best_model_found_mat_q_OIP), col ="red", type ="b", lty = 2)
plot(Save_iter$q_values, colMeans(Save_iter$iteration_thresholded_model_mat_q_PF), 
     ylim =  ylim3
     , type ="b", xlab ="Initial expected search size q", ylab ="", main= "Thresholded model: Mean iterations for varying q", pch=19, xaxt="n")
axis(1, at = Save_iter$q_values, las=1)
lines(Save_iter$q_values, colMeans(Save_iter$iteration_thresholded_model_mat_q_OIP), col ="red", type ="b", lty = 2)


plot(Save_iter$K_values, colMeans(Save_iter$iteration_best_model_found_mat_K_PF, na.rm = TRUE), 
     ylim =  ylim4, 
     type ="b", xlab ="Learning rate K", ylab ="", main= "Best model: Mean iterations for varying K", pch=19, xaxt="n")
axis(1, at = Save_iter$K_values, las=1)
lines(Save_iter$K_values, colMeans(Save_iter$iteration_best_model_found_mat_K_OIP), col ="red", type ="b", lty = 2)
plot(Save_iter$K_values, colMeans(Save_iter$iteration_thresholded_model_mat_K_PF), 
     ylim =  ylim4
     , type ="b", xlab ="Learning rate K", ylab ="", main= "Thresholded model: Mean iterations for varying K", pch=19, xaxt="n")
axis(1, at = Save_iter$K_values, las=1)
lines(Save_iter$K_values, colMeans(Save_iter$iteration_thresholded_model_mat_K_OIP), col ="red", type ="b", lty = 2)

par(mar=c(0,0,0,0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
text(1, 1, "Left: Mean iterations needed so that best AdaSub model equals S*", cex=1.4)

plot(1, type = "n", axes=FALSE, xlab="", ylab="")
text(1, 1, "Right: Mean iterations needed so that thresholded model contains S*", cex=1.4)

# Plot legend
par(mar=c(0,0,0,0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend(x = "center",inset = 0,
       legend = c("finite-sample PF          ", "minimal OIP"), 
       col=c("black", "red"), cex=1.7, box.lty=0,ncol = 2,
       pch=c(19,1), pt.bg = 'white', lty=c(1, 2))

plot(1, type = "n", axes=FALSE, xlab="", ylab="")
text(1, 1, "In each plot a single parameter is varied with the remaining ones constant at p = 2000, s* = 3, q = 10, K = 200.", cex=1.4)

par(mar=c(5.1, 4.1, 4.1, 2.1))
par(las=1, cex.axis = 1)




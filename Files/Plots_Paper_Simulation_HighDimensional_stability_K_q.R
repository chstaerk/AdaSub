


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
  
#setwd("C://Users//Staerk//sciebo//AdaSub paper for S&C//R simulation saves")

load("./Files/R_simulation_saves/Stability_EBIC1_Global00_p_increasing_600_2000_simulations20.RData")
Save_EBIC1_00 = Save1
ident = Count_ident(Save_EBIC1_00)
AdaSub_ident_EBIC1_00  = ident$AdaSub_ident
AdaSubBest_ident_EBIC1_00  = ident$AdaSubBest_ident
AdaSub_var_EBIC1_00  = ident$AdaSub_var
AdaSubBest_var_EBIC1_00  = ident$AdaSubBest_var

load("./Files/R_simulation_saves/Stability_EBIC06_Global00_p_increasing_600_2000_simulations20.RData")
Save_EBIC06_00 = Save1
ident = Count_ident(Save_EBIC06_00)
AdaSub_ident_EBIC06_00  = ident$AdaSub_ident
AdaSubBest_ident_EBIC06_00  = ident$AdaSubBest_ident
AdaSub_var_EBIC06_00  = ident$AdaSub_var
AdaSubBest_var_EBIC06_00  = ident$AdaSubBest_var


load("./Files/R_simulation_saves/Stability_EBIC1_Toeplitz09_p_increasing_600_2000_simulations20.RData")
Save_EBIC1_09 = Save1
ident = Count_ident(Save_EBIC1_09)
AdaSub_ident_EBIC1_09  = ident$AdaSub_ident
AdaSubBest_ident_EBIC1_09  = ident$AdaSubBest_ident
AdaSub_var_EBIC1_09  = ident$AdaSub_var
AdaSubBest_var_EBIC1_09  = ident$AdaSubBest_var

load("./Files/R_simulation_saves/Stability_EBIC06_Toeplitz09_p_increasing_600_2000_simulations20.RData")
Save_EBIC06_09 = Save1
ident = Count_ident(Save_EBIC06_09)
AdaSub_ident_EBIC06_09  = ident$AdaSub_ident
AdaSubBest_ident_EBIC06_09  = ident$AdaSubBest_ident
AdaSub_var_EBIC06_09  = ident$AdaSub_var
AdaSubBest_var_EBIC06_09  = ident$AdaSubBest_var

pchs = c(4,8,1,2)
cex.size = 1.4
ltys = c(1,4,3,2)
colors = c("black","gray40","orange","darkorange2")

font = 1.5
par(las=1,lwd=2,cex.main=font, cex.lab=font, cex.axis=font)

m <- matrix(c(1,2,3,4,5,5),nrow = 3,ncol = 2,byrow = TRUE)

layout(mat = m, heights = c(0.45,0.45,0.1))


xaxis = Save1$n.start  + Save1$stepsize * ((1:Save1$Scenarios)-1)
plot(xaxis,rowMeans(AdaSub_ident_EBIC06_00)/Save1$N_Rep,ylim = c(0,1),pch=pchs[1],cex=cex.size, xlab="n",ylab="",main="Indep. structure: Mean rel. freq. of agreement",xaxt="n",col=colors[1], type="b", lty = ltys[1])
points(xaxis,rowMeans(AdaSub_ident_EBIC1_00)/Save1$N_Rep,pch=pchs[2],cex=cex.size,col=colors[2], type="b", lty = ltys[2])
points(xaxis,rowMeans(AdaSubBest_ident_EBIC06_00)/Save1$N_Rep,pch=pchs[3],cex=cex.size,col=colors[3], type="b", lty = ltys[3])
points(xaxis,rowMeans(AdaSubBest_ident_EBIC1_00)/Save1$N_Rep,pch=pchs[4],cex=cex.size,col=colors[4], type="b", lty = ltys[4])
axis(1,at=xaxis)


plot(xaxis,rowMeans(AdaSub_ident_EBIC06_09)/Save1$N_Rep,ylim = c(0,1),pch=pchs[1],cex=cex.size, xlab="n",ylab="",main="Toeplitz structure: Mean rel. freq. of agreement",xaxt="n",col=colors[1], type="b", lty = ltys[1])
points(xaxis,rowMeans(AdaSub_ident_EBIC1_09)/Save1$N_Rep,pch=pchs[2],cex=cex.size,col=colors[2], type="b", lty = ltys[2])
points(xaxis,rowMeans(AdaSubBest_ident_EBIC06_09)/Save1$N_Rep,pch=pchs[3],cex=cex.size,col=colors[3], type="b", lty = ltys[3])
points(xaxis,rowMeans(AdaSubBest_ident_EBIC1_09)/Save1$N_Rep,pch=pchs[4],cex=cex.size,col=colors[4], type="b", lty = ltys[4])
axis(1,at=xaxis)



xaxis = Save1$n.start  + Save1$stepsize * ((1:Save1$Scenarios)-1)
plot(xaxis,rowMeans(AdaSub_var_EBIC06_00),ylim = c(0,1),pch=pchs[1],cex=cex.size, xlab="n",ylab="",main="Indep. structure: Mean variance of model size",xaxt="n",col=colors[1], type="b", lty = ltys[1])
points(xaxis,rowMeans(AdaSub_var_EBIC1_00),pch=pchs[2],cex=cex.size,col=colors[2], type="b", lty = ltys[2])
points(xaxis,rowMeans(AdaSubBest_var_EBIC06_00),pch=pchs[3],cex=cex.size,col=colors[3], type="b", lty = ltys[3])
points(xaxis,rowMeans(AdaSubBest_var_EBIC1_00),pch=pchs[4],cex=cex.size,col=colors[4], type="b", lty = ltys[4])
axis(1,at=xaxis)


plot(xaxis,rowMeans(AdaSub_var_EBIC06_09),ylim = c(0,1),pch=pchs[1],cex=cex.size, xlab="n",ylab="",main="Toeplitz structure: Mean variance of model size",xaxt="n",col=colors[1], type="b", lty = ltys[1])
points(xaxis,rowMeans(AdaSub_var_EBIC1_09),pch=pchs[2],cex=cex.size,col=colors[2], type="b", lty = ltys[2])
points(xaxis,rowMeans(AdaSubBest_var_EBIC06_09),pch=pchs[3],cex=cex.size,col=colors[3], type="b", lty = ltys[3])
points(xaxis,rowMeans(AdaSubBest_var_EBIC1_09),pch=pchs[4],cex=cex.size,col=colors[4], type="b", lty = ltys[4])
axis(1,at=xaxis)





legendnames = c(expression(paste("AdaSubThres ", gamma," = 0.6  ")), 
                expression(paste("AdaSubThres ", gamma," = 1  ")), 
                expression(paste("AdaSubBest ", gamma," = 0.6  ")), 
                expression(paste("AdaSubBest ", gamma," = 1 ")))

par(mar=c(0,0,0,0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
plot_colors <- colors
legend(x = "center",inset = 0,
       legend = legendnames, 
       col=plot_colors, cex=1.5, box.lty=1, ncol = 4,
       pch=pchs, pt.bg = 'white', lty=ltys) #,
par(mar=c(5.1, 4.1, 4.1, 2.1))







































### Plot results for high-dimensional simulation study for CV-moodels (Lasso, AdaLasso, SCAD)############################
#setwd("C://Users//Staerk//scieboBonn//AdaSub paper for EJS//R simulation saves")

#load("./Files/R_simulation_saves/CV_Global00_p_increasing_400_2000_simulations500.RData")
#N_Ex = Save1$N_Ex

ploteq =FALSE

PC = FALSE

corr.type = "block"
corr = 0.5

#####
if (corr.type =="global" & corr==0) {
  load("./Files/R_simulation_saves/EBIC06_Global00_p_increasing_400_2000_simulations500.RData")
  Save06 = Save1
  load("./Files/R_simulation_saves/EBIC1_Global00_p_increasing_400_2000_simulations500.RData")
  load("./Files/R_simulation_saves/CV_Global00_p_increasing_400_2000_simulations500.RData")
}
if (corr.type =="global" & corr==0.7) {
  load("./Files/R_simulation_saves/EBIC06_Global07_p_increasing_400_2000_simulations500.RData")
  Save06 = Save1
  load("./Files/R_simulation_saves/EBIC1_Global07_p_increasing_400_2000_simulations500.RData")
  load("./Files/R_simulation_saves/CV_Global07_p_increasing_400_2000_simulations500.RData")
}
if (corr.type =="toeplitz" & corr==0.9) {
  load("./Files/R_simulation_saves/EBIC06_Toeplitz09_p_increasing_400_2000_simulations500.RData")
  Save06 = Save1
  load("./Files/R_simulation_saves/EBIC1_Toeplitz09_p_increasing_400_2000_simulations500.RData")
  load("./Files/R_simulation_saves/CV_Toeplitz09_p_increasing_400_2000_simulations500.RData")
}
if (corr.type =="block" & corr==0.5) {
  load("./Files/R_simulation_saves/EBIC06_Blocks10_05_p_increasing_400_2000_simulations500.RData")
  Save06 = Save1
  load("./Files/R_simulation_saves/EBIC1_Blocks10_05_p_increasing_400_2000_simulations500.RData")
  load("./Files/R_simulation_saves/CV_Blocks10_05_p_increasing_400_2000_simulations500.RData")
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
color_AdaLasso = "deepskyblue"
color_SCAD = "purple"

legendnames = c(expression(paste("AdaSubThres ", gamma," = 0.6  ")), 
                expression(paste("AdaSubBest ", gamma," = 0.6  ")),
                expression(paste("AdaSubThres ", gamma," = 1  ")),
                expression(paste("AdaSubBest ", gamma," = 1  ")),
                "Lasso CV", "AdaLasso CV", "SCAD CV")

pchs = c(4,1,8,2,3,15,6)
cex.size = 1.5
ltys = c(1,3,2,4,5,6,7)

#dev.off()
#win.graph(width = 14, height = 8, pointsize = 8)

font = 2
par(las=1,lwd=2,cex.main=font, cex.lab=font, cex.axis=font)
par(mfcol=c(2,3))

xaxis = Save_CV$n.start  + Save_CV$stepsize * ((1:Save_CV$Scenarios)-1)

############################################################################
m <- matrix(c(1,2,3,4,5,3,6,7,3),nrow = 3,ncol = 3,byrow = FALSE)
layout(mat = m, heights = c(0.45,0.45,0.1))
############################################################################



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









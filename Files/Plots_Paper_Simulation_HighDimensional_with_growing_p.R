
### Plot results for high-dimensional simulation study (see Section 5.2 of paper) ############################

#setwd("C://Users//Chris//scieboBonn//AdaSub paper for EJS//R simulation saves")
load("./Files/R_simulation_saves/EBIC1_Global00_p_increasing_400_2000_simulations500.RData")

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
color_AdaLasso = "deepskyblue"
color_SCAD = "purple"
color_Tilting = "pink"
color_PC = "yellow"

legendnames = c("AdaSubThres", "AdaSubBest", "StabSel", "Forward", "Lasso", "AdaLasso", "SCAD", "Tilting")

pchs = c(4,8,1,2,3,15,6,0)
cex.size = 1.5
ltys = c(1,2,3,4,5,6,7,8)

#dev.off()
#win.graph(width = 14, height = 8, pointsize = 8)

font = 2
par(las=1,lwd=2,cex.main=font, cex.lab=font, cex.axis=font)
par(mfcol=c(2,3))

xaxis = Save1$n.start  + Save1$stepsize * ((1:Save1$Scenarios)-1)
xaxis.Tilting = xaxis[xaxis<=120]

############################################################################
m <- matrix(c(1,2,3,4,5,3,6,7,3),nrow = 3,ncol = 3,byrow = FALSE)
layout(mat = m, heights = c(0.45,0.45,0.1))
############################################################################

#m <- matrix(c(1,2,3,4,1,5,6,4,1,7,8,4),nrow = 4,ncol = 3,byrow = FALSE)
#layout(mat = m, heights = c(0.05,0.425,0.425,0.1))


#par(mar=c(0, 0, 2, 0))
#plot(1, type = "n", axes=FALSE, xlab="", ylab="", cex.main=1.5,
#     main=paste("High-dimensional simulations (p=10n) for EBIC with gamma =",Save1$const,"with",Save1$N_Ex, "simulated datasets for", Save1$Scenarios, "scenarios"))
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

if (Save1$corr==0 & Save1$const==1) values_false_pos = c(0,0.5)
#######################################################
#values_false_pos = c(0, 0.5) # for EBIC1 for independence correlation structure ("Global00")
#######################################################

plot(xaxis,Save1$AdaSub_summary["false_pos_mean",],pch=pchs[1],cex=cex.size,ylim=values_false_pos, xlab="n",ylab="",main="Mean false positives",xaxt="n",col=color_AdaSub, type="b", lty = ltys[1])
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




plot(xaxis,Save1$AdaSub_summary["equal_count_0",]/Save1$N_Ex,pch=pchs[1],cex=cex.size,ylim=c(0,0.5), xlab="n",ylab="",main="Rel. freq. true model selected",xaxt="n",col=color_AdaSub, type="b", lty = ltys[1])
axis(1,at=xaxis)
points(xaxis,Save1$AdaSubBest_summary["equal_count_0",]/Save1$N_Ex,pch=pchs[2],cex=cex.size,col=color_AdaSubBest, type="b", lty = ltys[2])
points(xaxis,Save1$Stability_summary["equal_count_0",]/Save1$N_Ex,pch=pchs[3],cex=cex.size,col=color_Stability, type="b", lty = ltys[3])
points(xaxis,Save1$Forward_summary["equal_count_0",]/Save1$N_Ex,pch=pchs[4],cex=cex.size,col=color_Forward, type="b", lty = ltys[4])
points(xaxis,Save1$Lasso_summary["equal_count_0",]/Save1$N_Ex,pch=pchs[5],cex=cex.size,col=color_Lasso, type="b", lty = ltys[5])
points(xaxis,Save1$AdaLasso_summary["equal_count_0",]/Save1$N_Ex,pch=pchs[6],cex=cex.size,col=color_AdaLasso, type="b", lty = ltys[6])
points(xaxis,Save1$SCAD_summary["equal_count_0",]/Save1$N_Ex,pch=pchs[7],cex=cex.size,col=color_SCAD, type="b", lty = ltys[7])
points(xaxis.Tilting,Save1$Tilting_summary["equal_count_0",]/Save1$N_Ex,pch=pchs[8],cex=cex.size,col=color_Tilting, type="b", lty = ltys[8])
if (PC){
points(xaxis,Save1$PC_summary["equal_count_0",]/Save1$N_Ex,pch=20,cex=2.5,col=color_PC)
}
points(xaxis,Save1$AdaSubBest_summary["equal_count_0",]/Save1$N_Ex,pch=pchs[2],cex=cex.size,col=color_AdaSubBest, type="b", lty = ltys[2])
points(xaxis,Save1$AdaSub_summary["equal_count_0",]/Save1$N_Ex,pch=pchs[1],cex=cex.size,col=color_AdaSub, type="b", lty = ltys[1])



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


























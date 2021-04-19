
### Plot results for high-dimensional simulation study (see Section 5.2 of paper) ############################

#setwd("C://Users//Staerk//scieboBonn//AdaSub paper for EJS//R simulation saves")

load("./Files/R_simulation_saves/All_info_BIC_Global00_p30_simulations100.RData")

Save_BIC_Global00 <- Save1

load("./Files/R_simulation_saves/All_info_BIC_Toeplitz09_p30_simulations100.RData")

Save_BIC_Toeplitz09 <- Save1

load("./Files/R_simulation_saves/All_info_EBIC06_Global00_p_increasing_400_2000_simulations500.RData")

Save_EBIC06_Global00 <- Save1

load("./Files/R_simulation_saves/All_info_EBIC06_Toeplitz09_p_increasing_400_2000_simulations500.RData")

Save_EBIC06_Toeplitz09 <- Save1

load("./Files/R_simulation_saves/All_info_EBIC1_Global00_p_increasing_400_2000_simulations500.RData")

Save_EBIC1_Global00 <- Save1

load("./Files/R_simulation_saves/All_info_EBIC1_Toeplitz09_p_increasing_400_2000_simulations500.RData")

Save_EBIC1_Toeplitz09 <- Save1

##################################


color_BIC = "red"
color_EBIC06 = "orange"
color_EBIC1 = "blue"

legendnames = c("BIC, c = 0", 
                "BIC, c = 0.9", 
                expression(paste("EBIC ", gamma," = 0.6, c = 0    ")), 
                expression(paste("EBIC ", gamma," = 0.6, c = 0.9    ")), 
                expression(paste("EBIC ", gamma," = 1, c = 0    ")), 
                expression(paste("EBIC ", gamma," = 1, c = 0.9    ")))

pchs = c(4,8,1,2,3,15,6,0)
cex.size = 1.5
ltys = c(1,2,3,4,5,6,7,8)

#dev.off()
#win.graph(width = 14, height = 8, pointsize = 8)

font = 1.5
par(las=1,lwd=2,cex.main=font, cex.lab=font, cex.axis=font)
par(mfcol=c(2,3))

xaxis = Save1$n.start  + Save1$stepsize * ((1:Save1$Scenarios)-1)
xaxis.Tilting = xaxis[xaxis<=120]

############################################################################
m <- matrix(c(1,2),nrow = 2,ncol = 1,byrow = TRUE)
layout(mat = m, heights = c(0.8,0.2))
############################################################################




plot(xaxis, rowMeans(Save_BIC_Global00$PF_count_rel, na.rm = TRUE), ylim =c(0,1), type ="b" , xaxt="n", xlab="n", ylab="", main="Mean proportion of variables in the 'best' model for which PF empirically holds", col = color_BIC, pch = pchs[1], lty = ltys[1])
points(xaxis,rowMeans(Save_BIC_Toeplitz09$PF_count_rel , na.rm = TRUE), type ="b", col = color_BIC, pch = pchs[2], lty = ltys[2])
points(xaxis,rowMeans(Save_EBIC06_Global00$PF_count_rel, na.rm = TRUE), type ="b", col = color_EBIC06, pch = pchs[3], lty = ltys[1])
points(xaxis,rowMeans(Save_EBIC06_Toeplitz09$PF_count_rel, na.rm = TRUE), type ="b", col = color_EBIC06, pch = pchs[4], lty = ltys[2])
points(xaxis,rowMeans(Save_EBIC1_Global00$PF_count_rel, na.rm = TRUE), type ="b", col = color_EBIC1, pch = pchs[5], lty = ltys[1])
points(xaxis,rowMeans(Save_EBIC1_Toeplitz09$PF_count_rel, na.rm = TRUE), type ="b", col = color_EBIC1, pch = pchs[6], lty = ltys[2])

axis(1,at=xaxis)




legend_order <- matrix(c(1:6),ncol=3, byrow = FALSE)
par(mar=c(0,0,0,0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
plot_colors <- c(color_BIC, color_BIC, color_EBIC06, color_EBIC06, color_EBIC1, color_EBIC1)
ltys <- rep(ltys[c(1,2)],3)
legend(x = "center",inset = 0,
       legend = legendnames[legend_order], 
       col=plot_colors[legend_order], cex=1.3, box.lty=1, ncol = 3,
       pch=pchs[legend_order], pt.bg = 'white', lty=ltys[legend_order]) #,
par(mar=c(5.1, 4.1, 4.1, 2.1))







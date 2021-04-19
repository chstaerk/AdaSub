mychoice <- menu( c("No","Yes"), graphics=TRUE, title="Install required libraries?" )

if (mychoice==2){ 
	install.packages("leaps")
	#install.packages("c060")
	install.packages("ncvreg")
	install.packages("glmnet")
	install.packages("MASS")
	install.packages("tilting")
	#source("http://bioconductor.org/biocLite.R")
	#biocLite("RBGL")
	#install.packages("pcalg")
	install.packages("stabs")
	install.packages("mvnfast")
} 

library("leaps")
#library("c060")
library("ncvreg")
library("glmnet")
library("MASS")
library("tilting")
#library("pcalg")
library("mvnfast")
library("stabs")

#savewd = getwd()
#setwd(paste0(getwd(), "/Files"))

source("./Files/AdaSub_main_functions.R")

choices <-c( "Low Dimensional Simulation ", "High Dimensional with growing p Simulation ", "Sensitivity analysis for effects of K and q", "Sensitivity analysis for algorithmic stability of AdaSub", "Illustrative Example (Appendix)", "PCR data", 
             "High Dimensional regularization methods with cross-validation", "Speed of convergence of AdaSub")
mychoiceM <- menu( choices , graphics=TRUE, title="Which setting \n do you want to run?" )

if (mychoiceM==1) source("./Files/Simulation_LowDimensional.R")
if (mychoiceM==2) source("./Files/Simulation_HighDimensional_with_growing_p.R")
if (mychoiceM==3) source("./Files/Simulation_HighDimensional_choices_K_q.R")
if (mychoiceM==4) source("./Files/Simulation_HighDimensional_stability_K_q.R")
if (mychoiceM==5) source("./Files/Illustrative_Example.R")
if (mychoiceM==6){  
	# Example on PCR data 
	# Please download the PCR data from JRSSB Datasets Vol. 77(5), Song and Liang (2015)
	# Download file from https://wol-prod-cdn.literatumonline.com/pb-assets/hub-assets/rss/Datasets/Vol_77_2015-1521879345620.zip
	# extract  files Xgenes.txt, gene_id.txt & Y3.txt from 77-5.Song.zip inside the zip file
	print( noquote("STEP 1: Please download the PCR data from JRSSB Datasets Vol. 77(5), Song and Liang (2015)" ))
	print( noquote("        from https://wol-prod-cdn.literatumonline.com/pb-assets/hub-assets/rss/Datasets/Vol_77_2015-1521879345620.zip" ))
	print( noquote("STEP 2: Extract  files Xgenes.txt, gene_id.txt & Y3.txt from 77-5.Song.zip from the zip file to the subfolder 'Files' "))
	print( noquote("        (path: '77-5.Song.zip/data_code.zip/data_and_code/file/data/PCR')" ))

	menu("Press 1 when you are ready")
	X = t( read.table("./Files/Xgenes.txt") )
	Y = scan("./Files/Y3.txt")
	Xnames = scan("./Files/gene_id.txt",what=character())
	source("./Files/AdaSub on PCR data.R")
} 
if (mychoiceM==7) source("./Files/Simulation_HighDimensional_with_growing_p_CV.R")
if (mychoiceM==8) source("./Files/Simulation_Speed_of_Convergence.R")


#setwd(savewd)




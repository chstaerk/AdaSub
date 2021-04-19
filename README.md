# AdaSub

Source files for the paper  

Staerk, C., Kateri, M., & Ntzoufras, I. (2021). High-dimensional variable selection via low-dimensional adaptive learning. Electronic Journal of Statistics, 15(1), 830-879.

Please make sure that a version of R >= 3.6.0 is installed.

In order to replicate the results of the paper, please start R, 
set the working directory (setwd(..)) to the unzipped folder
and run the command: 

`source("Master_file.R")`

You have the options to carry out the simulation study in the low-dimensional setting (Section 5.1) as well as in the high-dimensional setting (Section 5.2), the sensitivy analyses (Section 5.3), the illustrative example of AdaSub (Section A2), the PCR data example (Section 6). For line and colour coding of the graphs, please also refer to the paper. 

In order to run the real data example, please download the PCR data from JRSSB Datasets Vol. 77(5), Song and Liang (2015)
from https://wol-prod-cdn.literatumonline.com/pb-assets/hub-assets/rss/Datasets/Vol_77_2015-1521879345620.zip"
and extract the files Xgenes.txt, gene_id.txt & Y3.txt from the zip file 77-5.Song.zip to the subfolder "Files".

Please note that the computation time for a full low- or high-dimensional simulation study are quite long (considering 500 simulated datasets for each sample size). 
In order to obtain faster results, you also have the option to choose a lower number of simulated data examples per sample size
and a lower number of considered different sample sizes. Further note that the computation time for the real data examples are quite long,
since the number of iterations of AdaSub is chosen to be very large (see Section 6). 
In order to obtain faster results, you also have the option to choose a smaller value for the number of iterations.

The file "AdaSub_main_functions.R" includes functions for running the AdaSub algorithm 
and simulating data from a linear regression model with different correlation structures. 
These functions may be used in order to examine further simulated or real data examples with AdaSub. 

In particular, the function 

`AdaSub(...)`

can be used in order to run the AdaSub algorithm. 
For instructions, please see the beginning of the file "AdaSub_main_functions.R" (with description of input and output format).

# AdaSub

Source files for the Adaptive Subspace (AdaSub) method introduced in the paper   

Staerk, C., Kateri, M., & Ntzoufras, I. (2021). High-dimensional variable selection via low-dimensional adaptive learning. Electronic Journal of Statistics, 15(1), 830-879.  https://doi.org/10.1214/21-EJS1797

## Instructions 

Please make sure that a version of R >= 3.6.0 is installed.

In order to replicate the results of the paper, please start R, 
set the working directory `setwd(...)` to the unzipped folder
and run the command: 

`source("Master_file.R")`

You have the options to carry out the simulation study in the low-dimensional setting (Section 5.1) as well as in the high-dimensional setting (Section 5.2), the sensitivity analyses (Section 5.3), the illustrative example of AdaSub (Section A2) and the PCR data example (Section 6). 

In order to run the real data example, please download the PCR data from JRSSB Datasets Vol. 77(5), Song and Liang (2015) from https://wol-prod-cdn.literatumonline.com/pb-assets/hub-assets/rss/Datasets/Vol_77_2015-1521879345620.zip and extract the files `Xgenes.txt`, `gene_id.txt` & `Y3.txt` from the zip file `77-5.Song.zip` to the subfolder `Files`.

Please note that the computation time for a full low- or high-dimensional simulation study are quite long (considering 500 simulated datasets for each sample size). 
In order to obtain faster results, you also have the option to choose a lower number of simulated data examples per sample size and a lower number of considered different sample sizes. Further note that the computation times for the PCR data application are quite long, since the number of iterations of AdaSub is chosen to be very large (see Section 6). In order to obtain faster results, you also have the option to choose a smaller value for the number of iterations.

The file `AdaSub_main_functions.R` includes functions for running the AdaSub algorithm and simulating data from a linear regression model with different correlation structures. These functions may be used in order to examine further simulated or real data examples with AdaSub. 

In particular, the function 

`AdaSub(...)`

can be used in order to run the AdaSub algorithm. 
For instructions, please see the beginning of the file `AdaSub_main_functions.R` (with description of input and output format); for an illustrative application of AdaSub, see the file `Illustrative_Example.R`.  

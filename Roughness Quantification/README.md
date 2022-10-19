This folder contains the python code files that were used to generate, filter, process and analyse the data. 


Dataset_filelocations contains the folder structure of the datafolders of the 2 datasets (original and repeat) and the functions within are used by the other files to read the data in correctly.

Blinds_Core_Functions_0030 contains some essential functions that are used by the other files to help read in, and remove manually selected anomalous frames, as well as plotting functions. It contains the list of all individual frames that have been determined as anomalous. 

Blinds_Stats_and_otherplots_0203 contains statistical functions on the data. This is the main file used to generate the convergence graphs (2 functions are contained for this, with the latter being the latest one). Also contains code for normality/normal tests, finding the confidence overlap, convergence points, and statistical IQR filters on the data. 


Blinds_Method_1_Plotting_0105 contains all functions related to the area roughness data and analysis (Method 1). 

Blinds_Method_2_Plotting_0104 contains all functions related to the line roughness data and analysis (Method 2). 

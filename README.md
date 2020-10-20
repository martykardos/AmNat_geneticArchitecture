# AmNat_geneticArchitecture
This repository contains all code necessary to replicate analyses in “The genegtic architecture of fitness drives population viability during rapid environmental change”, by Marty Kardos and Gordon Luikart.

‘makeFigures’ folder contains an R script (rCode_figures_25June2020.R) to replicate all of the figures in our paper. You will first need to run the simulations for the figures you are interested in replicating (see below). 

‘constantP0’, ‘h2Pt4’, ‘h2Pt6’, ‘h2Pt8’, and ‘mutationSims’ folders contain the code to replicatge the simulations in the paper. To run therse, you can look in rCode_figures_25June2020.R to identify the subfolder where the code for the relevant simulations is located. Note that the directories pointed to in the script rCode_figures_25June2020.R are those located on Marty Kardos’s computer; you will need to change everything up to the /originalSims/ part of the directories to point yourself to the correct place on your own computer. To run a particular set of simualtions, just run the R script(s) located in the appropriate folder.

You can calculate selection limits for particular simulations using the two scripts located in the folder ‘selectionLimits’.

‘RPackages’ folder contains the R packages (quantPop and pedR) necessary to run the simulations are located in the folder ‘RPackages’. You can install these using install.packages() in R. Note that these were developed by Marty Kardos for R versions preceding 4.0. It is possible that the current version of the package will be incompatible with future versions of R.

Direct any questions to martin.kardos@noaa.gov

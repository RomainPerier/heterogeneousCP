Supplementary information / reproducible research files for the manuscript 


Authors: Périer, R., Blanchard, G., Durand, G.
Code was written by Périer, R.
In case of questions or comments please contact romain.perier@universite-paris-saclay.fr !

The code was written/evaluated in R with the following software versions:
R version 4.6.0 (2026-04-24)
Platform: x86_64-pc-linux-gnu
Running under: Linux Mint 21.3

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.10.0 
LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.10.0  LAPACK version 3.10.0

locale:
 [1] LC_CTYPE=fr_FR.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=fr_FR.UTF-8        LC_COLLATE=fr_FR.UTF-8    
 [5] LC_MONETARY=fr_FR.UTF-8    LC_MESSAGES=fr_FR.UTF-8   
 [7] LC_PAPER=fr_FR.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=fr_FR.UTF-8 LC_IDENTIFICATION=C       

time zone: Europe/Paris
tzcode source: system (glibc)

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] pals_1.10         ggplot2_4.0.3     ggmagnify_0.4.2   forcats_1.0.1    
 [5] stringr_1.6.0     purrr_1.2.2       tidyr_1.3.2       dplyr_1.2.1      
 [9] pracma_2.4.6      DiscreteFDR_2.1.1 sanssouci_0.16.2 

loaded via a namespace (and not attached):
 [1] Matrix_1.7-5           gtable_0.3.6           compiler_4.6.0        
 [4] maps_3.4.3             tidyselect_1.2.1       Rcpp_1.1.1-1.1        
 [7] dichromat_2.0-0.1      scales_1.4.0           lattice_0.22-9        
[10] R6_2.6.1               labeling_0.4.3         generics_0.1.4        
[13] mapproj_1.2.12         DiscreteTests_0.4.0    backports_1.5.1       
[16] checkmate_2.3.4        tibble_3.3.1           pillar_1.11.1         
[19] RColorBrewer_1.1-3     rlang_1.2.0            stringi_1.8.7         
[22] S7_0.2.2               cli_3.6.6              withr_3.0.2           
[25] magrittr_2.0.5         matrixTests_0.2.3.1    DiscreteDatasets_0.2.0
[28] grid_4.6.0             rstudioapi_0.18.0      lifecycle_1.0.5       
[31] vctrs_0.7.3            glue_1.8.1             farver_2.1.2          
[34] colorspace_2.1-2       matrixStats_1.5.0      tools_4.6.0           
[37] pkgconfig_2.0.3      

# install necessary packages that are available from Github 
> remotes::install_github("sanssouci-org/sanssouci")


This folder contains the following data and files that can be used to reproduce all analysis and figures of the manuscript.
It contains :
-A master.R file that calls the necessary files to produce figures 2-7 of the manuscript
-Four subfolders containing the following files:

./scripts/:
    calibration_functions.R
    An R script that contains the code for the computations of calibrators functions, shortcuts, etc...
    discrete_functions.R
    An R script that contains the code for some useful functions to generate/consider discrete data, and provides the code for the heterogeneous locally compound p-values.


./config_files/:
    config_shortcuts.R
    An R script that set the parameters of the simulations of the shorcut versions of closure principle based methods.
    config_full_closure.R
    An R script that set the parameters of the simulations of the true closure principle based methods.
   
./simus/
    simu_evalues_shortcuts.R
    An R script that generates the data for figure 1 and save it. The seed is set to 1.
    simu_evalues_shortcuts.R
    An R script that generates the data for simulations based on the true closure principle. The seed is set to 1.
    ./resData/PBD
        A subfolder that contains the data produced by the call of simu_evalues_shortcuts.R.

./figures/:
    script_figure_CE.R
    An R script that is used to generate and save figure 1 in the manuscript. It use simulated data stored in ./../simus/resData/PBD
    

Supplementary information / reproducible research files for the manuscript 


Authors: Périer, R., Blanchard, G., Durand, G.
Code was written by Périer, R.
In case of questions or comments please contact romain.perier@universite-paris-saclay.fr !
R version 4.3.3 (2024-02-29)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Linux Mint 22.3

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.12.0 
LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.12.0

locale:
 [1] LC_CTYPE=fr_FR.UTF-8       LC_NUMERIC=C               LC_TIME=fr_FR.UTF-8        LC_COLLATE=fr_FR.UTF-8     LC_MONETARY=fr_FR.UTF-8   
 [6] LC_MESSAGES=fr_FR.UTF-8    LC_PAPER=fr_FR.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=fr_FR.UTF-8 LC_IDENTIFICATION=C       

time zone: Europe/Paris
tzcode source: system (glibc)

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] eClosure_0.9.6    forcats_1.0.1     stringr_1.6.0     purrr_1.2.2       DiscreteFDR_2.1.1 dplyr_1.2.1       tidyr_1.3.2       ggplot2_4.0.3    

loaded via a namespace (and not attached):
 [1] vctrs_0.7.3            cli_3.6.6              rlang_1.3.0            stringi_1.8.9          generics_0.1.4         S7_0.2.2              
 [7] glue_1.8.1             backports_1.5.1        scales_1.4.0           grid_4.3.3             DiscreteDatasets_0.2.0 tibble_3.3.1          
[13] lifecycle_1.0.5        compiler_4.3.3         DiscreteTests_0.4.1    RColorBrewer_1.1-3     Rcpp_1.1.2             pkgconfig_2.0.3       
[19] farver_2.1.2           R6_2.6.1               dichromat_2.0-1        tidyselect_1.2.1       pillar_1.11.1          magrittr_2.0.5        
[25] checkmate_2.3.4        tools_4.3.3            withr_3.0.3            gtable_0.3.6          

This folder contains the following data and files that can be used to reproduce all analysis and figures of the manuscript.
It contains :
-A master.R file that calls the necessary files to produce figures 1-4 of the manuscript
-Four subfolders containing the following files:

./scripts/:
    calibration_functions.R
    An R script that codes the calibrators used in the paper.
    discrete_functions.R
    An R script that codes the functions to generate contingency tables and the functions related to discretness.
./config_files/:
    config_shortcuts.R
    An R script that sets the parameters of the simulations on the shortcuts.
   
./simus/
    simu_evalues_shortcuts.R
    An R script that generates the data for figure 1,2 and save it. The seed is set to 1.
    ./resData/PBD/
        A subfolder that contains the data produced by the call of simu_shortcut.R, rejection_sets_length_SC.RData.
    simu_evalues_shortcuts2.R
    An R script that generates the data for figure 3,4 and save it. The seed is set to 1.
    ./resData/PBD/
    A subfolder that contains the data produced by the call of simu_shortcut.R, rejection_sets_length_SC2.RData.

./figures/:
    script_figure.R
    An R script that is used to generate and save the figures 1 and 2 in the manuscript. It use simulated data stored in ./../simus/resData/PBD.
    script_figure_2.R
    An R script that is used to generate and save the figures 3 and 4 in the manuscript. It use simulated data stored in ./../simus/resData/PBD.
    

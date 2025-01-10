# importance_qRF

This repository contains the R scripts for running the experiments for importance ranking of modelling factors in quantile random forest-based spatial predictions (qRF predictions).

[] Step 1: dowload the AGB and OCS datasets from [de Bruin et al., 2022](https://doi.org/10.5281/zenodo.6513429).
[] Step 2: unzip utils folder (with scripts for extraction and generation of the spatial data, SIC and random)
[] Step 3: set up subfolders fro storing samples and results

Suggested structures:

|_utils
|_data_EUR
          |__OCS/
          |__AGB/
--> for storing the intermediate results (training and test samples)

./resu_EUR/

--> for storing the CRPS results after applying the random iterations

          |__OCS/
                |_CLUST1/
                |_NCLUST/
                |_RAND/     
                
          |__AGB
                |_CLUST1/
                |_NCLUST/
                |_RAND/   
                
--> for storing the PAWN sensitivity analysis results in the subfolders **CLUST1** (for experiments with one cluster), **NCLUST** (for experiments with mulitple clusters), and **RAND** (for experiments with randomly distributed observations)

The necessary R packages are:
- [ranger](https://doi.org/10.32614/CRAN.package.ranger) for training qRF models
- [terra](https://doi.org/10.32614/CRAN.package.terra), [raster](https://doi.org/10.32614/CRAN.package.raster), and [sf](https://doi.org/10.32614/CRAN.package.sf) for processing the spatial datasets
- [tidyverse](https://doi.org/10.32614/CRAN.package.tidyverse) for generic processing of the data
  

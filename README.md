# importance_qRF

This repository contains the R scripts for running the experiments for importance ranking of modelling factors in quantile random forest-based spatial predictions (qRF predictions) with Sparse, Imprecise, and CLustered observations (SIC).

- [] **Step 1**: dowload the *AGB* and *OCS* datasets from [de Bruin et al., 2022](https://doi.org/10.5281/zenodo.6513429).
- [] **Step 2**: unzip utils folder (with scripts for extraction and generation of the spatial data, SIC and random)
- [] **Step 3**: set up subfolders for storing samples and results

____________________________________________________________________________
Suggested structures:

/utils

/data_EUR/OCS --> for storing the intermediate results (training and test samples)

/resu_EUR/OCS --> for storing the PAWN sensitivity analysis for *OCS* results in the subfolders **CLUST1** (for experiments with one cluster), **NCLUST** (for experiments with mulitple clusters), and **RAND** (for experiments with randomly distributed observations)
____________________________________________________________________________

- [] **Step 4**: run the R scripts where
  - *run_PAWN_RAND_part1.R*  generates randomly distributed samples and compute the CRPS values
  - *run_PAWN_RAND_part2.R*  computes the PAWN sensitivity analysis results from the CRPS values
  - *run_PAWN_oneCluster_part1.R*  generates SIC samples with one cluster and compute the CRPS values
  - *run_PAWN_oneCluster_part2.R*  computes the PAWN sensitivity analysis results from the CRPS values
  - *run_PAWN_multiCluster_part1.R*  generates SIC samples with multiple clusters and compute the CRPS values
  - *run_PAWN_multiCluster_part2.R*  computes the PAWN sensitivity analysis results from the CRPS values

At the begining of each script set up the values of experiment's parameters:
- The choice in the environment variables **CASE0**, either *OCS* or *AGB* 
- Number of iterations: **IT**
- Number of observations: **NN**
- Number of Monte Carlo simulations: **nmc**
- Number of clusters: **Nc**
- Number of observations within the clusters: **N**
- The name of the folders:
  - **infolder** where the raw data are
  - **outfolder** where to store the samples
  - **resfolder** where to store the results

**NOTE** The *CASE.list_NcYY.txt* files contains the list of cluster indices to reproduce the results for YY clusters.

**OUTPUT** PAWN sensitivity measures and p-value of the KS test
- SS.m,PP.m: for absolute error
- SS.ca,PP.ca: for coverage indicator
- SS.crps,PP.crps: for CRPS

The necessary R packages are:
- [ranger](https://doi.org/10.32614/CRAN.package.ranger) for training qRF models
- [terra](https://doi.org/10.32614/CRAN.package.terra), [raster](https://doi.org/10.32614/CRAN.package.raster), and [sf](https://doi.org/10.32614/CRAN.package.sf) for processing the spatial datasets
- [tidyverse](https://doi.org/10.32614/CRAN.package.tidyverse) for generic processing of the data
  

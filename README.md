## Accompanying Code for the Chapter "Latent Gaussian Models and Computation for Large Spatial Data" in the Handbook of MCMC

Authors: Murali Haran, John Hughes, and Ben Seiyon Lee

This repository contains the accompanying code for the book chapter "Latent Gaussian Models and Computation for Large Spatial Data" in the Handbook of MCMC. We provide instructions for implementing the featured approaches for the following three examples: 

Code for the following examples:
1. Hermit Thrush Presence (n=4,937 areal units)
2. MODIS Total Precipitable Water - Massive (n=2,470,059 locations)
3. MODIS Total Precipitable Water - Small (n=20,000 locations)


<p float="left">
  <img src="/MODIS_precipitableWater/TotalPrecipitableWater.png" width="300" />
  <img src="/hermit_thrush_BSF/thrushdata.png" width="500" />
</p>


## Required Packages:
The code has been tested with R version 4.1.1, "Kick Things."  The following R packages must be installed before the code will run successfully:

### Hermit Thrush Example
- `cmdstanr`
- `RSpectra`
- `parallel`

### MODIS Example
- `fields`
- `mvtnorm`
- `classInt`
- `nimble`
- `rstan`

The following programming environments must be installed prior to running the code:
- `nimble`
- `stan`

## Instructions

Before running any code, make sure the required R packages and programming environments (nimble and stan) have been installed.  Set the R working directory to the location of this README file.

### Example 1: Hermit Thrush Example 
- Fit model using stan: `/hermit_thrush_BSF/thrush_stan.R`
- Data including thursh presence, cover, elevation, and locations (x and y): `/hermit_thrush_BSF/thrush_data.csv`
- Adjacency matrix: `/hermit_thrush_BSF/thrush_adjacency.txt`
- Helper functions saved in:
  + `/hermit_thrush_BSF/thrush.stan`

### Preliminaries for Examples 2 and 3 (MODIS Total Precipitable Water)
- Generate samples by running the following files in `/MODIS_precipitableWater/dataProcessing/`
  + `A_MCMCHandbookMODIS.py': Extracts data from hdf file
  + `B_MODISGenerateData.R`: Pre-processes data
  + `C_SampleGeneration_massive.R`: Generates training and test data for massive dataset (n=2,470,059 locations)
  + `C_SampleGeneration_small.R`: Generates training and test data for small dataset (n=20,000 locations)
  + `D_GenerateBisquareBases.R`: Builds bisquare basis functions
  + All datasets are saved in the folder `/MODIS_precipitableWater/samples/`

<p float="left">
  <img src="/MODIS_precipitableWater/source/figures/bisquareBasis.png" width="300" />
  <img src="/MODIS_precipitableWater/source/figures/BS_Basis_1.png" width="300" />
  </p>
  <p float="left">
  <img src="/MODIS_precipitableWater/source/figures/BS_Basis_2.png" width="300" />
  <img src="/MODIS_precipitableWater/source/figures/BS_Basis_3.png" width="300" />
</p>


- Helper files `/MODIS_precipitableWater/source/`
  + `batchmeans.r`: File to compute Montel Carlo standard errors (batch means)
  + `sharedFunctions.R`: Miscellaneous helper functions stored here.
  + `modelStan.stan`: stan file with relevant functions
  + `modelNimble.R`: nimble file with relevant functions

### Examples 2: MODIS Total Precipitable Water - Massive Dataset 
- Fit BR model on massive dataset via nimble
  + Directory with relevant files is `/MODIS_precipitableWater/run_massive/`
  + `A_MCMC.R`: Fits model using nimble
  + `B_analysis.R`: File for validation   

### Examples 3: MODIS Total Precipitable Water - Small Dataset 
- Fit BR model on smaller dataset via nimble, stan, and spNNGP
  + Directory with relevant files is `/MODIS_precipitableWater/run_small/`
  + `A_MCMC.R`: Fits BR model using nimble
  + `A_MCMC_stan.r`: Fits BR model using stan
  + `A_MCMC_spNNGP.R`: Fits SM model using spNNGP
  + `B_analysis.R`: File for validation


 
   

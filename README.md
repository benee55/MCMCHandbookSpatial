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

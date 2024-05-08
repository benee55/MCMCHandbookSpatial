## Accompanying Code for the Chapter "Latent Gaussian Models and Computation for Large Spatial Data" in the Handbook of MCMC

Authors: Murali Haran, John Hughes, and Ben Seiyon Lee


Code for the following examples:
1. MODIS Total Precipitable Water - Massive (n=2,470,059 locations)
2. MODIS Total Precipitable Water - Small (n=20,000 locations)
3. Hermit Thrush Presence (n=4,937 areal units)

<img src="/MODIS_precipitableWater/TotalPrecipitableWater.png" width="700">



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

# Paleoclim_bias_correction

This repository contains scripts for performing **bias correction** of paleoclimate model simulations, and calculating [**bioclim variables**](https://www.worldclim.org/data/bioclim.html) for various ecological studies. 

The workflow of **bias correction** is similar to the studies by [Lima-Ribeiro et al. (2015)](https://doi.org/10.17161/bi.v10i0.4955), [Ficks et al. (2017)](https://doi.org/10.1002/joc.5086), and [Brown et al. (2018)](https://doi.org/10.1038/sdata.2018.254). In general, it includes the following steps:
1. Compute the anomalies between the paleoclimate simulation and the baseline (e.g., 1960-2000 or pre-industry) climate simulation for each variable at its original spatial resolution. 
2. Regrid the anomalies to the desired spatial resolution of the observed baseline climate data. 
3. Apply regridded anomalies to the observed baseline climate data. Similar to the study by [Tang et al. (2011)](https://doi.org/10.5194/cp-7-847-2011), we apply absolute anomalies for daily maximum and minimum temperature. But for precipitation, we apply absolute anomalies when they are positive, and relative anomalies when they are negative.

The bias-corrected monthly averaged daily maximum temperature, daily minimum temperature and precipitation, are used for calculating **bioclim variables** using the “biovars” function from [dismo](https://cran.r-project.org/web/packages/dismo/index.html) package in R.

## Folder content
- **Liu_etal_2023**: scripts for generating bias-corrected climate data and bioclim variables for Early Miocene, Late Miocene and Pliocene, and interpolating the data to the fossil sites in China for redescripition mining. The downstream application can be found [**here**](https://github.com/zliobaite/redescription-asia-neogene) 

## Other relevant github repos


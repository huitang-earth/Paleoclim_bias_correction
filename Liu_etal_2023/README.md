# Scripts for generating Miocien climate data for redescription mining (Liu et al. 2023)

## Software required:
[NCL](https://www.ncl.ucar.edu) is used for general workflow, input/output management and interpolation. 
[ESMF](https://earthsystemmodeling.org) is used for regridding global climate data to different resolution.
[R](https://www.r-project.org/) is used for calculating bioclim variables.

Converting all the scripts into [Python](https://www.python.org/) is under development. 

## Workflow
1. Derive delta-changes (absolute & relative) from paleo-climate model simulations (**prepare_bioclim_###.ncl**). 
2. Regrid delta-changes from the paleoclimate model grids to the desired grids of observation using ESMF. 
    a. Create ScriptGrid file for paleo-climate model grids and observation data grid (**prepare_ScriptGrid_###.ncl**). 
    b. Create mapping file from paleo-climate model grids to observation data grid (using ESMF functionality: [**ESMF_RegridWeightGen**](https://earthsystemmodeling.org/docs/release/ESMF_8_0_1/ESMF_refdoc/node3.html))
    c. Regrid delta-changes from the paleoclimate model grids to the desired grids of observation (**prepare_bioclim_###.ncl**)
3. Apply delta-changes to observation data to derive bias-corrected paleoclimate data(**prepare_bioclim_###.ncl**).
4. Calculate bioclim variables using the bias corrected paleoclimate data (**bioclim_var_###.r**)  
5. Interpolate bioclim data to fossil localities (**prepare_bioclim_###.ncl**).
6. The output data of this workflow can be found [here](https://github.com/huitang-earth/redescription-asia-neogene/tree/main/data_past_org).

## Model data
The information of the paleoclimate data used in this study can be found [here](https://docs.google.com/spreadsheets/d/1BqNGzlA6BiCGrrfwGo2XAtFi_OBoReiTuChnNVq8hgg/edit?usp=sharing). The raw climate data used in these scripts are available either publically or on request to the responsible person. 

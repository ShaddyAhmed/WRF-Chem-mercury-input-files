# WRF-Chem-mercury-input-files
This repository contains scripts and input files to complete the WRF-Chem Arctic mercury run for spring 2020 from Ahmed et al., 2022. 

The WRF-Chem model version used for this run can be downloaded here: https://zenodo.org/record/7137482#.Y3eINuxKh_Q

The description below details the input files are scripts to run (in order) to prepare the WRF-Chem run.

## Running WPS

namelist.wps: The namelist to set up the WRF grid. NOTE: The WRF grid used in the mercury run is domain 2.
1. update_snow_depth.m: Matlab script to set the snow depth, albedo, and sea ice depth to a fixed variable and add to the met_em files.
2. regrid_seaice_depth.ncl: NCL script that regrids sea ice thickness data from the Arctic Ocean Physics Reanalysis dataset to the WRF grid in the met_em files.

## Running real, MEGAN, and WRF-Chem

namelist.input: The namelist to configure the WRF-Chem simulation for the mercury run.
3. regrid_hg_ic_bc_model_ensemble.ncl: NCL script to regrid the Hg(0) initial and boundary conditions to the CAM-Chem grid that will subsequently be regridded to the WRF grid with the mozbc preprocessor.
4. megan_bioemiss.inp: Bio emissions input file to be run with the MEGAN biogenic emissions preprocessor.
5. mozbc_saprc99_mosaic_mercury.inp: Input file to be run with the mozbc initial and boundary conditions preprocessor. This sets the initial and boundary conditions in WRF-Chem for chemical species derived from the CAM-Chem global data set (downloaded at: https://www.acom.ucar.edu/cam-chem/cam-chem.shtml and subset used: Lat = 0 to 90, Lon = 0 to 360, Dates = March 2020 to May 2020).
6. wesely.inp: Input file to be run with the wesely preprocessor.
7. exo_coldens.inp: Input file to be run with exo_coldens preprocessor.
8. fire_emis_saprc.inp: Fire emissions input file to be run with the fire emissions preprocessor. This uses data from the FINNv2.5 data from NCAR (https://rda.ucar.edu/datasets/ds312.9/index.html).
9. prep_anthro_emissions_saprcmosaic_eclipsev6_dms.m: MATLAB script to regrid ECLIPSEv6b anthropogenic emissions (https://previous.iiasa.ac.at/web/home/research/researchPrograms/air/ECLIPSEv6b.html) to the WRF grid. This uses matlab functions in the matlab_functions directory.
10. regrid_hg_emissions.ncl: NCL script to regrid emissions of Hg0, Hg2, HgP to the WRF grid and adds these variables to the wrfchemi_* files.

You are now ready to run WRF-Chem!

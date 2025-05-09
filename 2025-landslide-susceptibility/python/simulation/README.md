# Simulation code
The scripts in this folder are used to run the different simulations for the publication. Note that a separate installation of SlideforMAP is required. 
- __run_sensitivity_analysis_simulations_ha.py__: Runs the simulations for the general sensitivity analysis. Note that this generates a lot of data (>1TB for the 608 plots) and many small files (may be a problem in some filesystems). 
- __run_sensitivity_analysis_simulations_ha_forest.model3__: Runs the simulations for the sensitivity analysis with forest.  
- __run_optimization.py__: Runs the simulations for optimizing the soil physical parameters. The code needs to be modified based on the desired initial values (comment/uncomment the appropriate sections). Each optimization should be run in a separate folder. Otherwise, results are overwritten.
- __run_catchment_simulations.py__: Runs the simulations for the application in the catchments. The configs for the soil physical parameters are found in a subfolder.
- __run_catchment_postprocess.py__: Merges the simulation results from the catchments to absolute landslide probabilities based on diferent correction factors and generates inputs for the final analysis.
- __simulation_parameters.md__: This file documents the parameters used in the simulations.
- __configs__: This folder contains a config.toml file that served as template for the simulations configurations and four CSV files with the soil physical properties used in the catchment applications.
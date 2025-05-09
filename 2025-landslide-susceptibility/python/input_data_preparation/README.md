# Data preparation code
The scripts in this folder prepare the inputs for the simulations. Before running the scripts, the R script for predicting the failure depth has to be run. Then the scripts should be run in the following order:
- __generate_inputs.py__: Prepares the canton wide inputs used by the other scripts. Needs the soil depth raster predicted in R as input. 
- __Prepare_Slope_Stability_Plots.model3__: This QGIS model prepares some of the inputs used by the scripts to select and prepare the hectare plots and catchments. The process sources are enriched by the region and statistics on slides, slope, and soil type. The catchments are enriched by slope statistics. The model also gneerates the hectare plots based on the clistered slides in the inventors as well as the regular hectare grid.
- __prepare_inputs_ha.py__: Selects the hectare plots used for the sensitivity analysis and optimization. Inputs per plot for the simulations in SlideforMAP are prepared. 
- __prepare_inputs_catchments.py__: Prepares the inputs per catchment for the simulations in SlideforMAP. 

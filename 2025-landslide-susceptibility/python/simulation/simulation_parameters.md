# Simulation parameters

## General parameters for SlideforMAP simulations

The following general configuration parameters are used in the TOML configuration files for the SlideforMAP simulations, 

|**Parameter**  | **Value** | **Description**| 
|:---|:---|:---|
|demPath   | individual DEM raster for plot or catchment | Path to DEM Raster |
|flowAccPath   | individual flow accumulation raster for plot or catchment | Path to flow accumulation Raster|
|soilTypePath  | individual soil type raster for plot or catchment | Path to soil type Raster|
|soilThicknessPath  | Not used | Path to soil depth CSV|
|physPropPath  | CSV file with soil physical property values matching the soil types based on the respective simulation inputs | Path to the Physical Properties CSV|
|soilThicknessValueRaster | individual raster with predicted failure depth for plot or catchment | Path to soil thickness values raster|
|soilThicknessSDPercent | 0 | Soil thickness standard deviation in percent |
|landUseRaster | Not used | Path to land use raster. Leave empty if not available|
|landUseCSV | Not used | Path to land use CSV. Leave empty if land use raster not used|
|treeFileFormat  | 1 | Partial data (1) or full data (2)|
|treeFilePath | individual treefile for plot or catchment | Path to the tree data set|
|treeSpecies | "Picea abies" of "Fagus sylvatica" depending on the mixing degree in the plot or catchment | Name of tree species|
|constantLateralRR | 0 | Value of constant lateral root reinforcement (N/m)|
|outputPath | Depending on the simulation | Path to output directory|
|outputExtension | .tif | File extension for output rasters|
|outputTreeWeight | 0 | 0 no output; 1 output raster|
|outputSlope  | 1 | 0 no output; 1 output raster|
|rainfall | depending on the simulation inputs& Rainfall rate [mm/hr]|
|minLandslideArea | 10 | Minimum landslide area [m^2]|
|maxLandslideArea | 2600 | Maximum landslide area [m^2]|
|safetyFactorWithSoilCompression  | 1 | Compute safety factor with soil compression|
|exportLandslides | 1 for sensitivity analysis, 0 for optimisation and application | Export landslide data to file|
|runLargeWoodModel | 0 | Compute large wood recruitment|
|runSedimentModel | 0 | Compute sediment recruitment|
|rasterizerWidth | 4 | Half width over which to average landslide for raster creation (cells) [m]|
|landslideDensity | 1 | Number of landslides per square meter |
|streamSearchDistance | 200 | Maximum distance from landslide to stream for runout [m]|
|catchmentAreaThreshold  | 50000 | Minimum value of catchment area for stream search [m^2]|
|treeSearchRadius | 7 | Search radius to find nearby trees of raster cell [m]|
|slopeThreshold  | 50 | Threshold slope above which soil cohesion is set to 100 kPa [deg]|
|runoutAngle  | 23 | Value of angle used in Weibull distribution (scale factor) [deg]|
|saturatedThicknessFraction | -1 | Fraction of soil thickness that is water saturated [-]|
|passiveEarthPressureCorrection | 30 | Passive earth pressure correction in [%]|
|rainCorrectionFactor | 1 | Rainfall rate correction factor for small landslides. Set to 1 for no correction [-]|
|rainCatchmentAreaThreshold | 1000 | Upper threshold catchment area for rainfall correction [m^2]|
|invGammaShape | 0.784707594 | Landslide area inverse Gamma distribution shape factor [-] (Malamud  et al. 2004), Eq (3) variable rho|
|invGammaScale | 0.000170460 | Landslide area inverse Gamma distribution scale factor [km^2] (Malamud  et al. 2004), Eq (3) variable a|
|invGammaS |     2.11826309e-30 | Landslide area inverse Gamma distribution s parameter [km^2] (Malamud  et al. 2004), Eq (3) variable s|

## Soil physical properties
The following table shows the parameters for the physical properties of the different soil texture classes. The *Expert* columns contain default values based on expert opinion. The *Min* and *Max* columns indicate the minimum and maximum bounds of the value range applied during optimisation.

|**Soil ID** | **Soil Name** | **Porosity**| **Field capacity**| **Friction angle mean [deg]**|**Friction angle mean [deg]**|**Friction angle mean [deg]**| **Friction angle std dev [deg]**| **Cohesion mean [kPa]**|**Cohesion mean [kPa]**|**Cohesion mean [kPa]**| **Cohesion std dev [kPa]**| **Ks mean [m/d]** |**Ks mean [m/d]**|**Ks mean [m/d]**| **Ks std dev [m/d]** |
|---:|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:| 
| |  | *Expert* |  *Expert* |  *Expert* |  *Min* |  *Max* |  *Expert* |  *Expert* |  *Min* |  *Max* |  *Expert* |  *Expert* |  *Min* |  *Max* |  *Expert*|
|1 | Sand > SP | 0.38 | 0.10 | 36 |  |  | 0 | 0 |  |  | 0 | 8600 |  |  | 0 | 
|2 | Loamy Sand > SM | 0.37 | 0.16 | 34 |  |  | 0 | 0 |  |  | 0 | 20 |  |  | 0 | 
|3 | Sandy Loam > SM | 0.37 | 0.21 | 34 | 29 | 34 | 0 | 0.5 | 0 | 5 | 0 | 20 | 10 | 100 | 0 | 
|4 | Sandy Clay Loam > SC | 0.40 | 0.36 | 32 | 28 | 36 | 0 | 1.5 | 0 | 5 | 0 | 20 | 0 | 80 | 0 | 
|5 | Sandy Clay > SC | 0.40 | 0.32 | 32 |  |  | 0 | 0 |  |  | 0 | 20 |  |  | 0 | 
|6 | Loam > CL | 0.41 | 0.27 | 27 |  |  | 0 | 2 |  |  | 0 | 20 |  |  | 0 | 
|7 | Silt Loam > ML | 0.47 | 0.3 | 33 | 28 | 34 | 0 | 0 | 0 | 5 | 0 | 20 | 0 | 80 | 0 | 
|8 | Silt > ML | 0.47 | 0.29 | 33 | 29 | 37 | 0 | 0 | 0 | 5 | 0 | 20 | 0 | 80 | 0 | 
|9 | Clay Loam > CL | 0.41 | 0.29 | 27 |  |  | 0 | 2 |  |  | 0 | 20 |  |  | 0 | 
|10 | Silty Clay Loam > CL | 0.41 | 0.28 | 27 |  |  | 0 | 2 |  |  | 0 | 20 |  |  | 0 | 
|11 | Clay > CH | 0.56 | 0.22 | 22 |  |  | 0 | 2.5 |  |  | 0 | 20 |  |  | 0 | 
|12 | Silty Clay > CH | 0.56 | 0.20 | 22 |  |  | 0 | 2.5 |  |  | 0 | 20 |  |  | 0 | 
|13 | NoSlide | 0 | 0.17 | 90 |  |  | 0 | 999 |  |  | 0 | 20000 |  |  | 0 | 
|21 | Sand > SP in Forest | 0.38 | 0.1 | 36 |  |  | 0 | 0 |  |  | 0 | 9300 |  |  | 0 | 
|22 | Loamy Sand > SM in Forest | 0.37 | 0.16 | 34 |  |  | 0 | 0 |  |  | 0 | 720 |  |  | 0 | 
|23 | Sandy Loam > SM in Forest | 0.37 | 0.21 | 34\* |  |  | 0 | 0.5\* |  |  | 0 | 720\* |  |  | 0 | 
|24 | Sandy Clay Loam > SC in Forest | 0.40 | 0.36 | 32\* |  |  | 0 | 1.5\* |  |  | 0 | 720\* |  |  | 0 | 
|25 | Sandy Clay > SC in Forest | 0.40 | 0.32 | 32 |  |  | 0 | 0 |  |  | 0 | 720 |  |  | 0 | 
|26 | Loam > CL in Forest | 0.41 | 0.27 | 27 |  |  | 0 | 2 |  |  | 1 | 720 |  |  | 0 | 
|27 | Silt Loam > ML in Forest | 0.47 | 0.30 | 33\* |  |  | 0 | 0\* |  |  | 0 | 720\* |  |  | 0 | 
|28 | Silt > ML in Forest | 0.47 | 0.29 | 33\* |  |  | 0 | 0\* |  |  | 0 | 720\* |  |  | 0 | 
|29 | Clay Loam > CL in Forest | 0.41 | 0.29 | 27 |  |  | 0 | 2 |  |  | 1 | 720 |  |  | 0 | 
|30 | Silty Clay Loam > CL in Forest | 0.41 | 0.28 | 27 |  |  | 0 | 2 |  |  | 1 | 720 |  |  | 0 | 
|31 | Clay > CH in Forest | 0.56 | 0.22 | 22 |  |  | 0 | 2.5 |  |  | 1 | 720 |  |  | 0 | 
|32 | Silty Clay > CH in Forest | 0.56 | 0.2 | 22 |  |  | 0 | 2.5 |  |  | 1 | 720 |  |  | 0 | 
|33 | NoSlide in Forest | 0 | 0.17 | 90 |  |  | 0 | 999 |  |  | 0 | 20000 |  |  | 0 | 
\* Parameter value is coupled with corresponding soil type outside forest.

The following table shows the parameter values for mean friction angle, mean cohesion, and mean hydraulic conductivity (Ks) resulting from selected optimisation runs. 

| **Optimisation run** | **Soil ID** | **Friction angle mean (deg)** | **Cohesion mean (kPa)** | **Ks mean (m/d)** | 
|:---|---:|---:|---:|---:|
|A) Minimum initial, epsfcn=0.1 | 3 | 32.0830 | 1.6398 | 11.4039 |
| | 4 | 28.2739 | 1.6514 | 8.1972 |
| | 7 | 28.9375 | 2.0803 | 63.1700 |
| | 8 | 30.4841 | 2.0771 | 2.1287 | 
|B) Maximum initial, epsfcn=0.25 | 3 | 33.8007 | 4.8627 | 99.8485 |
| | 4 | 31.3516 | 1.1501 | 8.9047 |
| | 7 | 34.0000 | 5.0000 | 80.0000 |
| | 8 | 35.5693 | 0.4980 | 28.8937 | 
|C) Middle initial, epsfcn=0.25 | 3 | 30.6080 | 4.3540 | 29.9669 |
| | 4 | 32.2517 | 2.9184 | 26.6446 |
| | 7 | 28.5362 | 1.1349 | 0.4471 |
| | 8 | 32.9823 | 2.3055 | 26.0249 | 
|D) Expert input, epsfcn=0.1 | 3 | 29.8969 | 1.7703 | 20.5889 |
| | 4 | 31.4550 | 1.7736 | 14.9045 |
| | 7 | 28.8020 | 0.0617 | 0.4380 |
| | 8 | 32.2398 | 2.3383 | 4.6839 | 
|E) Expert input initial values | 3 | 34.0000 | 0.5000 | 20.0000 |
| | 4 | 32.0000 | 1.5000 | 20.0000 |
| | 7 | 33.0000 | 0.0000 | 20.0000 |
| | 8 | 33.0000 | 0.0000 | 20.0000 |

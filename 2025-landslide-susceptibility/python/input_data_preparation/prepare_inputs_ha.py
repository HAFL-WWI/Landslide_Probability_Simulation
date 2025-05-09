######################################################################
# Copyright (C) 2025 BFH
#
# Script for preparing the inpurs for the SlideforMAP simulations 
# per plot for the sensitivity analysis and the optimization. 
# Requires the canton wide inputs and the outputs of the 
# Prepare_Slope_Stability_Plots.model3 QGIS model (catchmants,
# process sources, plots) as basis.
#
# Author: Christoph Schaller, BFH-HAFL, March 2025
######################################################################

#
# Imports
#

from osgeo import gdal
from osgeo.gdalconst import *
from osgeo.gdalnumeric import *
from osgeo_utils import gdal_calc, gdal_polygonize
from osgeo import gdal_array
gdal.UseExceptions()

import rasterio as rio

from shapely import geometry

import numpy as np
import pandas as pd
import geopandas as gpd
import os
import sys
import math

import tomlkit

from SALib.sample.sobol import sample

from sklearn.model_selection import train_test_split

# Adding the path to the utilities.py to the path and import
utility_path = "E:/GIS_Projekte/Paper_3/code/utilities"
sys.path.append(utility_path)
import utilities

#
# Configs and paths
#

slide_path = "E:/GIS_Projekte/Paper_3/data/BE_SL_all_2025_T01/BE_SL_all_2025_T01.shp"
slide_path = "E:/GIS_Projekte/Paper_3/data/BE_SL_all_2025_T01/slides_with_soiltype.gpkg"
slide_layer = "slides"

hades_base_path = "E:/GIS_Projekte/Paper_2/data/HADES_B04"
hades_files = {
      2: "xspace_data_for_hades_24h_2J_0975_.tif",
     10: "xspace_data_for_hades_24h_10J_0975_.tif",
     30: "xspace_data_for_hades_24h_30J_0975_.tif",
    100: "xspace_data_for_hades_24h_100J_0975_.tif",
    300: "xspace_data_for_hades_24h_300J_0975_.tif",
    }

runoff_coefficients = {
      2: 0.10,
     10: 0.12,
     30: 0.15,
    100: 0.20,
    300: 0.30,
    }

fintch_db_path = "F:/fint-ch/AP10__Validierung_der_EBD_Analysen/FINT_v2_BE.gpkg"

raster_resolution = 2

canton = "BE"

input_base_path = "E:/GIS_Projekte/Paper_3/data/slide_inputs"
output_base_path = "E:/GIS_Projekte/Paper_3/data/optimization_inputs"
plot_name_template = "plot_{0}"

perimeter_buffer = 50

#
# Sample parameter combinations for sensitivity analysis
#

utilities.ensure_dir(output_base_path)

# Combinations for the general analysis

# Create a problem dictionary. 
# Here we supply the number of variables, the names of each variable, and the bounds of the variables.
problem = {
    "num_vars": 6,
    "names": ["c_mean", "fc_mean", "ks_mean", "por_mean", "r_mean","rainfall"],
    "bounds": [[0, 5], # cohesion
               [0.1, 0.4], # field capacity
               [0.000864, 8640], #transmissivity ks
               [0.25,0.75], # porosity
               [20, 45], # internal angle of friction r
               [0.1, 5], # maybe add significantly higher maximum 
               ]
}

# Generate parmeter values using the sobol.sample function (Saltelli algo)
param_values = sample(problem, 1024)

# Save combinations to textfile
param_out_path = os.path.join(output_base_path,"sensitivity_param_values.txt")
np.savetxt(param_out_path, param_values)

# Combinations for the ganalysis with forest

#Create a problem dictionary. 
problem = {
    "num_vars": 7,
    "names": ["c_mean", "fc_mean", "ks_mean", "por_mean", "r_mean","rainfall","forest"],
    "bounds": [[0, 5], # cohesion
               [0.1, 0.4], # field capacity
               [0.000864, 8640], #transmissivity ks
               [0.25,0.75], # porosity
               [20, 45], # internal angle of friction r
               [0.1, 5], # maybe add significantly higher maximum 
               [0,1], # simulate without (<0.5) or with (>=0.5) forest
               ]
}
# Generate parmeter values using the sobol.sample function (Saltelli algo)
param_values = sample(problem, 1024)

# Round forest variable to get 0/1 values
param_values[:,6] = np.round(param_values[:,6],0)

# Save combinations to textfile
param_out_path = os.path.join(output_base_path,"sensitivity_param_values_with_forest.txt")
np.savetxt(param_out_path, param_values)


#
# Plot preparation functions
#

# Clip a raster and multiply with a 0/1 mask
def clip_and_mask(in_raster,out_raster,mask_raster,bounds):
    tmp_raster = out_raster.replace(".tif","_tmp.tif")
    with gdal.Warp(destNameOrDestDS=tmp_raster,srcDSOrSrcDSTab=in_raster, format="GTiff", creationOptions={"COMPRESS":"DEFLATE", "BIGTIFF":"YES"}, outputBounds=bounds) as ds:
        ds=None

    with gdal_calc.Calc(calc="A*B",A=tmp_raster, B=mask_raster, outfile=out_raster, creation_options = ["COMPRESS=DEFLATE", "BIGTIFF=YES"] , overwrite=True) as ds:
        ds=None
    
    utilities.delete_raster(tmp_raster)

# Prepare the inputs for the specified plots
def prepare_inputs(plot_ids,plots_gpkg,point_buff_layer,sensitivity=False):

    # Prepare the hades rasters for sampling by creating the source objects
    hades_rasters = {}
    for k in hades_files:
        hades_raster_path = os.path.join(hades_base_path,hades_files[k])
        src = rio.open(hades_raster_path, "r")
        hades_rasters[k] = src

    # Read all plots
    all_plots_df = gpd.read_file(plots_gpkg, layer="plots", engine="pyogrio", fid_as_index=True)
    all_plots_df.reset_index(inplace=True)

    # Process each selected plot
    for fid in plot_ids:
        # Get the process source id
        ps_id = int(all_plots_df.loc[all_plots_df["fid"]==fid].iloc[0]["ps_psid"])

        # Create output folder for plot
        plot_folder = plot_name_template.format(ps_id)
        plot_path = os.path.join(output_base_path,plot_folder)
        utilities.ensure_dir(plot_path)

        # Create temporary geopackage containing only the current plot
        perimeter_tmp_gpkg = os.path.join(plot_path,"perimeter_tmp.gpkg")
        perimeter_gpkg = os.path.join(plot_path,"perimeter.gpkg")
        perimeter_layer = "perimeter"
        perimeter_buffered_layer = "perimeter_buffered"
        with gdal.VectorTranslate(destNameOrDestDS=perimeter_tmp_gpkg,  srcDS=plots_gpkg, layerName=perimeter_layer, SQLStatement="SELECT * FROM {0} WHERE fid={1}".format(ps_layer,fid),  accessMode="overwrite") as ds:
            ds = None

        # Read plot geometry and get extent
        plot_df = gpd.read_file(perimeter_tmp_gpkg, layer=perimeter_layer, engine="pyogrio", fid_=True)
        extent_gpkg = utilities.get_gpkg_extent(perimeter_tmp_gpkg,perimeter_layer)
        utilities.delete_gpkg(perimeter_tmp_gpkg)

        # Determine plot center coordinates
        x_center = extent_gpkg[0]+(extent_gpkg[1]-extent_gpkg[0])/2
        y_center = extent_gpkg[2]+(extent_gpkg[3]-extent_gpkg[2])/2

        # Round the extent to resolution and buffer by 1 pixel
        x_min = raster_resolution*(math.floor(extent_gpkg[0]/raster_resolution)) - raster_resolution
        x_max = raster_resolution*(math.ceil(extent_gpkg[1]/raster_resolution)) + raster_resolution
        y_min = raster_resolution*(math.floor(extent_gpkg[2]/raster_resolution)) - raster_resolution
        y_max = raster_resolution*(math.ceil(extent_gpkg[3]/raster_resolution)) + raster_resolution

        # Write perimeter geopackage from rounded extent
        perimeter_geom = geometry.Polygon(((x_min,y_min), (x_min,y_max), (x_max,y_max), (x_max,y_min)))
        df = pd.DataFrame({"id": [1], "geom": [perimeter_geom]})
        gdf = gpd.GeoDataFrame(df, geometry="geom", crs="EPSG:2056")
        gdf.to_file(perimeter_gpkg, layer=perimeter_layer, driver="GPKG")

        # Copy slides in plot to geopackage
        with gdal.VectorTranslate(destNameOrDestDS=perimeter_gpkg, srcDS=slide_path, layerName=slide_layer, spatFilter=[x_min,y_min,x_max,y_max], accessMode="overwrite") as ds:
            ds = None

        # Buffer extent by 50m
        x_min -= perimeter_buffer
        y_min -= perimeter_buffer
        x_max += perimeter_buffer
        y_max += perimeter_buffer

        # Create mask of normal perimeter at 2m resolution
        perimeter_mask_path = os.path.join(plot_path,"perimeter_mask.tif")
        with gdal.Rasterize(destNameOrDestDS=perimeter_mask_path, srcDS=perimeter_gpkg, layers=perimeter_layer, format="GTiff", burnValues = 1, initValues=-9999, noData =-9999, outputType=gdalconst.GDT_Int16, creationOptions={"COMPRESS":"DEFLATE", "BIGTIFF":"YES"}, xRes=raster_resolution, yRes=raster_resolution, outputBounds=[x_min,y_min,x_max,y_max]) as ds:
            ds=None

        # Get numbers of cells in perimeter
        perimeter_mask_array = gdal_array.LoadFile(perimeter_mask_path)
        perimeter_cell_count = np.sum(perimeter_mask_array[perimeter_mask_array!=-9999])

        # Add layer with buffered perimeter to geopackage
        perimeter_buffered_geom = geometry.Polygon(((x_min,y_min), (x_min,y_max), (x_max,y_max), (x_max,y_min)))
        df = pd.DataFrame({"id": [1], "geom": [perimeter_buffered_geom]})
        gdf = gpd.GeoDataFrame(df, geometry="geom", crs="EPSG:2056")
        gdf.to_file(perimeter_gpkg, layer=perimeter_buffered_layer, driver="GPKG")

        # Create mask of buffered perimeter at 2m resolution. This is the basis for cutting the simulation inputs.
        perimeter_mask_buffered_path = os.path.join(plot_path,"perimeter_buffered_mask.tif")
        with gdal.Rasterize(destNameOrDestDS=perimeter_mask_buffered_path, srcDS=perimeter_gpkg, layers=perimeter_buffered_layer, format="GTiff", burnValues = 1, initValues=-9999, noData =-9999, outputType=gdalconst.GDT_Int16, creationOptions={"COMPRESS":"DEFLATE", "BIGTIFF":"YES"}, xRes=raster_resolution, yRes=raster_resolution, outputBounds=[x_min,y_min,x_max,y_max]) as ds:
            ds=None

        # Get number of cells in buffered perimeter
        perimeter_buffered_mask_array = gdal_array.LoadFile(perimeter_mask_buffered_path)
        perimeter_buffered_cell_count = np.sum(perimeter_buffered_mask_array[perimeter_buffered_mask_array!=-9999])

        # Clip and mask the DEM
        dem_2m_path = os.path.join(input_base_path,"dem_2_LV95_BE_filled.tif") 
        dem_path = os.path.join(plot_path, "dem.tif")
        clip_and_mask(dem_2m_path,dem_path,perimeter_mask_buffered_path,[x_min,y_min,x_max,y_max])

        # Clip and mask the contributing area
        contributing_area_2m_name = "contributing_area_"+canton
        contributing_area_2m_path = os.path.join(input_base_path,contributing_area_2m_name+".sdat")
        contributing_area_path = os.path.join(plot_path,"contributing_area.tif")
        with gdal.Warp(destNameOrDestDS=contributing_area_path,srcDSOrSrcDSTab=contributing_area_2m_path, format="GTiff", creationOptions={"COMPRESS":"DEFLATE", "BIGTIFF":"YES"}, outputBounds=[x_min,y_min,x_max,y_max]) as ds:
            ds=None

        # Clip and mask the failure depth
        soil_depth_2m_path = os.path.join(input_base_path,"soil_depth.tif")
        soil_depth_path = os.path.join(plot_path,"soil_depth.tif")
        clip_and_mask(soil_depth_2m_path,soil_depth_path,perimeter_mask_buffered_path,[x_min,y_min,x_max,y_max])
        
        # Clip and mask the constant soil type for sensitivity analysis
        soil_type_name = "soil_type_"+canton
        soil_type_constant_2m_path = os.path.join(input_base_path,soil_type_name+".tif")
        soil_type_constant_path = os.path.join(plot_path,"soil_type_constant.tif")
        clip_and_mask(soil_type_constant_2m_path,soil_type_constant_path,perimeter_mask_buffered_path,[x_min,y_min,x_max,y_max])

        # Clip and mask the "USDA" soil type for optimization 
        soil_type_slide_2m_path = os.path.join(input_base_path,"soil_type_slide.tif")
        soil_type_slide_path = os.path.join(plot_path,"soil_type_slide.tif")
        clip_and_mask(soil_type_slide_2m_path,soil_type_slide_path,perimeter_mask_buffered_path,[x_min,y_min,x_max,y_max])
        
        # Clip the mask with the unbuffered slides
        slide_mask_2m_path = os.path.join(input_base_path,"slide_mask.tif") 
        slide_mask_path = os.path.join(plot_path,"slide_mask.tif")
        clip_and_mask(slide_mask_2m_path,slide_mask_path,perimeter_mask_path,[x_min,y_min,x_max,y_max])
   
        # Clip the mask with the 5m buffered slides
        slide_buff_mask_path = os.path.join(plot_path,"slide_buff{0}_mask.tif".format(buffer_dist))
        clip_and_mask(point_buff_layer,slide_buff_mask_path,perimeter_mask_path,[x_min,y_min,x_max,y_max])

        # Clip and mask the tree type raster
        tree_type_2m_name = "tree_type_"+canton
        tree_type_2m_path = os.path.join(input_base_path,tree_type_2m_name+".tif")
        tree_type_path = os.path.join(plot_path,"tree_type.tif")
        clip_and_mask(tree_type_2m_path,tree_type_path,perimeter_mask_buffered_path,[x_min,y_min,x_max,y_max])
        
        # Determine the majority of the treetype to decide the species for the simulation
        tree_type_array = gdal_array.LoadFile(tree_type_path)
        fagus_count = np.sum(tree_type_array==2)
        picea_count = np.sum(tree_type_array==1)
        tree_type = "Picea abies" if picea_count>fagus_count else "Fagus sylvatica"

        # Generate a treefile for the trees within the plot
        treefile_gpk_tmp = os.path.join(plot_path, "treefile_tmp.gpkg")
        treefile_path = os.path.join(plot_path,"treefile.txt")
        tree_layer = "tree_layer"
        with gdal.VectorTranslate(destNameOrDestDS=treefile_gpk_tmp, srcDS=fintch_db_path, format="GPKG", layerName=tree_layer, layers="be_processed_tree", accessMode="overwrite", spatFilter=[x_min,y_min,x_max,y_max]) as ds :
            ds = None

        df_out = gpd.GeoDataFrame.from_file(treefile_gpk_tmp,layer=tree_layer,engine="pyogrio",)

        x_min = raster_resolution*(math.floor(extent_gpkg[0]/raster_resolution)) - raster_resolution
        x_max = raster_resolution*(math.ceil(extent_gpkg[1]/raster_resolution)) + raster_resolution
        y_min = raster_resolution*(math.floor(extent_gpkg[2]/raster_resolution)) - raster_resolution
        y_max = raster_resolution*(math.ceil(extent_gpkg[3]/raster_resolution)) + raster_resolution

        # Write treefile if there are any trees
        has_forest = False
        forest_cell_count = 0
        if len(df_out)>0:
            has_forest = True
            with open(treefile_path,"w+") as ofile:
                fmt = "{:>15.3f}{:>15.3f}{:>15.3f}"
                fmt = "%15.3f %15.3f %15.3f"
                np.savetxt(ofile, df_out[["x","y","bhd"]].values, fmt=fmt) #FINT-CH
    
            # Write forest mask from 10m buffered trees
            tree_buff_dist = 10
            tree_buff_layer = "trees_buffered"
            with gdal.VectorTranslate(destNameOrDestDS=treefile_gpk_tmp, srcDS=treefile_gpk_tmp, format="GPKG", layerName=tree_buff_layer, SQLStatement="select ST_buffer(geom, {0}) as geom FROM {1}".format(tree_buff_dist,tree_layer), SQLDialect="SQLite", accessMode="overwrite") as ds :
                ds = None

            forest_mask_rast = os.path.join(plot_path, "forest_mask.tif")
            with gdal.Rasterize(destNameOrDestDS=forest_mask_rast, srcDS=treefile_gpk_tmp, layers=tree_buff_layer, format="GTiff", outputType=gdalconst.GDT_Byte, creationOptions={"COMPRESS":"DEFLATE", "BIGTIFF":"YES"}, xRes=raster_resolution, yRes=raster_resolution, outputBounds=[x_min,y_min,x_max,y_max], burnValues=1, initValues=0, noData=0) as ds:
                ds=None

            # Get number of forest cells to determine forest cover of plot
            forest_mask_array = gdal_array.LoadFile(forest_mask_rast)
            forest_cell_count = np.sum(forest_mask_array[forest_mask_array==1])

            # Generate polygon version of forest mask
            forest_mask_gpkg = os.path.join(plot_path, "forest_mask.gpkg")
            forest_mask_layer = "forest_mask"
            gdal_polygonize.gdal_polygonize(src_filename=forest_mask_rast, band_number=1, dst_filename=forest_mask_gpkg, dst_layername=forest_mask_layer, driver_name="GPKG")
            utilities.delete_raster(forest_mask_rast)
        utilities.delete_gpkg(treefile_gpk_tmp)


        # Determine percipitation from HADES
        coord_list = [(x_center, y_center)]

        perc_amounts = {}
        for k in hades_files:
            src = hades_rasters[k]
            perc_sample = [x for x in src.sample(coord_list)]

            perc = round(float(perc_sample[0]/24*runoff_coefficients[k]),4)
            perc_amounts["perc_{0}".format(k)] = perc
            
        # Write config file by modifying the "default_config.toml" template
        toml_template_path = os.path.join(input_base_path,"default_config.toml")
        toml_output_path = os.path.join(plot_path,"default_config.toml")
        with open(toml_template_path, mode="rt") as fp:
            config = tomlkit.load(fp)

            config["plot_sensitivity"] = sensitivity

            config["perimeter_cell_count"] = int(perimeter_cell_count)
            config["perimeter_buffered_cell_count"] = int(perimeter_buffered_cell_count)
            config["forest_cell_count"] = int(forest_cell_count)
            config["x_min"] = x_min
            config["y_min"] = y_min
            config["x_max"] = x_max
            config["y_max"] = y_max

            for k in perc_amounts:
                config[k] = perc_amounts[k]

            config["ps_slide_count"] = plot_df.iloc[0]["ps_slide_count"]
            config["catchment_slide_count"] = plot_df.iloc[0]["ps_ezg_slide_count"]
            config["catchment_slide_per_ha"] = plot_df.iloc[0]["ps_ezg_slides_per_ha"]
            config["slide_count"] = int(plot_df.iloc[0]["slide_count"])

            config["slide_in_forest_count"] = int(plot_df.iloc[0]["slide_in_forest_count"]) if "slide_in_forest_count" in plot_df.columns else 0

            config["demPath"] = "../{0}".format(os.path.basename(dem_path))
            config["flowAccPath"] = "../{0}".format(os.path.basename(contributing_area_path))
            config["soilTypeRaster"] =  "../{0}".format(os.path.basename(soil_type_constant_path))    
            config["soilThicknessValueRaster"] = "../{0}".format(os.path.basename(soil_depth_path)) 

            # Don't write this information for plots without forest
            if has_forest:
                config["treeFilePath"] = "../{0}".format(os.path.basename(treefile_path))   
                config["treeSpecies"] = tree_type

            with open(toml_output_path, mode="wt") as fp:
                tomlkit.dump(config, fp)


#
# Sample plots with slides for sensitivity analysis
#

# Read process sources that were enriched with the region, statistics about the slides, slope at 10m, and soil type 
ps_gpkg = os.path.join(input_base_path,"process_sources_with_slides.gpkg")
ps_layer = "plots"
ps_df = gpd.read_file(ps_gpkg, layer=ps_layer, engine="pyogrio", fid_as_index=True)

# Read hectare plots around slide clusters 
plot_gpkg = os.path.join(input_base_path,"plots_with_slides.gpkg")
plot_layer = "plots"
plots_df = gpd.read_file(plot_gpkg, layer=plot_layer, engine="pyogrio", fid_as_index=True)

# Read hectare plots for regular grid
grid_gpkg = os.path.join(input_base_path,"grid_plots.gpkg")
grid_layer = "plots"
grid_df = gpd.read_file(grid_gpkg, layer=grid_layer, engine="pyogrio", fid_as_index=True)

# Parameters for selecting the plots
fraction = 0.1 # Percentage of all process sources with slides
random_seed = 42
min_slope_limit = 15

# Add index/FID as field
plots_df.reset_index(inplace=True)
ps_df.reset_index(inplace=True)

# Add boolean flag for observed slides
ps_df["has_slides"] = ps_df["slide_count"]>0
# Add group label based on region (tecto) and slide flag
ps_df["group"] = ps_df["tecto"]+ps_df["has_slides"].astype("str")

#Include only plots that contain at least one pixel with slope_10m>=15°
ps_df = ps_df.loc[ps_df["slope10_max"]>=min_slope_limit]

ps_ids = []

# Get inputs for selection only from process sources with observed slides
X = ps_df.loc[ps_df["has_slides"]]["tecto"]
y = ps_df.loc[ps_df["has_slides"]]["fid"]
groups = ps_df.loc[ps_df["has_slides"]]["group"]

# Do selection through stratified train test split
X_train, X_test, y_train, y_test= train_test_split(X,y, train_size=fraction, random_state=random_seed, stratify=groups)

# Get IDs of selected plots
ps_ids = y_train.values

# Write a geopackage containing the selected process sources
ps_sensitivity_df = ps_df[ps_df["fid"].isin(ps_ids)]
ps_sensitivity_gpkg = os.path.join(input_base_path,"ps_sensitivity.gpkg")
ps_sensitivity_df.to_file(ps_sensitivity_gpkg, layer=ps_layer, driver="GPKG")

# For all selected process sources, select the plot with the most observed slides that is above the slope_10m>15§ limit
plot_ids = []
for fid in ps_ids:
    psid = ps_df.loc[ps_df["fid"]==fid]["psid"].iloc[0]
    sel = plots_df.loc[plots_df["ps_psid"]==psid]
    plot_id = None

    # Check plots for limit by descending slide_count
    for index, row in sel.sort_values(by=["slide_count"],ascending=False).iterrows():
        if row["slope10_max"]>=min_slope_limit:
            plot_id = row["fid"]
            break # Plot found -> end selection
    # Add to list of plots to be processes
    if plot_id != None:
        plot_ids.append(plot_id)

n_sensitivity = len(plot_ids)

# Create geopackage with selected plots
plots_sensitivity_df = plots_df[plots_df["fid"].isin(plot_ids)]
plots_sensitivity_gpkg = os.path.join(input_base_path,"plots_sensitivity.gpkg")
plots_sensitivity_df.to_file(plots_sensitivity_gpkg, layer=ps_layer, driver="GPKG")

# Generate the inputs for the simulations
buffer_dist = 5
slide_buff_mask_2m_path = os.path.join(input_base_path,"slide_buff{0}_mask.tif".format(buffer_dist))
prepare_inputs(plot_ids,plots_sensitivity_gpkg,slide_buff_mask_2m_path,sensitivity=True)

#
# Sample plots without slides for optimization
#

# Determine eligible Process sources and plots
grid_df.reset_index(inplace=True)
print(len(grid_df))
grid_df = grid_df[np.invert(np.isnan(grid_df["ps_psid"]))] # Exclude plots without process source
print(len(grid_df))
grid_df = grid_df[grid_df["slide_count"]==0] # Only plots without slides
print(len(grid_df))
grid_df = grid_df[grid_df["slope10_max"]>=min_slope_limit] # Only plots with at leas one pixel slope_10m>=15°
print(len(grid_df))
grid_df = grid_df[grid_df["soil_type_majority"]!=13] # Only plots who don't mainly lie in No Slide soil type
print(len(grid_df))
grid_df = grid_df[grid_df["ps_psid"].isin(ps_df[np.invert(ps_df["has_slides"])]["psid"])] # Only plots in process sources without slides
print(len(grid_df))

# Get inputs for selection only from eligible plots
psid_eligible = grid_df["ps_psid"].unique()

ps_ids = []

X = ps_df.loc[ps_df["psid"].isin(psid_eligible)]["tecto"]
y = ps_df.loc[ps_df["psid"].isin(psid_eligible)]["fid"]
groups = ps_df.loc[ps_df["psid"].isin(psid_eligible)]["group"]

# Do selection through stratified train test split
X_train, X_test, y_train, y_test= train_test_split(X,y, train_size=int(n_sensitivity), random_state=random_seed, stratify=groups)

# Get IDs of selected plots
ps_ids = y_train.values

# Write geopackage with selected process sources
ps_optimization_df = ps_df[ps_df["fid"].isin(ps_ids)]
ps_optimization_gpkg = os.path.join(input_base_path,"ps_optimization.gpkg")
ps_optimization_df.to_file(ps_optimization_gpkg, layer=ps_layer, driver="GPKG")

# For all selected process sources, select a random plot that is above the slope_10m>15§ limit
plot_ids = []
for fid in ps_ids:
    psid = ps_df.loc[ps_df["fid"]==fid]["psid"].iloc[0]
    sel = grid_df.loc[grid_df["ps_psid"]==psid] # Only plots of the current process source
    sel = sel[sel["slope10_max"]>=min_slope_limit] # Only plots with at leas one pixel slope_10m>=15°
    plot_id = None

    row = sel.sample(n=1,random_state=(random_seed)).iloc[0] # Select random plot from candidates
    plot_id = row["fid"]

    plot_ids.append(plot_id)

# Create geopackage with selected plots
plots_optimization_df = grid_df[grid_df["fid"].isin(plot_ids)]

plots_optimization_gpkg = os.path.join(input_base_path,"plots_optimization.gpkg")
plots_optimization_df.to_file(plots_optimization_gpkg, layer=ps_layer, driver="GPKG")

soil_type_slide_2m_path = os.path.join(input_base_path,"soil_type_slide.tif")
soil_type_src = rio.open(soil_type_slide_2m_path, "r")


# Generate random "no slide" points in selected plots 
noslide_points = []
random.seed(42)
npts = 2
for index, row in plots_optimization_df.iterrows():
    [x_min,y_min,x_max, y_max] = row["geometry"].bounds

    j = 0
    while j<npts:
        x = random.uniform(x_min,x_max)
        y = random.uniform(y_min,y_max)

        
        soil_type_sample = [v for v in soil_type_src.sample([(x,y)])]

        if soil_type_sample[0] == 13 or soil_type_sample[0] == 33:
            continue 

        noslide_points.append({"x":x,"y":y})
        j+=1

# Save "no slide" points  to a geopackage
noslide_points_df = pd.DataFrame(noslide_points)
noslide_points_df = gpd.GeoDataFrame(noslide_points_df, geometry=gpd.points_from_xy(x=noslide_points_df.x, y=noslide_points_df.y), crs="EPSG:2056")
noslide_optimization_gpkg = os.path.join(input_base_path,"noslide_optimization.gpkg")
noslide_layer = "noslide_points"
noslide_points_df.to_file(noslide_optimization_gpkg, layer=noslide_layer, driver="GPKG")

# Create Mask with "no slide" points buffered by 5m at 2m resolution 
dem_2m_path = os.path.join(input_base_path,"dem_2_LV95_BE_filled.tif") 
[x_min,y_min,x_max,y_max] = utilities.get_tif_extent(dem_2m_path)

noslides_gpkg = os.path.join(input_base_path,"noslides.gpkg")
noslides_layer = "noslides"
noslides_buffered_layer = "noslides_buffered"
buffer_dist = 5
with gdal.VectorTranslate(destNameOrDestDS=noslides_gpkg, srcDS=noslide_optimization_gpkg, format="GPKG", layerName=noslides_layer, accessMode="overwrite") as ds:
    gdal.VectorTranslate(destNameOrDestDS=ds, srcDS=ds, format="GPKG", layerName=noslides_buffered_layer , SQLStatement="select st_buffer(geom,{1}) as geom from {0} ".format(noslides_layer,buffer_dist), accessMode="append")
    ds = None

noslide_buff_mask_path = os.path.join(output_base_path,"noslide_buff{0}_mask.tif".format(buffer_dist))
with gdal.Rasterize(destNameOrDestDS=noslide_buff_mask_path, srcDS=noslides_gpkg, layers=noslides_buffered_layer, format="GTiff", outputType=gdalconst.GDT_Byte, creationOptions={"COMPRESS":"DEFLATE", "BIGTIFF":"YES"}, xRes=raster_resolution, yRes=raster_resolution, outputBounds=[x_min,y_min,x_max,y_max], burnValues=1, initValues=0, noData=None, allTouched=True) as ds:
    ds=None

utilities.unset_nodata(noslide_buff_mask_path)
utilities.delete_gpkg(noslides_gpkg)

# Generate the inputs for the simulations
prepare_inputs(plot_ids,plots_optimization_gpkg,noslide_buff_mask_path,sensitivity=False)



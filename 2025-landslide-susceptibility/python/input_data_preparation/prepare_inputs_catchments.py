######################################################################
# Copyright (C) 2025 BFH
#
# Script for preparing the inputs for the SlideforMAP simulations 
# per catchment. Requires the canton wide inputs and the outputs 
# of the Prepare_Slope_Stability_Plots.model3 QGIS model (catchmants,
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
from osgeo_utils import gdal_calc
from osgeo import gdal_array
gdal.UseExceptions()


from rasterstats import zonal_stats

import fiona

import numpy as np
import geopandas as gpd
import os
import sys
import math

import tomlkit

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

hades_rasters = {}
for k in hades_files:
    hades_raster_path = os.path.join(hades_base_path,hades_files[k])
    hades_rasters[k] = hades_raster_path

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
output_base_path = "E:/GIS_Projekte/Paper_3/data/ezg_inputs"
plot_name_template = "ezg_{0}"

perimeter_buffer = 50

# Read process sources that were enriched with the region, statistics about the slides, slope at 10m, and soil type 
ps_gpkg_path = "E:/GIS_Projekte/Paper_3/data/slide_inputs/process_sources_with_slides.gpkg"
ps_layer = "plots"
ps_df = gpd.read_file(ps_gpkg_path, layer=ps_layer, engine="pyogrio", fid_as_index=True)

# IDs of the selecte catchments
catchment_ezgnr = [107220, 102117,103925,102117, 111731]
# Read catchments that were enriched with the statistics about the slides 
ezg_gpkg_path = "E:/GIS_Projekte/Paper_3/data/slide_inputs/ezg_with_slides_be.gpkg"
ezg_layer = "plots"
ezg_df = gpd.read_file(ezg_gpkg_path, layer=ezg_layer, engine="pyogrio", fid_as_index=True)

# Helper functions that clips a raster to an extent and multiplies it with a 0/1 mask
def clip_and_mask(in_raster,out_raster,mask_raster,bounds,nodata=-9999):
    tmp_raster = out_raster.replace(".tif","_tmp.tif")
    with gdal.Warp(destNameOrDestDS=tmp_raster,srcDSOrSrcDSTab=in_raster, format="GTiff", creationOptions={"COMPRESS":"DEFLATE", "BIGTIFF":"YES"}, outputBounds=bounds) as ds:
        ds=None

    with gdal_calc.Calc(calc="A*(B == 1)+(B == 0)*{0}".format(nodata),A=tmp_raster, B=mask_raster, outfile=out_raster, creation_options = ["COMPRESS=DEFLATE", "BIGTIFF=YES"] , NoDataValue=nodata, overwrite=True) as ds:
        ds=None
    
    utilities.delete_raster(tmp_raster)


utilities.ensure_dir(output_base_path)

# Process each catchment
for ezgnr in catchment_ezgnr:
    ezg_row = ezg_df.loc[ezg_df["EZGNR"]==ezgnr].iloc[0]
    ezg_teilezgnr = int(ezg_row["TEILEZGNR"])
    slide_count = int(ezg_row["slide_count"])
    slides_per_ha = float(ezg_row["slides_per_ha"])

    # Create output path for catchment
    ezg_folder = plot_name_template.format(ezgnr)
    ezg_path = os.path.join(output_base_path,ezg_folder)
    utilities.ensure_dir(ezg_path)

    # Create geopackage containing only the catchment as perimeter
    perimeter_tmp_gpkg = os.path.join(ezg_path,"perimeter_tmp.gpkg")
    perimeter_gpkg = os.path.join(ezg_path,"perimeter.gpkg")
    perimeter_layer = "perimeter"
    perimeter_buffered_layer = "perimeter_buffered"
    with gdal.VectorTranslate(destNameOrDestDS=perimeter_gpkg,  srcDS=ezg_gpkg_path, layerName=perimeter_layer, SQLStatement="SELECT * FROM {0} WHERE EZGNR={1}".format(ezg_layer,ezgnr),  accessMode="overwrite") as ds:
        ds = None

    # Get extent of unbuffered perimeter/catchment
    extent_gpkg = utilities.get_gpkg_extent(perimeter_gpkg,perimeter_layer)

    x_center = extent_gpkg[0]+(extent_gpkg[1]-extent_gpkg[0])/2
    y_center = extent_gpkg[2]+(extent_gpkg[3]-extent_gpkg[2])/2

    x_min = perimeter_buffer*(math.floor(extent_gpkg[0]/perimeter_buffer)) - perimeter_buffer
    x_max = perimeter_buffer*(math.ceil(extent_gpkg[1]/perimeter_buffer)) + perimeter_buffer
    y_min = perimeter_buffer*(math.floor(extent_gpkg[2]/perimeter_buffer)) - perimeter_buffer
    y_max = perimeter_buffer*(math.ceil(extent_gpkg[3]/perimeter_buffer)) + perimeter_buffer

    x_min -= perimeter_buffer
    y_min -= perimeter_buffer
    x_max += perimeter_buffer
    y_max += perimeter_buffer

    # Copy slides within extent to geopackage and only retain the ones within the catchment.
    with gdal.VectorTranslate(destNameOrDestDS=perimeter_gpkg, srcDS=slide_path, layerName=slide_layer, spatFilter=[x_min,y_min,x_max,y_max], accessMode="overwrite") as ds:
        # Delete slides outside the catchment
        gdal.VectorInfo(ds=ds, SQLStatement="DELETE FROM {0} WHERE fid NOT IN (SELECT s.fid FROM {0} AS s, {1} AS p WHERE ST_Intersects(s.geom,p.geom))".format(slide_layer,perimeter_layer))
        ds = None

    # Create 0/1 mask of the regular perimeter. This will be used to mask the simulation outputs and remove edge effects within the buffer.
    perimeter_mask_path = os.path.join(ezg_path,"perimeter_mask.tif")
    with gdal.Rasterize(destNameOrDestDS=perimeter_mask_path, srcDS=perimeter_gpkg, layers=perimeter_layer, format="GTiff", burnValues = 1, initValues=0, noData =0, outputType=gdalconst.GDT_Int16, creationOptions={"COMPRESS":"DEFLATE", "BIGTIFF":"YES"}, xRes=raster_resolution, yRes=raster_resolution, outputBounds=[x_min,y_min,x_max,y_max]) as ds:
        ds=None
    utilities.unset_nodata(perimeter_mask_path)
    perimeter_mask_array = gdal_array.LoadFile(perimeter_mask_path)
    perimeter_cell_count = np.sum(perimeter_mask_array[perimeter_mask_array!=-9999])

    # Create 0/1 mask of the regular perimeter at 50m resolution. This is used during result processing to check if a cell is in the perimeter
    perimeter_mask50_path = os.path.join(ezg_path,"perimeter_mask50.tif")
    with gdal.Rasterize(destNameOrDestDS=perimeter_mask50_path, srcDS=perimeter_gpkg, layers=perimeter_layer, format="GTiff", burnValues = 1, initValues=0, noData =0, outputType=gdalconst.GDT_Int16, creationOptions={"COMPRESS":"DEFLATE", "BIGTIFF":"YES"}, xRes=50, yRes=50, outputBounds=[x_min,y_min,x_max,y_max]) as ds:
        ds=None
    utilities.unset_nodata(perimeter_mask50_path)

    # Create a buffered version of the perimeter. This is used for extracting the simulation inputs.
    with gdal.VectorTranslate(destNameOrDestDS=perimeter_gpkg,  srcDS=ezg_gpkg_path, layerName=perimeter_buffered_layer, SQLStatement="SELECT ST_BUFFER(geom,{2}) geom FROM {0} WHERE EZGNR={1}".format(ezg_layer,ezgnr, perimeter_buffer),  accessMode="append") as ds:
        ds = None

    # Create 0/1 mask of the buffered perimeter. This will be used to mask the simulation inputs to save computation time.
    perimeter_mask_buffered_path = os.path.join(ezg_path,"perimeter_buffered_mask.tif")
    with gdal.Rasterize(destNameOrDestDS=perimeter_mask_buffered_path, srcDS=perimeter_gpkg, layers=perimeter_buffered_layer, format="GTiff", burnValues = 1, initValues=0, noData =0, outputType=gdalconst.GDT_Int16, creationOptions={"COMPRESS":"DEFLATE", "BIGTIFF":"YES"}, xRes=raster_resolution, yRes=raster_resolution, outputBounds=[x_min,y_min,x_max,y_max]) as ds:
        ds=None
    utilities.unset_nodata(perimeter_mask_buffered_path)

    # Clip and mask the DEM
    dem_2m_path = os.path.join(input_base_path,"dem_2_LV95_BE_filled.tif") 
    dem_path = os.path.join(ezg_path, "dem.tif")
    clip_and_mask(dem_2m_path,dem_path,perimeter_mask_buffered_path,[x_min,y_min,x_max,y_max])

    # Clip and mask the contributing area
    contributing_area_2m_name = "contributing_area_"+canton
    contributing_area_2m_path = os.path.join(input_base_path,contributing_area_2m_name+".sdat")
    contributing_area_path = os.path.join(ezg_path,"contributing_area.tif")
    with gdal.Warp(destNameOrDestDS=contributing_area_path,srcDSOrSrcDSTab=contributing_area_2m_path, format="GTiff", creationOptions={"COMPRESS":"DEFLATE", "BIGTIFF":"YES"}, outputBounds=[x_min,y_min,x_max,y_max]) as ds:
        ds=None

    # Clip and mask the failure depth
    soil_depth_2m_path = os.path.join(input_base_path,"soil_depth.tif")
    soil_depth_path = os.path.join(ezg_path,"soil_depth.tif")
    clip_and_mask(soil_depth_2m_path,soil_depth_path,perimeter_mask_buffered_path,[x_min,y_min,x_max,y_max])
   
    # Clip and mask the soil type
    soil_type_slide_2m_path = os.path.join(input_base_path,"soil_type_slide.tif")
    soil_type_slide_path = os.path.join(ezg_path,"soil_type_slide.tif")
    clip_and_mask(soil_type_slide_2m_path,soil_type_slide_path,perimeter_mask_buffered_path,[x_min,y_min,x_max,y_max])
    
    # Create a mask with the slides in the catchment buffered with 5m at 2m resolution
    slide_tmp_gpkg = os.path.join(ezg_path,"slide_tmp.gpkg")
    with gdal.VectorTranslate(destNameOrDestDS=slide_tmp_gpkg,  srcDS=perimeter_gpkg, layerName=slide_layer, SQLStatement="SELECT ST_BUFFER(geom,{1}) geom FROM {0} ".format(slide_layer, 5),  accessMode="append") as ds:
        ds = None

    slide_buff_mask_path = os.path.join(ezg_path,"slide_buff5_mask.tif")
    with gdal.Rasterize(destNameOrDestDS=slide_buff_mask_path, srcDS=slide_tmp_gpkg, layers=slide_layer, format="GTiff", burnValues = 1, initValues=0, noData = None, outputType=gdalconst.GDT_Int16, creationOptions={"COMPRESS":"DEFLATE", "BIGTIFF":"YES"}, xRes=raster_resolution, yRes=raster_resolution, outputBounds=[x_min,y_min,x_max,y_max]) as ds:
        ds=None

    utilities.delete_gpkg(slide_tmp_gpkg)

    # Clip and mask the tree type raster
    tree_type_2m_name = "tree_type_"+canton
    tree_type_2m_path = os.path.join(input_base_path,tree_type_2m_name+".tif")
    tree_type_path = os.path.join(ezg_path,"tree_type.tif")
    clip_and_mask(tree_type_2m_path,tree_type_path,perimeter_mask_buffered_path,[x_min,y_min,x_max,y_max])

    # Determine the majority of the treetype to decide the species for the simulation
    tree_type_array = gdal_array.LoadFile(tree_type_path)
    fagus_count = np.sum(tree_type_array==2)
    picea_count = np.sum(tree_type_array==1)
    tree_type = "Picea abies" if picea_count>fagus_count else "Fagus sylvatica"

    # Generate a treefile for the trees within the catchment
    treefile_gpk = os.path.join(input_base_path,"treefile.gpkg")
    treefile_gpk_tmp = os.path.join(ezg_path, "treefile_tmp.gpkg")
    treefile_buffered_gpk_tmp = os.path.join(ezg_path, "treefile_buffered_tmp.gpkg")
    treefile_path = os.path.join(ezg_path,"treefile.txt")

    # Copy trees within the extent and the perimeter to a temporary geopackage
    tree_layer = "tree_layer"
    with gdal.VectorTranslate(destNameOrDestDS=treefile_gpk_tmp, srcDS=fintch_db_path, format="GPKG", layerName=tree_layer, layers="be_processed_tree", accessMode="overwrite", spatFilter=[x_min,y_min,x_max,y_max]) as ds :
        gdal.VectorTranslate(destNameOrDestDS=ds, srcDS=perimeter_gpkg, format="GPKG", layerName=perimeter_buffered_layer, layers=perimeter_buffered_layer, accessMode="overwrite", spatFilter=[x_min,y_min,x_max,y_max])
        ds = None

    # Read only trees within the catchment
    df_out = gpd.GeoDataFrame.from_file(treefile_gpk_tmp,sql="SELECT tl.* FROM {0} AS tl, {1} AS p WHERE ST_Intersects(tl.geom,p.geom)".format(tree_layer,perimeter_buffered_layer),engine="pyogrio",)

    # Write treefile if there are any trees
    has_forest = False
    if len(df_out)>0:
        has_forest = True
        with open(treefile_path,"w+") as ofile:
            fmt = "{:>15.3f}{:>15.3f}{:>15.3f}"
            fmt = "%15.3f %15.3f %15.3f"
            np.savetxt(ofile, df_out[["x","y","bhd"]].values, fmt=fmt) #FINT-CH
            
        # Write forest mask from 10m buffered trees (not used in catchments)
        tree_buff_dist = 10
        tree_buff_layer = "trees_buffered"
        with gdal.VectorTranslate(destNameOrDestDS=treefile_buffered_gpk_tmp, srcDS=treefile_gpk_tmp, format="GPKG", layerName=tree_buff_layer, SQLStatement="select ST_buffer(geom, {0}) as geom FROM {1}".format(tree_buff_dist,tree_layer), SQLDialect="SQLite", accessMode="overwrite") as ds :
            ds = None

        forest_mask_rast = os.path.join(ezg_path, "forest_mask.tif")
        with gdal.Rasterize(destNameOrDestDS=forest_mask_rast, srcDS=treefile_buffered_gpk_tmp, layers=tree_buff_layer, format="GTiff", outputType=gdalconst.GDT_Byte, creationOptions={"COMPRESS":"DEFLATE", "BIGTIFF":"YES"}, xRes=raster_resolution, yRes=raster_resolution, outputBounds=[x_min,y_min,x_max,y_max], burnValues=1, initValues=0, noData=0) as ds:
            ds=None

    utilities.delete_gpkg(treefile_buffered_gpk_tmp)
    utilities.delete_gpkg(treefile_gpk_tmp)

    # Determine percipitation from HADES
    perc_amounts = {}

    with fiona.open(perimeter_gpkg,layer=perimeter_layer) as geom:
        for k in hades_files:
            src = hades_rasters[k]
            stats = zonal_stats(geom, src, stats="mean max")[0]
            perc = round(float(stats["mean"]/24*runoff_coefficients[k]),4)
            perc_amounts["perc_{0}".format(k)] = perc
    
    # Get slide count from process source with highest slide count in catchment
    ps_selection = ps_df[ps_df["ezg_TEILEZGNR"]==ezg_teilezgnr]
    ps_row = ps_selection.sort_values(by=["slide_count"],ascending=False).iloc[0]
    ps_id = int(ps_row["psid"])
    ps_slide_count = int(ps_row["slide_count"])
    ps_area = ps_row.geometry.area
    ps_slides_per_ha = ps_slide_count/(ps_area/10000)

    ps_layer = "plots"
    ps_out_layer = "process_sources"
    with gdal.VectorTranslate(destNameOrDestDS=perimeter_gpkg,  srcDS=ps_gpkg_path, layerName=ps_out_layer, SQLStatement="SELECT * FROM {0} WHERE ezg_TEILEZGNR={1}".format(ps_layer,ezg_teilezgnr),  accessMode="overwrite") as ds:
        ds = None


    # Write config file by modifying the "default_config.toml" template
    toml_template_path = os.path.join(input_base_path,"default_config.toml")
    toml_output_path = os.path.join(ezg_path,"default_config.toml")
    with open(toml_template_path, mode="rt") as fp:
        config = tomlkit.load(fp)

        config["perimeter_cell_count"] = int(perimeter_cell_count)
        config["x_min"] = x_min
        config["y_min"] = y_min
        config["x_max"] = x_max
        config["y_max"] = y_max

        for k in perc_amounts:
            config[k] = perc_amounts[k]

        config["catchment_slide_count"] = slide_count
        config["catchment_slide_per_ha"] = slides_per_ha
        config["ps_slide_count"] = ps_slide_count
        config["ps_slide_per_ha"] = ps_slides_per_ha        
        config["ps_psid"] = ps_id

        config["demPath"] = "../{0}".format(os.path.basename(dem_path))
        config["flowAccPath"] = "../{0}".format(os.path.basename(contributing_area_path))
        config["soilTypeRaster"] =  "../{0}".format(os.path.basename(soil_type_slide_path))    
        config["soilThicknessValueRaster"] = "../{0}".format(os.path.basename(soil_depth_path)) 

        if has_forest:
            config["treeFilePath"] = "../{0}".format(os.path.basename(treefile_path))   
            config["treeSpecies"] = tree_type

        with open(toml_output_path, mode="wt") as fp:
            tomlkit.dump(config, fp)


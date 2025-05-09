######################################################################
# Copyright (C) 2025 BFH
#
# Script for preparing all necessary inputs for the SlideforMAP
# simulations across the entire canton of Bern.
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
gdal.UseExceptions()

import numpy as np
import geopandas as gpd
import os
import sys 

import math

# Adding the path to the utilities.py to the path and import
utility_path = "E:/GIS_Projekte/Paper_3/code/utilities"
sys.path.append(utility_path)
import utilities

#
# Configs and paths
#

saga_path = "C:/OSGeo4W64/apps/saga-ltr"

dhm_rast = "E:/GIS_Projekte/Geodaten/swissalti3D_BE/dem_2_LV95.tif"
dhm_eu = "E:/GIS_Projekte/Geodaten/EU-DEM_v1_1/eu_dem_v11_E40N20.TIF"

mixing_degree_rast = "E:/GIS_Projekte/FINT-CH/MG2020_0306_oneRF_S1S2_AS_DEM_LV95.tif"
fintch_db_path = "F:/fint-ch/AP10__Validierung_der_EBD_Analysen/FINT_v2_BE.gpkg"
swissBoundaries_gpkg = "E:/GIS_Projekte/Geodaten/swissboundaries3d_2024-01/swissboundaries3d_2024-01_2056_5728.gpkg/swissBOUNDARIES3D_1_5_LV95_LN02.gpkg"

tlm_gpkg = "E:/GIS_Projekte/Geodaten/swisstlm3d_2024-03_2056_5728/SWISSTLM3D_2024_LV95_LN02.gpkg"

slide_path = "E:/GIS_Projekte/Paper_3/data/BE_SL_all_2025_T01/BE_SL_all_2025_T01.shp"
slides_layer = "slides"

geocover_gpkg = "E:/GIS_Projekte/Geodaten/Geocover/geocover_merged.gpkg"
deposit_layer_name = "Unconsolidated_Deposits_PLG"
bedrock_layer_name = "Bedrock_PLG"

soil_type_mapping_path = "E:/GIS_Projekte/GHK/BE/Soil_Type_Mappings_BE.csv"
soil_type_rockfall_path = "E:/GIS_Projekte/GHK/BE/Adelboden/Ground_Properties/soiltype.tif"

fintch_layer_name = "processed_trees"
tbk_layer_name = "tbk_bestandeskarte"

raster_resolution = 2
raster_resolution_l2 = 5
raster_resolution_l3 = 10

slide_buffer = 1200

output_base_path = "E:/GIS_Projekte/Paper_3/data/slide_inputs_2"

soil_depth_raster = "E:/GIS_Projekte/Paper_2/data/geomorph/ch_tiles/BE_Perimeter/pred_rf_model_hmdb.tif"

utilities.ensure_dir(output_base_path)

#
# Buffering the perimeter
#     

border_full_file_name = os.path.split(swissBoundaries_gpkg)[1]
border_file_name = border_full_file_name.replace(".gpkg","")

canton_table = "tlm_kantonsgebiet"

canton_number = 2
canton = "BE"
perimeter_buffer = 0

perimeter_file_name = "perimeter_"+canton
perimeter_path = os.path.join(output_base_path,perimeter_file_name+".gpkg")

with gdal.VectorTranslate(destNameOrDestDS=perimeter_path, srcDS=swissBoundaries_gpkg, format="GPKG", layerName=perimeter_file_name, SQLStatement="select ST_buffer(geom, {0}) as geom FROM {1} WHERE Kantonsnummer={2}".format(perimeter_buffer,canton_table,canton_number), SQLDialect="SQLite", accessMode="overwrite") as ds :
    ds = None

# Get the extent of the buffered perimeter
extent_gpkg = utilities.get_gpkg_extent(perimeter_path,perimeter_file_name)
print(extent_gpkg)
tile_buffer = 1000

# Round extent to coarsest resolution and add a buffer (overlap when using tiles)
x_min = raster_resolution_l3*(math.floor(extent_gpkg[0]/raster_resolution_l3)) - raster_resolution_l3 - tile_buffer
x_max = raster_resolution_l3*(math.ceil(extent_gpkg[1]/raster_resolution_l3)) + raster_resolution_l3 + tile_buffer
y_min = raster_resolution_l3*(math.floor(extent_gpkg[2]/raster_resolution_l3)) - raster_resolution_l3 - tile_buffer
y_max = raster_resolution_l3*(math.ceil(extent_gpkg[3]/raster_resolution_l3)) + raster_resolution_l3 + tile_buffer

print([x_min,x_max,y_min, y_max])

#
# Clip DHM 
#

dhm_rast_full_file_name = os.path.split(dhm_rast)[1]
dhm_file_name = dhm_rast_full_file_name.replace(".tif","")

dhm_clip_name = dhm_file_name+"_raw_"+canton
dhm_clip_raw_path = os.path.join(output_base_path,dhm_clip_name+".tif")

# Clip DHM to extent
with gdal.Warp(destNameOrDestDS=dhm_clip_raw_path,srcDSOrSrcDSTab=dhm_rast, format="GTiff", creationOptions={"COMPRESS":"DEFLATE", "BIGTIFF":"YES"}, outputBounds=[x_min,y_min,x_max,y_max]) as ds:
    ds=None

# Clip soil depth to extent
soil_depth_clip_path = os.path.join(output_base_path,"soil_depth.tif")
with gdal.Warp(destNameOrDestDS=soil_depth_clip_path,srcDSOrSrcDSTab=soil_depth_raster, xRes=raster_resolution, yRes=raster_resolution,  resampleAlg=gdal.GRA_NearestNeighbour, format="GTiff", creationOptions={"COMPRESS":"DEFLATE", "BIGTIFF":"YES"}, outputBounds=[x_min,y_min,x_max,y_max]) as ds:
    ds=None


#
# Clip EU DHM
#

# swissTLMRegio might be a better alternative with 10m resolution

[x_min,y_min,x_max,y_max] = utilities.get_tif_extent(dhm_clip_raw_path)

dhm_rast_eu_full_file_name = os.path.split(dhm_eu)[1]
dhm_eu_file_name = dhm_rast_eu_full_file_name.lower().replace(".tif","")

dhm_eu_clip_name = dhm_eu_file_name+"_"+canton
dhm_eu_clip_path = os.path.join(output_base_path,dhm_eu_clip_name+".tif")
with gdal.Warp(destNameOrDestDS=dhm_eu_clip_path,srcDSOrSrcDSTab=dhm_eu, xRes=raster_resolution, yRes=raster_resolution, dstSRS="EPSG:2056", resampleAlg=gdal.GRA_Bilinear, format="GTiff", creationOptions={"COMPRESS":"DEFLATE", "BIGTIFF":"YES", "TILED":"YES"}, outputBounds =[x_min,y_min,x_max,y_max]) as ds:
    ds=None

#
# Fill CH DEM voids
#

# Not strictly necessary for Bern but border cantons would need this

utilities.unset_nodata(dhm_clip_raw_path)
utilities.unset_nodata(dhm_eu_clip_path)

dhm_clip_name = dhm_file_name+"_"+canton
dhm_clip_path = os.path.join(output_base_path,dhm_clip_name+".tif")

# Keep A where >-1 (i.e. not NoData), else take B
with gdal_calc.Calc(calc="A*(A>-1)+B*(A<-1)",A=dhm_clip_raw_path, B=dhm_eu_clip_path, outfile=dhm_clip_path, creation_options = ["COMPRESS=DEFLATE", "BIGTIFF=YES"] , overwrite=True) as ds:
    ds=None


#
# Calculate contributing area
#

# trenslate to SAGA grid
dhm_clip_sgrd_path = dhm_clip_path.replace(".tif",".sdat") 
with gdal.Translate(destName=dhm_clip_sgrd_path, srcDS=dhm_clip_path, format="SAGA") as ds:
    ds=None

# Fill sinks/depressions using Wang & Liu XXL (slightly simplified but a lot faster than regular Wang & Liu on big rasters)
dhm_clip_sgrd_filled_path = dhm_clip_sgrd_path.replace(".sdat","_filled.sdat") 
fill_sinks_cmd = "{0}/saga_cmd ta_preprocessor 5 -ELEV {1} -FILLED {2}" #Wang & Liu XXL
dhm_fill_cmd = fill_sinks_cmd.format(saga_path, dhm_clip_sgrd_path.replace(".sdat",".sgrd"), dhm_clip_sgrd_filled_path.replace(".sdat",".sgrd"))
os.system(dhm_fill_cmd)

# Calculate contributing area usint Deterministic 8 algorithm
contributung_area_name = "contributing_area_"+canton
contributung_area_path = os.path.join(output_base_path,contributung_area_name+".sdat")

catchment_area_cmd = "{0}/saga_cmd ta_hydrology 0 -ELEVATION {1} -METHOD 0 -FLOW {2}" #Deterministic 8
dhm_ca_cmd = catchment_area_cmd.format(saga_path, dhm_clip_sgrd_filled_path.replace(".sdat",".sgrd"),contributung_area_path)
os.system(dhm_ca_cmd)

#Translate outputs to GeoTIFF
dhm_clip_filled_path = dhm_clip_sgrd_filled_path.replace(".sdat",".tif")
with gdal.Translate(destName=dhm_clip_filled_path, srcDS=dhm_clip_sgrd_filled_path, format="GTiff", creationOptions={"COMPRESS":"DEFLATE", "BIGTIFF":"YES"}) as ds:
    ds=None


# Resample DHM to 5m

[x_min,y_min,x_max,y_max] = utilities.get_tif_extent(dhm_clip_filled_path)

# Mainly necessary as derivative for later runout simulations
dhm_clip_sgrd_filled_5m_path = dhm_clip_sgrd_filled_path.replace(".sdat","_5m.tif")
with gdal.Warp(destNameOrDestDS=dhm_clip_sgrd_filled_5m_path,srcDSOrSrcDSTab=dhm_clip_filled_path, xRes=raster_resolution_l2,yRes=raster_resolution_l2, dstSRS="EPSG:2056", outputBounds=[x_min,y_min,x_max,y_max], resampleAlg=gdal.GRA_Average,  format="GTiff", creationOptions={"COMPRESS":"DEFLATE", "BIGTIFF":"YES", "TILED":"YES"}) as ds:
    ds=None

# Mainly necessary to calculate Slope at 10m
dhm_clip_sgrd_filled_10m_path = dhm_clip_sgrd_filled_path.replace(".sdat","_10m.tif")
with gdal.Warp(destNameOrDestDS=dhm_clip_sgrd_filled_10m_path,srcDSOrSrcDSTab=dhm_clip_filled_path, xRes=raster_resolution_l3,yRes=raster_resolution_l3, dstSRS="EPSG:2056", outputBounds=[x_min,y_min,x_max,y_max], resampleAlg=gdal.GRA_Average,  format="GTiff", creationOptions={"COMPRESS":"DEFLATE", "BIGTIFF":"YES", "TILED":"YES"}) as ds:
    ds=None

# Calculate slope for the different resolutions
dhm_slope_name = dhm_file_name+"_slope_"+canton
dhm_slope_path = os.path.join(output_base_path,dhm_slope_name+".tif")
with gdal.DEMProcessing(destName=dhm_slope_path, srcDS=dhm_clip_filled_path, processing="slope", slopeFormat="degree", creationOptions={"COMPRESS":"DEFLATE", "BIGTIFF":"YES"}) as ds:
    ds=None        

dhm_slope_5m_name = dhm_file_name+"_5m_slope_"+canton
dhm_slope_5m_path = os.path.join(output_base_path,dhm_slope_5m_name+".tif")
with gdal.DEMProcessing(destName=dhm_slope_5m_path, srcDS=dhm_clip_sgrd_filled_5m_path, processing="slope", slopeFormat="degree", creationOptions={"COMPRESS":"DEFLATE", "BIGTIFF":"YES"}) as ds:
    ds=None        

dhm_slope_10m_name = dhm_file_name+"_10m_slope_"+canton
dhm_slope_10m_path = os.path.join(output_base_path,dhm_slope_10m_name+".tif")
with gdal.DEMProcessing(destName=dhm_slope_10m_path, srcDS=dhm_clip_sgrd_filled_10m_path, processing="slope", slopeFormat="degree", creationOptions={"COMPRESS":"DEFLATE", "BIGTIFF":"YES"}) as ds:
    ds=None   

# Slope mask >50° for 5m resolution
slope_mask_path = os.path.join(output_base_path,"slope_gt50_mask.tif")
with gdal_calc.Calc(calc="(A>50)*1",A=dhm_slope_path, outfile=slope_mask_path, creation_options = ["COMPRESS=DEFLATE", "BIGTIFF=YES"]) as ds:
    ds=None


#
# Generate binary raster with slide locations
#

# Rasterizing "as is" at 2m resolution
slide_mask_path = os.path.join(output_base_path,"slide_mask.tif")
with gdal.Rasterize(destNameOrDestDS=slide_mask_path, srcDS=slide_path, format="GTiff", outputType=gdalconst.GDT_Byte, creationOptions={"COMPRESS":"DEFLATE", "BIGTIFF":"YES"}, xRes=raster_resolution, yRes=raster_resolution, outputBounds=[x_min,y_min,x_max,y_max], burnValues=1, initValues=0) as ds:
    ds=None

# Rasterizing slides buffered by 5m at 2m resolution
slides_gpkg = os.path.join(output_base_path,"slides.gpkg")
slides_layer = "slides"
slides_buffered_layer = "slides_buffered"
buffer_dist = 9.8
with gdal.VectorTranslate(destNameOrDestDS=slides_gpkg, srcDS=slide_path, format="GPKG", layerName=slides_layer, accessMode="overwrite") as ds:
    gdal.VectorTranslate(destNameOrDestDS=ds, srcDS=ds, format="GPKG", layerName=slides_buffered_layer , SQLStatement="select st_buffer(geom,{1}) as geom from {0} ".format(slides_layer,buffer_dist), accessMode="append")
    ds = None

slide_buff_mask_path = os.path.join(output_base_path,"slide_buff{0}_mask.tif".format(buffer_dist))
with gdal.Rasterize(destNameOrDestDS=slide_buff_mask_path, srcDS=slides_gpkg, layers=slides_buffered_layer, format="GTiff", outputType=gdalconst.GDT_Byte, creationOptions={"COMPRESS":"DEFLATE", "BIGTIFF":"YES"}, xRes=raster_resolution, yRes=raster_resolution, outputBounds=[x_min,y_min,x_max,y_max], burnValues=1, initValues=0, noData=None, allTouched=True) as ds:
    ds=None

utilities.unset_nodata(slide_buff_mask_path)

#
# TLM Forest mask
#

# Extract forest polygons from ground cover
forest_gpkg = os.path.join(output_base_path,"tlm_forest.gpkg")
tlm_groundcover_layer = "tlm_bb_bodenbedeckung"
forest_layer = "forest"
with gdal.VectorTranslate(destNameOrDestDS=forest_gpkg, srcDS=tlm_gpkg, format="GPKG", layerName=forest_layer , SQLStatement="select * from {0} WHERE (OBJEKTART = 'Wald' OR  OBJEKTART = 'Wald offen')".format(tlm_groundcover_layer), spatFilter=[x_min,y_min,x_max,y_max]) as ds :
    ds = None

# Rasterize at 2m resolution
forest_mask_path = os.path.join(output_base_path,"forest_mask.tif")
with gdal.Rasterize(destNameOrDestDS=forest_mask_path, srcDS=forest_gpkg, layers=forest_layer, format="GTiff", outputType=gdalconst.GDT_Byte, creationOptions={"COMPRESS":"DEFLATE", "BIGTIFF":"YES"}, xRes=raster_resolution, yRes=raster_resolution, outputBounds=[x_min,y_min,x_max,y_max], burnValues=1, initValues=0) as ds:
    ds=None

#
# TLM Fels Lockergestein mask
#

# Extract rock and loose rock polygons from ground cover
temp_gpkg = os.path.join(output_base_path,"tlm_tmp.gpkg")
tlm_groundcover_layer = "tlm_bb_bodenbedeckung"
rock_deposit_layer = "rock_deposit"
with gdal.VectorTranslate(destNameOrDestDS=temp_gpkg, srcDS=tlm_gpkg, format="GPKG", layerName=rock_deposit_layer , SQLStatement="select * from {0} WHERE (OBJEKTART = 'Fels' OR  OBJEKTART = 'Lockergestein')".format(tlm_groundcover_layer), spatFilter=[x_min,y_min,x_max,y_max]) as ds :
    ds = None

# Rasterize at 2m resolution
rock_deposit_path = os.path.join(output_base_path,"{0}.tif".format(rock_deposit_layer))
with gdal.Rasterize(destNameOrDestDS=rock_deposit_path, srcDS=temp_gpkg, layers=rock_deposit_layer, format="GTiff", outputType=gdalconst.GDT_Byte, creationOptions={"COMPRESS":"DEFLATE", "BIGTIFF":"YES"}, xRes=raster_resolution, yRes=raster_resolution, outputBounds=[x_min,y_min,x_max,y_max], burnValues=1, initValues=0) as ds:
    ds=None


#
# Generate uniform soil type raster for sensitivity analysis
#

soil_type_value = 1

soil_type_name = "soil_type_"+canton
soil_type_path = os.path.join(output_base_path,soil_type_name+".tif")

with gdal_calc.Calc(calc="(A>=0)*{0}".format(soil_type_value),A=dhm_clip_path, type=gdalconst.GDT_Byte, outfile=soil_type_path, creation_options = ["COMPRESS=DEFLATE", "BIGTIFF=YES"], overwrite=True) as ds:
    ds=None


#
# Generate Treefile
#
input_treefile = os.path.join(output_base_path,"treefile.txt")
treefile_gpk = os.path.join(output_base_path,"treefile.gpkg")

# Clip FINT-CH result to extent
treefile_gpk_tmp = os.path.join(output_base_path, "treefile_tmp.gpkg")
tree_layer = "tree_layer"
with gdal.VectorTranslate(destNameOrDestDS=treefile_gpk_tmp, srcDS=fintch_db_path, format="GPKG", layerName=tree_layer, layers="be_processed_tree", accessMode="overwrite", spatFilter=[x_min,y_min,x_max,y_max]) as ds : #FINT-CH
    ds = None

# Read input data
df_out = gpd.GeoDataFrame.from_file(treefile_gpk_tmp,layer=tree_layer)

# Write cantonwide file
with open(input_treefile,"w+") as ofile:
    fmt = "{:>15.3f}{:>15.3f}{:>15.3f}"
    fmt = "%15.3f %15.3f %15.3f"
    np.savetxt(ofile, df_out[["x","y","bhd"]].values, fmt=fmt) #FINT-CH

utilities.delete_gpkg(treefile_gpk_tmp)


#
# Caclulate Soil type from soil property data
#

kobo_folder_path = "E:/GIS_Projekte/Paper_2/data/020_HAFL_Phd_Projekt_Rutschungen/Daten"

sand_0_name = "Soil_sand_depth_0_30.tif"
sand_30_name = "Soil_sand_depth_30_60.tif"
sand_60_name = "Soil_sand_depth_60_120.tif"
silt_0_name = "Soil_silt_depth_0_30.tif"
silt_30_name = "Soil_silt_depth_30_60.tif"
silt_60_name = "Soil_silt_depth_60_120.tif"
clay_0_name = "Soil_clay_depth_0_30.tif"
clay_30_name = "Soil_clay_depth_30_60.tif"
clay_60_name = "Soil_clay_depth_60_120.tif"

# Convert soil proprty raters to 2m resolution
kobo_files = [sand_0_name,sand_30_name,sand_60_name,silt_0_name,silt_30_name,silt_60_name,clay_0_name,clay_30_name,clay_60_name]
kobo_paths = []
for f in kobo_files:
    in_path = os.path.join(kobo_folder_path,f)
    out_path = os.path.join(output_base_path,f)
    kobo_paths.append(out_path)
    with gdal.Warp(destNameOrDestDS=out_path,srcDSOrSrcDSTab=in_path, format="GTiff", xRes=raster_resolution, yRes=raster_resolution, resampleAlg=gdal.GRA_NearestNeighbour, creationOptions={"COMPRESS":"DEFLATE", "BIGTIFF":"YES", "TILED":"YES"} , outputBounds=[x_min,y_min,x_max,y_max]) as ds:
        ds=None

# Function generating a raster that selects the value from th soil property rasters basd on a depth raster
def select_values_by_depth(soil_depth,values_0,values_30,values_60,out_raster):
    with gdal_calc.Calc(calc="(A <= 0.30) * B + (A > 0.30) * (A <= 0.60) * C + (A > 0.60) * D",A=soil_depth, B=values_0, C=values_30, D=values_60, outfile=out_raster, type=gdalconst.GDT_Float32, creation_options = ["COMPRESS=DEFLATE", "BIGTIFF=YES"], overwrite=True) as ds:
        ds=None

# Select sand, silt, and clay contents based on predicted failure depth
sand_raster = os.path.join(output_base_path,"sand_selected.tif")
silt_raster = os.path.join(output_base_path,"silt_selected.tif")
clay_raster = os.path.join(output_base_path,"clay_selected.tif")
select_values_by_depth(soil_depth_clip_path, kobo_paths[0], kobo_paths[1], kobo_paths[2], sand_raster)
select_values_by_depth(soil_depth_clip_path, kobo_paths[3], kobo_paths[4], kobo_paths[5], silt_raster)
select_values_by_depth(soil_depth_clip_path, kobo_paths[6], kobo_paths[7], kobo_paths[8], clay_raster)

# Function for deriving USDA soil texture types from sand, silt, and clay contents
def generate_usda_soil_types(sand_raster, silt_raster, clay_raster, out_raster, nodata_value = 13):
    # Logic from soil texture calculator
    # ((silt + 1.5*clay) < 15)	SAND
    # ((silt + 1.5*clay >= 15) && (silt + 2*clay < 30))	LOAMY SAND
    # ((clay >= 7 && clay < 20) && (sand > 52) && ((silt + 2*clay) >= 30) || (clay < 7 && silt < 50 && (silt+2*clay)>=30))	SANDY LOAM
    # ((clay >= 7 && clay < 27) && (silt >= 28 && silt < 50) && (sand <= 52))	LOAM
    # ((silt >= 50 && (clay >= 12 && clay < 27)) || ((silt >= 50 && silt < 80) && clay < 12))	SILT LOAM
    # (silt >= 80 && clay < 12)	SILT
    # ((clay >= 20 && clay < 35) && (silt < 28) && (sand > 45)) 	SANDY CLAY LOAM
    # ((clay >= 27 && clay < 40) && (sand > 20 && sand <= 45))	CLAY LOAM
    # ((clay >= 27 && clay < 40) && (sand  <= 20))	SILTY CLAY LOAM
    # (clay >= 35 && sand > 45)	SANDY CLAY
    # (clay >= 40 && silt >= 40)	SILTY CLAY
    # (clay >= 40 && sand <= 45 && silt < 40)	CLAY

    out_s1_path = out_raster.replace(".tif","_s1.tif")
    calc_s1 = "maximum(  ((B + 1.5*C) < 15) * 1,maximum(  (((B + 1.5*C) >= 15) * ((B + 2*C) < 30)) * 2,maximum(  ((((C >= 7) * (C < 20)) * (A > 52) * ((B + 2*C) >= 30)) + ((C < 7) * (B < 50) * ((B+2*C)>=30))) * 3,maximum(  (((C >= 7) * (C < 27)) * ((B >= 28) * (B < 50)) * (A <= 52)) * 4,maximum(  (((B >= 50) * ((C >= 12) * (C < 27))) + (((B >= 50) * (B < 80)) * (C < 12))) * 5,maximum(  ((B >= 80) * (C < 12)) * 6,maximum(  (((C >= 20) * (C < 35)) * (B < 28) * (A > 45)) * 7,maximum(  (((C >= 27) * (C < 40)) * ((A > 20) * (A <= 45))) * 8,maximum(  (((C >= 27) * (C < 40)) * (A  <= 20)) * 9,maximum(  ((C >= 35) * (A > 45)) * 10,maximum(  ((C >= 40) * (B >= 40)) * 11,maximum(  ((C >= 40) * (A <= 45) * (B < 40)) * 12, 0))))))))))))"
    with gdal_calc.Calc(calc=calc_s1,A=sand_raster, B=silt_raster, C=clay_raster, outfile=out_s1_path, type=gdalconst.GDT_Byte, creation_options = ["COMPRESS=DEFLATE", "BIGTIFF=YES"], NoDataValue = nodata_value, overwrite=True) as ds:
        ds=None

    out_s2_path = out_raster.replace(".tif","_s2.tif")
    calc_s2 = "(D > 0) * D + (D == 0) * maximum(  (((100-A-C) + 1.5*C) < 15) * 1,maximum(  ((((100-A-C) + 1.5*C) >= 15) * (((100-A-C) + 2*C) < 30)) * 2,maximum(  ((((C >= 7) * (C < 20)) * (A > 52) * (((100-A-C) + 2*C) >= 30)) + ((C < 7) * ((100-A-C) < 50) * (((100-A-C)+2*C)>=30))) * 3,maximum(  (((C >= 7) * (C < 27)) * (((100-A-C) >= 28) * ((100-A-C) < 50)) * (A <= 52)) * 4,maximum(  ((((100-A-C) >= 50) * ((C >= 12) * (C < 27))) + ((((100-A-C) >= 50) * ((100-A-C) < 80)) * (C < 12))) * 5,maximum(  (((100-A-C) >= 80) * (C < 12)) * 6,maximum(  (((C >= 20) * (C < 35)) * ((100-A-C) < 28) * (A > 45)) * 7,maximum(  (((C >= 27) * (C < 40)) * ((A > 20) * (A <= 45))) * 8,maximum(  (((C >= 27) * (C < 40)) * (A  <= 20)) * 9,maximum(  ((C >= 35) * (A > 45)) * 10,maximum(  ((C >= 40) * ((100-A-C) >= 40)) * 11,maximum(  ((C >= 40) * (A <= 45) * ((100-A-C) < 40)) * 12, 0))))))))))))"
    with gdal_calc.Calc(calc=calc_s2,A=sand_raster, B=silt_raster, C=clay_raster, D=out_s1_path, outfile=out_raster, type=gdalconst.GDT_Byte, creation_options = ["COMPRESS=DEFLATE", "BIGTIFF=YES"], NoDataValue = nodata_value, overwrite=True) as ds:
        ds=None
    utilities.unset_nodata(out_s1_path)
    utilities.unset_nodata(out_raster)


soil_type_kobo_depth_path = os.path.join(output_base_path, "soil_type_usda.tif") 
generate_usda_soil_types(sand_raster,silt_raster,clay_raster,soil_type_kobo_depth_path)

# Convert to SAGA gris
soil_type_kobo_depth_sgrd_path = soil_type_kobo_depth_path.replace(".tif",".sdat") 
with gdal.Translate(destName=soil_type_kobo_depth_sgrd_path, srcDS=soil_type_kobo_depth_path, format="SAGA") as ds:
    ds=None

# Smoothen with Mode filter of 15m radius
majority_filter_cmd = "{0}/saga_cmd grid_filter  6 -INPUT {1} -RESULT {2} -MODE 1 -RADIUS {3}" 
soil_type_kobo_depth_sgrd_filtered_path = soil_type_kobo_depth_sgrd_path.replace(".sdat","_majority.sdat") 

soiltype_majority_filter_cmd = majority_filter_cmd.format(saga_path, soil_type_kobo_depth_sgrd_path.replace(".sdat",".sgrd"), soil_type_kobo_depth_sgrd_filtered_path.replace(".sdat",".sgrd"),15)
os.system(soiltype_majority_filter_cmd)

# Translate SAGA grid to GeoTIFF
soil_type_kobo_depth_filtered_path = soil_type_kobo_depth_sgrd_filtered_path.replace(".sdat",".tif") 
with gdal.Translate(destName=soil_type_kobo_depth_filtered_path, srcDS=soil_type_kobo_depth_sgrd_filtered_path, format="GTiff") as ds:
    ds=None

# Combine with forest mask to derive soilt types in forest
soil_type_kobo_depth_filtered_forest_path = soil_type_kobo_depth_sgrd_filtered_path.replace(".sdat","_forest.tif") 
with gdal_calc.Calc(calc="A + 20 * B",A=soil_type_kobo_depth_filtered_path, B=forest_mask_path, outfile=soil_type_kobo_depth_filtered_forest_path, type=gdalconst.GDT_Byte, creation_options = ["COMPRESS=DEFLATE", "BIGTIFF=YES"], overwrite=True) as ds:
    ds=None

# Mask No Slide Areas
soil_type_slide_path = os.path.join(output_base_path, "soil_type_slide.tif")
calc_mask = "((B + C) > 0)*13+((B + C) <= 0) * A"
with gdal_calc.Calc(calc=calc_mask,A=soil_type_kobo_depth_filtered_forest_path, B=slope_mask_path, C=rock_deposit_path, outfile=soil_type_slide_path, type=gdalconst.GDT_Byte, creation_options = ["COMPRESS=DEFLATE", "BIGTIFF=YES"], NoDataValue = 13, overwrite=True) as ds:
    ds=None

utilities.unset_nodata(soil_type_slide_path)

#
# Generate tree type raster
#

mg_clip_name = "mg_clip_"+canton
mg_clip_path = os.path.join(output_base_path,mg_clip_name+".tif")

# Convert to 2m resolution
with gdal.Warp(destNameOrDestDS=mg_clip_path,srcDSOrSrcDSTab=mixing_degree_rast, xRes=raster_resolution, yRes=raster_resolution,  format="GTiff", creationOptions={"COMPRESS":"DEFLATE", "BIGTIFF":"YES", "TILED":"YES"}, outputBounds=[x_min,y_min,x_max,y_max]) as ds:
    ds=None
utilities.unset_nodata(mg_clip_path)

# Generate tree type raster. Coniferous proportion >=50 -> 1 else 2
tree_type_name = "tree_type_"+canton
tree_type_path = os.path.join(output_base_path,tree_type_name+".tif")
with gdal_calc.Calc(calc="2*(A>-1)*(A<5000)+1*(A>-1)*(A>=5000)",A=mg_clip_path, B=forest_mask_path, outfile=tree_type_path, type=gdalconst.GDT_Byte, creation_options = ["COMPRESS=DEFLATE", "BIGTIFF=YES"]) as ds:
    ds=None

#
# Cleanup of unneeded rasters
#
utilities.delete_raster(dhm_eu_clip_path)
utilities.delete_raster(dhm_clip_raw_path)
utilities.delete_raster(dhm_clip_sgrd_path)
utilities.delete_raster(dhm_clip_sgrd_filled_path)

utilities.delete_raster(rock_deposit_path)
utilities.delete_raster(slope_mask_path)

utilities.delete_raster(sand_raster)
utilities.delete_raster(silt_raster)
utilities.delete_raster(clay_raster)

utilities.delete_raster(soil_type_kobo_depth_sgrd_path)
utilities.delete_raster(soil_type_kobo_depth_sgrd_filtered_path)

for kobo_path in kobo_paths:
    utilities.delete_raster(kobo_path)

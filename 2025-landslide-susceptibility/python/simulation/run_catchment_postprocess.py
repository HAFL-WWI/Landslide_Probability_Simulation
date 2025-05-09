######################################################################
# Copyright (C) 2025 BFH
#
# Script for postprocessing the simulation results in the catchments 
# for the different return periods into absolute landslide probabilities 
# based on different criteria.
#
# Author: Christoph Schaller, BFH-HAFL, March 2025
######################################################################

#
# Imports
#

import os
import math
import glob

from osgeo import gdal
from osgeo import ogr
from osgeo.gdalconst import *
from osgeo.gdalnumeric import *
from osgeo_utils import gdal_calc
from osgeo import gdal_array
gdal.UseExceptions()

import csv

import tomlkit

import rasterio as rio
import numpy as np
import pandas as pd
import geopandas as gpd

import rasterstats

from sklearn.metrics import roc_auc_score

# Function to ensure that a directory exists
def ensure_dir (path):
    if not os.path.isdir(path):
        try:
            os.mkdir(path)
            return path
        except FileExistsError:
            return -1

# Function that returns the extent of a GeoTIFF raster
def get_tif_extent (raster_path):
    inDataSource = gdal.Open(raster_path, gdal.GA_ReadOnly)
    geoTransform = inDataSource.GetGeoTransform()
    minx = geoTransform[0]
    maxy = geoTransform[3]
    maxx = minx + geoTransform[1] * inDataSource.RasterXSize
    miny = maxy + geoTransform[5] * inDataSource.RasterYSize
    extent = [minx, miny, maxx, maxy]
    inDataSource = None
    return extent

# Function that deletes an existing Geotiff
def delete_raster(raster):
    data = gdal.Open(raster, gdal.GA_ReadOnly)
    driver = data.GetDriver()
    data = None
    if os.path.exists(raster):
        driver.Delete(raster)   

# Functions for deleting vector files
def delete_vector_file(path,driver_name):
    DriverName = driver_name
    driver = ogr.GetDriverByName(DriverName)
    if os.path.exists(path):
        driver.DeleteDataSource(path)

# Function to delete a shapefile inxluding all partial files
def delete_shapefile(path):
    delete_vector_file(path = path, driver_name = "ESRI Shapefile")

# Function to delete a geopackage
def delete_gpkg (path):
    delete_vector_file(path = path, driver_name = "GPKG")

#
# Business logic
#

# Setup the records describing the work units for processing
def process_setup(cfg, num_processes = 1):
    
    # Read the parameter values
    result_base_path = cfg["result_base_path"]

    # Loop all catchment folders
    records = []
    for plot_path in next(os.walk(result_base_path))[1]:
        plot_name = os.path.basename(plot_path)
        plot_number = int(plot_name.replace("ezg_",""))

        # Create record for each parameter set that was simulated
        for soil_properties_path in glob.glob(os.path.join(result_base_path,"SoilPhysicalProperties_*.csv")):
            soil_properties_name = os.path.basename(soil_properties_path)
            output_name = soil_properties_name.replace(".csv","").replace("SoilPhysicalProperties_","")

            rec = {
                "soil_properties_path": soil_properties_path,
                "output_name":output_name,
                "plot_path":os.path.join(result_base_path,plot_name),
                "plot_id": plot_number,
            }
            records.append(rec)

    return records


# Variant of the function for merging propabilities from multiple rainfall scenarios into one absolute landslide probability. 
# The correction factor is calculated using the sunction specified by the variant.
def merge_probabilities(output_folder,sim_folder,return_periods,observation_period,variant,perimeter_path,ezg_layer,ps_layer,slide_layer):
    """
    This function processes susceptibility TIFF files for different return periods to obtain an absolute release probability map. 
    (1) Each return period susceptibility is first weighted with the occurence probability that an event for that return period will occur during the period T (observation_perion).
    (2) Susceptibilities are weighted using the occurrence probability of the return period and then summed
    (3) Result from (2) is then scaled using the correction factor calculated by the specified variants function based on the ratio of observed landslides divided by the summed probability from (2) in a reference area
    (4) Result from (3) is then aggregated to 50m and 100m and sume summary statistic are calculated and returned 

    Author: Denis Cohen-Corticchiato, Christoph Schaller
    """
    # Define the base directory containing the adelboden* directories
    base_dir = output_folder

    # Directories to merge
    rp_dirs = [os.path.join(base_dir,sim_folder,"rp_{0}".format(rp)) for rp in return_periods]

    # Initialize variables to hold the sum of weighted TIFF files
    weighted_sum_occurence = None   # RLOP
    weighted_sum_release = None     # AbsLSprob

    # Loop through each directory
    for i, dir_path in enumerate(rp_dirs):
        weight = 1.0 - (1.0 - 1.0/return_periods[i])**(observation_period)
        if i < len(return_periods) - 1:
            weight_release = 1/return_periods[i] - 1/return_periods[i+1]
        else:
            weight_release = 1/return_periods[i]

        # Construct the path to the failure probability tif file
        tiff_path = os.path.join(dir_path, 'probability_failure.tif')
        
        # Open the TIFF file
        with rio.open(tiff_path) as src:
            # Read the data
            tiff_data = src.read(1)  # Read the first band
            
            # If weighted_sum_occurence is not initialized, initialize it with zeros
            if weighted_sum_occurence is None:
                weighted_sum_occurence = np.zeros_like(tiff_data, dtype=np.float64)
            # Add the weighted data to the weighted_sum_occurence
            weighted_sum_occurence += weight * tiff_data

            # If weighted_sum_release is not initialized, initialize it with zeros
            if weighted_sum_release is None:
                weighted_sum_release = np.zeros_like(tiff_data, dtype=np.float64)
            # Add the weighted data to the weighted_sum_occurence
            weighted_sum_release += weight_release * tiff_data

    # Save the occurence probability to TIFF files
    epsg = 2056
    occ_prob_file_name = "RLOP_sim_tot"
    occ_prob_file_path = os.path.join(base_dir, sim_folder, occ_prob_file_name + ".tif")
    with rio.open(
        occ_prob_file_path,
        'w',
        driver='GTiff',
        height=weighted_sum_occurence.shape[0],
        width=weighted_sum_occurence.shape[1],
        count=1,
        dtype=weighted_sum_occurence.dtype,
        crs=src.crs,
        transform=src.transform,
        compress="DEFLATE",
    ) as dst:
        dst.write(weighted_sum_occurence, 1)

    # Calculate the correction factor using the respective function in the variant record
    corr_factor = variant["corr_function"](output_folder,sim_folder,perimeter_path,ezg_layer,ps_layer,slide_layer,occ_prob_file_path)
    print(corr_factor)
    weighted_sum_release *= corr_factor

    with rio.open(occ_prob_file_path) as src:
        # Read the data
        tiff_data = src.read(1)  # Read the first band

    # Save the release probability to TIFF files
    rel_prob_file_name = "AbsLSprob_{0}.tif".format(variant["name"])
    rel_prob_file_path = os.path.join(base_dir, sim_folder, rel_prob_file_name)
    epsg = 2056
    with rio.open(
        rel_prob_file_path,
        'w',
        driver='GTiff',
        height=weighted_sum_release.shape[0],
        width=weighted_sum_release.shape[1],
        count=1,
        dtype=weighted_sum_release.dtype,
        crs=src.crs,
        transform=src.transform,
        compress="DEFLATE",
    ) as dst:
        dst.write(weighted_sum_release, 1)

    # Aggregate Result to 50m
    input_raster_path = rel_prob_file_path 
    resolution_50 = 50
    output_raster_50_path = os.path.join(base_dir, sim_folder,rel_prob_file_name.replace(".tif","_50m.tif"))
    with gdal.Warp(destNameOrDestDS=output_raster_50_path,srcDSOrSrcDSTab=input_raster_path, dstSRS="EPSG:{0}".format(epsg), srcSRS="EPSG:{0}".format(epsg), xRes=resolution_50, yRes=resolution_50,  resampleAlg=gdalconst.GRA_Sum, format="GTiff", creationOptions={"COMPRESS=DEFLATE", "BIGTIFF=YES", "TILED=YES"}) as ds:
        ds = None
    
    # Aggregate Result to 100m
    resolution_100 = 100

    # Get bound rounded to 100m
    [x_min,y_min,x_max,y_max] = get_tif_extent(input_raster_path)
    x_min = resolution_100*(math.floor(x_min/resolution_100))
    x_max = resolution_100*(math.ceil(x_max/resolution_100)) + resolution_100
    y_min = resolution_100*(math.floor(y_min/resolution_100))
    y_max = resolution_100*(math.ceil(y_max/resolution_100)) + resolution_100

    output_raster_100_path = os.path.join(base_dir, sim_folder,rel_prob_file_name.replace(".tif","_100m.tif"))
    with gdal.Warp(destNameOrDestDS=output_raster_100_path,srcDSOrSrcDSTab=input_raster_path, xRes=resolution_100, yRes=resolution_100,  resampleAlg=gdalconst.GRA_Sum, format="GTiff", outputBounds=[x_min,y_min,x_max,y_max], creationOptions={"COMPRESS=DEFLATE", "BIGTIFF=YES", "TILED=YES"}) as ds:
        ds = None

    # Calculate summary stats for this variant
    variant_stat = {}
    variant_stat["name"] = variant["name"]
    variant_stat["sim"] = sim_folder
    variant_stat["ezg"] = os.path.basename(os.path.dirname(perimeter_path))
    variant_stat["correction_factor"] = corr_factor

    prob_array = gdal_array.LoadFile(rel_prob_file_path)
    prob_50_array = gdal_array.LoadFile(output_raster_50_path)
    prob_100_array = gdal_array.LoadFile(output_raster_100_path)
    variant_stat["sum"] = np.nansum(prob_array)
    variant_stat["max_50"] = np.nanmax(prob_50_array)
    variant_stat["max_100"] = np.nanmax(prob_100_array)

    return variant_stat

# Function for calculating a range of statistics for analysis from a combined probability raster 
def generate_results(output_folder,sim_folder, perimeter_path,ps_layer, slide_layer, variant):
    # Paths and configs
    plot_path = output_folder
    base_path = plot_path
    sim_path = output_path


    rel_prob_filename = "AbsLSprob_{0}.tif".format(variant["name"]) 
    AbsLSprob_path = os.path.join(sim_path,rel_prob_filename)
    AbsLSprob50_path = os.path.join(sim_path, rel_prob_filename.replace(".tif","_50m.tif"))
    AbsLSprob100_path = os.path.join(sim_path, rel_prob_filename.replace(".tif","_100m.tif"))
    perimeter_mask50_path = os.path.join(base_path,"perimeter_mask50.tif")
    perimeter_mask_path = os.path.join(base_path,"perimeter_mask.tif")

    # Determine stats per process source
    ps_df = gpd.read_file(perimeter_path,layer=ps_layer, engine="pyogrio", fid_as_index=True)
    ps_df["psid"] = ps_df["psid"].astype("int32")

    stats = rasterstats.zonal_stats(perimeter_path, AbsLSprob_path, layer=ps_layer, stats="count sum min max mean median percentile_95", geojson_out=True)
    stat_recs = []
    for geom in stats:
        rec = {"psid":int(geom["properties"]["psid"]),
        "count":geom["properties"]["count"],
        "sum":geom["properties"]["sum"],
        "min":geom["properties"]["min"],
        "max":geom["properties"]["max"],
        "mean":geom["properties"]["mean"],
        "median":geom["properties"]["median"],
        "percentile_95":geom["properties"]["percentile_95"]
        }
        stat_recs.append(rec)
    stat_df = pd.DataFrame(stat_recs)
    ps_df = ps_df.merge(stat_df,right_on="psid",left_on="psid",how="left",suffixes=["","_stat"])

    ps_stats_gpkg = os.path.join(sim_path,"ps_stats_{0}.gpkg".format(variant["name"]))
    delete_gpkg(ps_stats_gpkg)
    ps_df.to_file(ps_stats_gpkg, layer=ps_layer, driver="GPKG")


    # Determine AUC overall (i.e. for all observed landslides at once)
    slide_buff_mask_path = os.path.join(base_path,"slide_buff5_mask.tif")

    slide_actual = gdal_array.LoadFile(slide_buff_mask_path).flatten()
    slide_prediction = gdal_array.LoadFile(AbsLSprob_path).flatten()
    slide_prediction[np.isnan(slide_prediction)] = 0

    auc = roc_auc_score(slide_actual,slide_prediction)
    res = 1-auc

    comb_path = os.path.join(sim_path,"auc_all_{0}.csv".format(variant["name"]))
    if os.path.exists(comb_path):
        os.remove(comb_path)

    with open(comb_path, "w") as file:
        writer = csv.writer(file, delimiter=",")

        writer.writerow(["sim","auc"] )    
        writer.writerow([output_name, auc])  

    # Determine probability per slide
    slide_df = gpd.read_file(perimeter_path,layer=slide_layer, engine="pyogrio", fid_as_index=True)
    slide_df.reset_index(inplace=True)
    coord_list = [(x, y) for x, y in zip(slide_df["geometry"].x, slide_df["geometry"].y)]

    with rio.open(AbsLSprob_path, "r") as src:
        slide_df["prob"] = [x for x in src.sample(coord_list)]
        # slide_df["prob"] = slide_df["prob"].apply(lambda x:  x[0] if isinstance(x, list) else x)
        slide_df["prob"] = slide_df["prob"].apply(lambda x:  x[0])

    with rio.open(AbsLSprob50_path, "r") as src:
        slide_df["prob50"] = [x for x in src.sample(coord_list)]
        # slide_df["prob50"] = slide_df["prob50"].apply(lambda x: x[0] if isinstance(x, list) else x)
        slide_df["prob50"] = slide_df["prob50"].apply(lambda x: x[0])

    with rio.open(AbsLSprob100_path, "r") as src:
        slide_df["prob100"] = [x for x in src.sample(coord_list)]
        # slide_df["prob100"] = slide_df["prob100"].apply(lambda x: x[0] if isinstance(x, list) else x)
        slide_df["prob100"] = slide_df["prob100"].apply(lambda x: x[0])

    # Determine AUC per slide
    raster_resolution = cfg["raster_resolution"]
    
    silde_stat_path = os.path.join(sim_path,"auc_slides_{0}.csv".format(variant["name"]))
    if os.path.exists(silde_stat_path):
        os.remove(silde_stat_path)

    with open(silde_stat_path, "w") as file:
        writer = csv.writer(file, delimiter=",")

        writer.writerow(["sim","fid","rsid","prob","prob50","prob100","min","max","mean","p5","p10","q1","median","q2","p90","p95","auc"] )    

        # For each slide
        for index, row in slide_df.iterrows():
            x = row["geometry"].x
            y = row["geometry"].y

            # Create geopackage with just one slide buffered by 5m
            tmp_slide_gpkg = os.path.join(sim_path,"tmp_slide.gpkg")
            tmp_slide_layer = "slides"
            tmp_mask_path = os.path.join(sim_path,"tmp_slide.tif")
            tmp_prob_path = os.path.join(sim_path,"tmp_prob.tif")

            with gdal.VectorTranslate(destNameOrDestDS=tmp_slide_gpkg,  srcDS=perimeter_path, layerName=tmp_slide_layer, SQLStatement="SELECT ST_BUFFER(geom,{2}) geom FROM {0} WHERE fid={1}".format(slide_layer,row["fid"], 5),  accessMode="append") as ds:
                ds = None

            # Get extent for rasterization
            plot_buffer = 10
            x_min = x-plot_buffer
            x_max = x+plot_buffer
            y_min = y-plot_buffer
            y_max = y+plot_buffer

            x_min = raster_resolution*(math.floor(x_min/raster_resolution))
            x_max = raster_resolution*(math.ceil(x_max/raster_resolution)) 
            y_min = raster_resolution*(math.floor(y_min/raster_resolution)) 
            y_max = raster_resolution*(math.ceil(y_max/raster_resolution))

            # Rasterize buffered slide
            with gdal.Rasterize(destNameOrDestDS=tmp_mask_path, srcDS=tmp_slide_gpkg, layers=tmp_slide_layer, format="GTiff", burnValues = 1, initValues=0, noData = None, outputType=gdalconst.GDT_Int16, creationOptions={"COMPRESS":"DEFLATE", "BIGTIFF":"YES"}, xRes=raster_resolution, yRes=raster_resolution, outputBounds=[x_min,y_min,x_max,y_max]) as ds:
                ds=None

            # Clip probability raster to extent 
            with gdal.Warp(destNameOrDestDS=tmp_prob_path,srcDSOrSrcDSTab=AbsLSprob_path, format="GTiff", creationOptions={"COMPRESS":"DEFLATE", "BIGTIFF":"YES"}, outputBounds=[x_min,y_min,x_max,y_max]) as ds:
                ds=None

            # Calculate AUC
            slide_actual = gdal_array.LoadFile(tmp_mask_path).flatten()
            slide_prediction = gdal_array.LoadFile(tmp_prob_path).flatten()

            auc = roc_auc_score(slide_actual,slide_prediction)

            # Calculate stats within slide area
            prob_selection = slide_prediction[slide_actual==1]
            pf_min = np.min(prob_selection) 
            pf_max = np.max(prob_selection) 
            pf_mean = np.mean(prob_selection) 
            [pf_5,p_10,pf_q1, pf_median, pf_q2,pf_90,pf_95] = np.nanpercentile(prob_selection,[5,10,25,50,75,90,95]) 
            writer.writerow([output_name, row["fid"], row["RS_ID"],row["prob"],row["prob50"],row["prob100"],pf_min, pf_max, pf_mean, pf_5,p_10,pf_q1, pf_median, pf_q2,pf_90,pf_95, auc])  

            # Remove temp files
            delete_raster(tmp_mask_path)
            delete_raster(tmp_prob_path)
            delete_gpkg(tmp_slide_gpkg)    


    # Determine stats per slide in 50m cells
    [x_min,y_min,x_max,y_max] = get_tif_extent(AbsLSprob50_path)

    resolution = 50
    
    slide_layer = "slides"
    slide_df = gpd.read_file(perimeter_path,layer=slide_layer, engine="pyogrio", fid_as_index=True)
    slide_df.reset_index(inplace=True)
    slide_df["x"] = slide_df["geometry"].x
    slide_df["y"] = slide_df["geometry"].y

    slide_df["x_corner"] = slide_df["x"].apply(lambda x: int(math.floor(x/resolution)*resolution))
    slide_df["y_corner"] = slide_df["y"].apply(lambda x: int(math.floor(x/resolution)*resolution))
    
    prob_array = gdal_array.LoadFile(AbsLSprob50_path)
    mask_array = gdal_array.LoadFile(perimeter_mask50_path)

    rast_records = []
    slide_records = []

    y_curr = y_max
    # For each 50m cell
    while y_curr>y_min:
        x_curr = x_min
        while x_curr<x_max:
            # Get indices from coordinates
            row = int((y_max-y_curr)/resolution)
            col = int((x_curr-x_min)/resolution)

            x_corner = x_curr
            y_corner = y_curr-resolution

            # If cell is in catchment
            if mask_array[row][col] == 1:
                # Add record to list with raster cell probabilities
                rast_records.append({"x":x_corner,"y":y_corner,"prob50":prob_array[row][col]})

                # Determine stats for slides within the cell
                slides_cell = slide_df.loc[(slide_df["x_corner"]==x_corner) & (slide_df["y_corner"]==y_corner)]
                prob50 =  prob_array[row][col]
                prob_slide = 0
                n = 0
                year_min = 0
                year_max = 0
                year_distinct = 0
                if len(slides_cell)>0:
                    n = len(slides_cell)
                    year_max = int(slides_cell["year"].max())
                    year_min = int(slides_cell["year"].min())
                    year_distinct = int(len(slides_cell["year"].unique()))
                    prob_slide = year_distinct/observation_period
                slide_records.append({"x":x_corner,"y":y_corner,"prob50":prob50,"prob_slide":prob_slide,"n":n,"year_min":year_min,"year_max":year_max,"year_distinct":year_distinct})

            x_curr+=resolution
        y_curr-=resolution

    # Remove CSV if exists 
    slide_stat_path = os.path.join(sim_path,"prob50_slide_{0}.csv".format(variant["name"]))
    if os.path.exists(slide_stat_path):
        os.remove(slide_stat_path)

    # Save slide stats to CSV
    with open(slide_stat_path, "w") as file:
        writer = csv.writer(file, delimiter=",")

        writer.writerow(["x","y","prob50","prob_slide","n","year_min","year_max","year_distinct"] )    

        for rec in slide_records:
            writer.writerow([rec["x"],rec["y"],rec["prob50"],rec["prob_slide"],rec["n"],rec["year_min"],rec["year_max"],rec["year_distinct"]])

    # Determine stats per slide in 100m cells
    [x_min,y_min,x_max,y_max] = get_tif_extent(AbsLSprob100_path)

    resolution = 100
    
    slide_layer = "slides"
    slide_df = gpd.read_file(perimeter_path,layer=slide_layer, engine="pyogrio", fid_as_index=True)
    slide_df.reset_index(inplace=True)
    slide_df["x"] = slide_df["geometry"].x
    slide_df["y"] = slide_df["geometry"].y

    slide_df["x_corner"] = slide_df["x"].apply(lambda x: int(math.floor(x/resolution)*resolution))
    slide_df["y_corner"] = slide_df["y"].apply(lambda x: int(math.floor(x/resolution)*resolution))
    
    perimeter_mask100_path = perimeter_mask_path.replace(".tif","100.tif")
    with gdal.Warp(destNameOrDestDS=perimeter_mask100_path,srcDSOrSrcDSTab=perimeter_mask_path, xRes=resolution, yRes=resolution,  resampleAlg=gdalconst.GRA_Max, format="GTiff", outputBounds=[x_min,y_min,x_max,y_max], creationOptions={"COMPRESS=DEFLATE", "BIGTIFF=YES", "TILED=YES"}) as ds:
        ds = None

    prob_array = gdal_array.LoadFile(AbsLSprob100_path)
    mask_array = gdal_array.LoadFile(perimeter_mask100_path)

    rast_records = []
    slide_records = []

    y_curr = y_max
    # For each 100m cell
    while y_curr>y_min:
        x_curr = x_min
        while x_curr<x_max:
            # Get indices from coordinates
            row = int((y_max-y_curr)/resolution)
            col = int((x_curr-x_min)/resolution)

            x_corner = x_curr
            y_corner = y_curr-resolution

            # print(col,row,x_corner,y_corner)
            if mask_array[row][col] == 1:
                # Add record to list with raster cell probabilities
                rast_records.append({"x":x_corner,"y":y_corner,"prob50":prob_array[row][col]})

                # Determine stats for slides within the cell
                slides_cell = slide_df.loc[(slide_df["x_corner"]==x_corner) & (slide_df["y_corner"]==y_corner)]
                prob100 =  prob_array[row][col]
                prob_slide = 0
                n = 0
                year_min = 0
                year_max = 0
                year_distinct = 0
                if len(slides_cell)>0:
                    n = len(slides_cell)
                    year_max = int(slides_cell["year"].max())
                    year_min = int(slides_cell["year"].min())
                    year_distinct = int(len(slides_cell["year"].unique()))
                    prob_slide = year_distinct/observation_period
                slide_records.append({"x":x_corner,"y":y_corner,"prob100":prob100,"prob_slide":prob_slide,"n":n,"year_min":year_min,"year_max":year_max,"year_distinct":year_distinct})

            x_curr+=resolution
        y_curr-=resolution

    # Remove CSV if exists 
    slide_stat_path = os.path.join(sim_path,"prob100_slide_{0}.csv".format(variant["name"]))
    if os.path.exists(slide_stat_path):
        os.remove(slide_stat_path)

    # Save slide stats to CSV
    with open(slide_stat_path, "w") as file:
        writer = csv.writer(file, delimiter=",")

        writer.writerow(["x","y","prob100","prob_slide","n","year_min","year_max","year_distinct"] )    

        for rec in slide_records:
            writer.writerow([rec["x"],rec["y"],rec["prob100"],rec["prob_slide"],rec["n"],rec["year_min"],rec["year_max"],rec["year_distinct"]])

# Function that calculates the correction factor using the 40 km2 catchment as reference
def get_correction_factor_catchment(output_folder,sim_folder,perimeter_path,ezg_layer,ps_layer,slide_layer,prob_path):
    # Get number of observed slides in catchment
    ezg_df = gpd.read_file(perimeter_path,layer=ezg_layer, engine="pyogrio", fid_as_index=True)
    slide_count = ezg_df.iloc[0]["slide_count"]

    # Calculate sum of probability in catchment (prob raster is masked to catchment)
    prob_array = gdal_array.LoadFile(prob_path)
    cur_rlop_sum = np.nansum(prob_array)
    corr_factor = slide_count / cur_rlop_sum
    
    return corr_factor    

# Function that calculates the correction factor using the process source with the highest number of observed landslides with the highest summed probability
def get_correction_factor_ps_fixed(output_folder,sim_folder,perimeter_path,ezg_layer,ps_layer,slide_layer,prob_path):
    # Calculate sum of the probability in the process sources using zonal statistics
    stats = rasterstats.zonal_stats(perimeter_path, prob_path, layer=ps_layer, stats="count sum", geojson_out=True)

    cur_rlop_sum = 0
    cur_event_count = 0
    # Loop all process sources
    for geom in stats:
        if geom["properties"]["slide_count"] == None: # No observed slides -> ignore and continue
            continue
        if geom["properties"]["slide_count"]>cur_event_count: # Process source has more obsrve slides thatn previously encoutered ones -> copy values
            print("a",geom["properties"]["psid"],geom["properties"]["slide_count"],geom["properties"]["sum"],)
            cur_event_count = geom["properties"]["slide_count"]
            cur_rlop_sum = geom["properties"]["sum"]
        elif geom["properties"]["slide_count"]==cur_event_count and geom["properties"]["sum"]>cur_rlop_sum: # Process source has same count as current max but higher probability -> copy values
            print("b",geom["properties"]["psid"],geom["properties"]["slide_count"],geom["properties"]["sum"],)
            cur_event_count = geom["properties"]["slide_count"]
            cur_rlop_sum = geom["properties"]["sum"]

    if cur_rlop_sum>0:
        corr_factor = cur_event_count / cur_rlop_sum
    else: # Use default correction factor is most active source has no probability
        corr_factor = 0.001
        print("Warning: default correction factor",output_folder,sim_folder,corr_factor)

    return corr_factor

# Function that calculates the correction factor using the process source with the highest number of observed landslides that has sliding probability
def get_correction_factor_ps(output_folder,sim_folder,perimeter_path,ezg_layer,ps_layer,slide_layer,prob_path):
    # Calculate sum of the probability in the process sources using zonal statistics
    stats = rasterstats.zonal_stats(perimeter_path, prob_path, layer=ps_layer, stats="count sum", geojson_out=True)

    cur_rlop_sum = 0
    cur_event_count = 0
    # Loop all process sources
    for geom in stats:
        if geom["properties"]["slide_count"] == None: # No observed slides -> ignore and continue
            continue
        if geom["properties"]["slide_count"]>cur_event_count and geom["properties"]["sum"]>0: # Process source has a higher observed landslide count and probability -> copy values 
            print("a",geom["properties"]["psid"],geom["properties"]["slide_count"],geom["properties"]["sum"],)
            cur_event_count = geom["properties"]["slide_count"]
            cur_rlop_sum = geom["properties"]["sum"]
        elif geom["properties"]["slide_count"]==cur_event_count and geom["properties"]["sum"]>cur_rlop_sum: # Process source has same count as current max but higher probability -> copy values
            print("b",geom["properties"]["psid"],geom["properties"]["slide_count"],geom["properties"]["sum"],)
            cur_event_count = geom["properties"]["slide_count"]
            cur_rlop_sum = geom["properties"]["sum"]

    if cur_rlop_sum>0:
        corr_factor = cur_event_count / cur_rlop_sum
    else: # Use default correction factor is most active source has no probability (unlikely to happen here)
        corr_factor = 0.001
        print("Warning: default correction factor",output_folder,sim_folder,corr_factor)

    return corr_factor

# Function that calculates the correction factor using the hectare with the highest number of observed landslides with the highest summed probability
def get_correction_factor_ha_fixed(output_folder,sim_folder,perimeter_path,ezg_layer,ps_layer,slide_layer,prob_path):
    # Determine rounded extent for hectare grid
    resolution = 100

    [x_min,y_min,x_max,y_max] = get_tif_extent(prob_path)
    x_min = resolution*(math.floor(x_min/resolution))
    x_max = resolution*(math.ceil(x_max/resolution)) + resolution
    y_min = resolution*(math.floor(y_min/resolution))
    y_max = resolution*(math.ceil(y_max/resolution)) + resolution

    # Create aggrecated probability raster matching the grid
    prob_100_path = prob_path.replace(".tif","_100.tif")
    with gdal.Warp(destNameOrDestDS=prob_100_path,srcDSOrSrcDSTab=prob_path, xRes=resolution, yRes=resolution,  resampleAlg=gdalconst.GRA_Sum, format="GTiff", outputBounds=[x_min,y_min,x_max,y_max], creationOptions={"COMPRESS=DEFLATE", "BIGTIFF=YES", "TILED=YES"}) as ds:
        ds = None

    # Read probability raster
    prob_array = gdal_array.LoadFile(prob_100_path)

    # Read observed slides
    slide_df = gpd.read_file(perimeter_path,layer=slide_layer, engine="pyogrio", fid_as_index=True)
    slide_df.reset_index(inplace=True)

    # "Grid" slides to lower left corner of containing hectare cell
    slide_df["x"] = slide_df["geometry"].x
    slide_df["y"] = slide_df["geometry"].y
    slide_df["x_corner"] = slide_df["x"].apply(lambda x: int(math.floor(x/resolution)*resolution))
    slide_df["y_corner"] = slide_df["y"].apply(lambda x: int(math.floor(x/resolution)*resolution))

    cur_rlop_sum = 0
    cur_event_count = 0

    # Loop all cells
    y_curr = y_max
    while y_curr>y_min:
        x_curr = x_min
        while x_curr<x_max:
            # Get array indices from coordinates
            row = int((y_max-y_curr)/resolution)
            col = int((x_curr-x_min)/resolution)

            # Get slides matching the coorner coordinates of the current cell
            x_corner = x_curr
            y_corner = y_curr-resolution
            slides_cell = slide_df.loc[(slide_df["x_corner"]==x_corner) & (slide_df["y_corner"]==y_corner)]
            slide_count = len(slides_cell)

            # Get probability in cell
            prob_sum = prob_array[row][col]

            if slide_count>cur_event_count: # Hectare has more observed slides -> copy values
                print("a",slide_count,prob_sum)
                cur_event_count = slide_count
                cur_rlop_sum = prob_sum
            elif slide_count==cur_event_count and prob_sum>cur_rlop_sum: # Hectare has same observed slides as current max but higher probability -> copy values
                print("b",slide_count,prob_sum)
                cur_event_count = slide_count
                cur_rlop_sum = prob_sum

            x_curr+=resolution
        y_curr-=resolution

    if cur_rlop_sum>0:
        corr_factor = cur_event_count / cur_rlop_sum
    else:
        corr_factor = 0.001 # Use default correction factor is most active source has no probability
        print("Warning: default correction factor",output_folder,sim_folder,corr_factor)

    # Remove aggregated raster
    delete_raster(prob_100_path)

    return corr_factor

# Function that calculates the correction factor using the hectare with the highest number of observed landslides that has sliding probability
def get_correction_factor_ha(output_folder,sim_folder,perimeter_path,ezg_layer,ps_layer,slide_layer,prob_path):
    # Determine rounded extent for hectare grid
    resolution = 100

    [x_min,y_min,x_max,y_max] = get_tif_extent(prob_path)
    x_min = resolution*(math.floor(x_min/resolution))
    x_max = resolution*(math.ceil(x_max/resolution)) + resolution
    y_min = resolution*(math.floor(y_min/resolution))
    y_max = resolution*(math.ceil(y_max/resolution)) + resolution

    # Create aggrecated probability raster matching the grid
    prob_100_path = prob_path.replace(".tif","_100.tif")
    with gdal.Warp(destNameOrDestDS=prob_100_path,srcDSOrSrcDSTab=prob_path, xRes=resolution, yRes=resolution,  resampleAlg=gdalconst.GRA_Sum, format="GTiff", outputBounds=[x_min,y_min,x_max,y_max], creationOptions={"COMPRESS=DEFLATE", "BIGTIFF=YES", "TILED=YES"}) as ds:
        ds = None

    # Read probability raster
    prob_array = gdal_array.LoadFile(prob_100_path)

    # Read observed slides
    slide_df = gpd.read_file(perimeter_path,layer=slide_layer, engine="pyogrio", fid_as_index=True)
    slide_df.reset_index(inplace=True)

    # "Grid" slides to lower left corner of containing hectare cell
    slide_df["x"] = slide_df["geometry"].x
    slide_df["y"] = slide_df["geometry"].y
    slide_df["x_corner"] = slide_df["x"].apply(lambda x: int(math.floor(x/resolution)*resolution))
    slide_df["y_corner"] = slide_df["y"].apply(lambda x: int(math.floor(x/resolution)*resolution))

    cur_rlop_sum = 0
    cur_event_count = 0

    # Loop all cells
    y_curr = y_max
    while y_curr>y_min:
        x_curr = x_min
        while x_curr<x_max:
            # Get array indices from coordinates
            row = int((y_max-y_curr)/resolution)
            col = int((x_curr-x_min)/resolution)

            # Get slides matching the coorner coordinates of the current cell
            x_corner = x_curr
            y_corner = y_curr-resolution
            slides_cell = slide_df.loc[(slide_df["x_corner"]==x_corner) & (slide_df["y_corner"]==y_corner)]
            slide_count = len(slides_cell)

            # Get probability in cell
            prob_sum = prob_array[row][col]

            if slide_count==cur_event_count and prob_sum>cur_rlop_sum: # Hectare has same observed slides as current max but higher probability -> copy values
                print("a",slide_count,prob_sum)
                cur_event_count = slide_count
                cur_rlop_sum = prob_sum
            elif slide_count>cur_event_count and prob_sum>0: # Hectare has more observed slides as current max and has probability -> copy values
                print("b",slide_count,prob_sum)
                cur_event_count = slide_count
                cur_rlop_sum = prob_sum

            x_curr+=resolution
        y_curr-=resolution

    corr_factor = cur_event_count / cur_rlop_sum

    # Remove aggregated raster
    delete_raster(prob_100_path)

    return corr_factor

# Default entry point
if __name__ == "__main__":

    # Main configuration for processing
    cfg = {
        "result_base_path": "/mnt/data/ezg_inputs", # path where the simulations will be written
        # "result_base_path": "E:/GIS_Projekte/Paper_3/data/ezg_inputs/ezg_inputs", # path where the simulations will be written
        "sfm_path": "/home/administrator/slideformap/slideformap_v2.2.2/SlideforMAP", # path to SlideforMAP executable
        "output_name_template": "opt_{0}", # prefix for simulation output folders
        "num_processes": 1, # number of parallel simulations
        "num_cores": 2, # number of cores per instance
        "execute_simulations": True, # Whether to actually run the simulation or just generate the setup
        "remove_outputs": False, # whether to delete "unneeded" simulation outputs  
        "observation_period": 62, # observation period for actual slides
        "slide_buffer_size": 5, # Buffer to use on historical slides
        "slide_threshold": 0.0001, # observation period for actual slides
        "raster_resolution": 2, # resolution of the probability raster
    }

    # Variants for calculating the correction factor to test
    variants = [
        {"name":"catchment","corr_function":get_correction_factor_catchment},
        # {"name":"process_source_fixed","corr_function":get_correction_factor_ps_fixed},
        {"name":"process_source","corr_function":get_correction_factor_ps},
        # {"name":"hectare_fixed","corr_function":get_correction_factor_ha_fixed},
        {"name":"hectare","corr_function":get_correction_factor_ha},
    ]

    result_base_path = cfg["result_base_path"]

    # Get the records describing the catchments and simulations to calculate the results for
    sim_records = process_setup(cfg, num_processes = 1)
    print(len(sim_records))

    variant_stats = []

    # Execute for each catchment and simulation
    for sim_record in sim_records:
        # Get parameters for simulation from config
        output_name = sim_record["output_name"]  
        
        plot_path = sim_record["plot_path"]
        soil_properties_path = sim_record["soil_properties_path"]
        output_path = os.path.join(plot_path,output_name)

        perimeter_mask_path = os.path.join(plot_path,"perimeter_mask.tif") 

        observation_period = cfg["observation_period"]

        # Create output directory
        ensure_dir(output_path)

        # Merge scenario probabilities
        
        base_path = plot_path
        sim_path = output_path

        slide_buff_mask_path = os.path.join(base_path,"slide_buff5_mask.tif")

        # Get information on catcment stored in the config template 
        toml_template_path = os.path.join(plot_path,"default_config.toml")       
        catchment_slide_count = 0
        catchment_slide_per_ha = 0.00001
        ps_id = -1
        ps_slide_count = 0 
        ps_slide_per_ha = -1
        perimeter_cell_count = 0

        with open(toml_template_path, mode="rt") as fp:
            config = tomlkit.load(fp)

            perimeter_cell_count = config["perimeter_cell_count"] 
            catchment_slide_count = config["catchment_slide_count"]
            catchment_slide_per_ha = config["catchment_slide_per_ha"]
            ps_slide_count = config["ps_slide_count"] 
            ps_slide_per_ha = config["ps_slide_per_ha"]      
            ps_id = config["ps_psid"]


        # Calculate Zonal statistics for occurence probability on process sources        
        return_periods = [2,10,30,100,300]
        print(sim_path)

        # Calculate for each variant
        for variant in variants: 
            print(variant)

            # Remove combined probability raster if exists
            AbsLSprob_path = os.path.join(sim_path,"AbsLSprob_{0}.tif".format(variant["name"]))
            if os.path.exists(AbsLSprob_path):
                os.remove(AbsLSprob_path) 

            # Merge the scenarios to an absolute probability
            perimeter_gpkg = os.path.join(plot_path,"perimeter.gpkg")
            ps_layer = "process_sources"
            slide_layer = "slides"
            ezg_layer = "perimeter"
            res = merge_probabilities(plot_path,output_name,return_periods,observation_period,variant,perimeter_gpkg,ezg_layer,ps_layer,slide_layer)
            
            # Generate additional summary statistics
            generate_results(plot_path,output_name,perimeter_gpkg,ps_layer,slide_layer,variant)
            variant_stats.append(res)
    
    # Merge summary statistics from variants to a single CSV
    variant_stat_csv = os.path.join(result_base_path,"variant_stats.csv")
    variant_stats_df = pd.DataFrame(variant_stats)
    variant_stats_df.to_csv(variant_stat_csv,sep=",",index=False)
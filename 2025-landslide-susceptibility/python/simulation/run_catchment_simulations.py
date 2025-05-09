######################################################################
# Copyright (C) 2025 BFH
#
# Script for running simulations for different return periods
# in the catchments.
#
# Author: Christoph Schaller, BFH-HAFL, March 2025
######################################################################

#
# Imports
#

import os
import math
import glob

import re

from osgeo import gdal
from osgeo import ogr
from osgeo.gdalconst import *
from osgeo.gdalnumeric import *
from osgeo_utils import gdal_calc
from osgeo import gdal_array
gdal.UseExceptions()

import csv

from datetime import datetime, date, time, timedelta
import time

from multiprocessing import Process, Pool, Queue, JoinableQueue, Manager, current_process, freeze_support
from queue import Empty

import tomlkit

import rasterio as rio
import numpy as np
import geopandas as gpd

from sklearn.metrics import roc_auc_score

#
# Helper functions
#

# Function to ensure that a directory exists
def ensure_dir (path):
    if not os.path.isdir(path):

        try:
            os.mkdir(path)
            return path
        except FileExistsError:
            return -1

# Function that deletes an existing raster using gdal
def delete_raster(raster):
    data = gdal.Open(raster, gdal.GA_ReadOnly)
    driver = data.GetDriver()
    data = None
    if os.path.exists(raster):
        driver.Delete(raster)   

# Function that deletes an existing 
def delete_vector_file(path,driver_name):
    DriverName = driver_name
    driver = ogr.GetDriverByName(DriverName)
    if os.path.exists(path):
        driver.DeleteDataSource(path)

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
# Parallel processing functions
#

# Worker function that will handle the processing of the individual records
def worker(q, work_function, cfg):
    
    while True:
        #Consume work as long as there are items in the queue
        try:
            job_record = q.get()
            if job_record == None:
                q.task_done()
                print("Queue End")
                break

            res = work_function(job_record, cfg)

            q.task_done()
        except Empty:
            # print("Queue empty")
            break
    #No more work available
    print("Exit:",current_process())
    return

# Process records in parallel fashion 
def process_records(process_records, process_function, cfg, num_processes = 1):
    # Create queues
    perimeter_queue = JoinableQueue()

    i = 0
    #Insert records into queue
    for r in process_records:
        perimeter_queue.put(r)
        i+=1
        # if i>1000:
        #     break

    #Create and start worker processes
    processes = [] 
    for i in range(num_processes):
        perimeter_queue.put(None)
        proc = Process(target=worker, args=(perimeter_queue,process_function,cfg,))
        processes.append(proc)
        print("Start: ",proc)
        proc.start()

    perimeter_queue.join()

    # for p in processes:
    #     if p.exitcode == None:
    #         p.terminate()

    print("Processing finished")

# Process records in a linear fashion -> good for debugging the script
def process_records_linear(records, process_function, cfg,  num_processes = 1):
    # Create queues
    job_queue = JoinableQueue()
    result_queue = Queue()

    #Insert records into queue
    for r in records:
        job_queue.put(r)

    job_queue.put(None)

    print("Start:")
    worker(job_queue,process_function,cfg)

    print("Processing finished")
    print("Queue empty:",result_queue.empty())
    res = []
    if not result_queue.empty():
        for r in records:
            item = result_queue.get()
            res.extend(item)

    return res

#
# Business logic
#

# Setup the records describing the work units for parallel processing
def process_setup(cfg, num_processes = 1):
    
    # Read the parameter values
    result_base_path = cfg["result_base_path"]


    records = []
    # For each catchment
    for plot_path in next(os.walk(result_base_path))[1]:
        plot_name = os.path.basename(plot_path)
        plot_number = int(plot_name.replace("ezg_",""))

        # For each set of soil property parameters
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

# Worker function to setup and run simulations
def process_simulation(sim_record, cfg):

    # Get parameters for simulation from config
    output_name = sim_record["output_name"]  
    
    plot_path = sim_record["plot_path"]
    soil_properties_path = sim_record["soil_properties_path"]
    output_path = os.path.join(plot_path,output_name)

    perimeter_mask_path = os.path.join(plot_path,"perimeter_mask.tif") 

    # Create output directory
    ensure_dir(output_path)

    return_periods = [2,10,30,100,300]
    catchment_slide_count = 0
    catchment_slide_per_ha = 0

    # Run for each return period/rainfall scenario
    for rp in return_periods:
        sim_name = "rp_{0}".format(rp)
        sim_path = os.path.join(output_path,sim_name)

        prob_path = os.path.join(sim_path,"probability_failure.tif")
        
        ensure_dir(sim_path)

        # Write config file by modifying the "default_config.toml" template
        toml_template_path = os.path.join(plot_path,"default_config.toml")
        toml_output_path = os.path.join(sim_path,"artificial.toml")
        [x_min,y_min,x_max,y_max] = [0,0,0,0]
        with open(toml_template_path, mode="rt") as fp:
            config = tomlkit.load(fp)

            perimeter_cell_count = config["perimeter_cell_count"] 
            [x_min,y_min,x_max,y_max] = [config["x_min"],config["y_min"],config["x_max"],config["y_max"]]
            catchment_slide_count = catchment_slide_count if isnan(config["catchment_slide_count"]) else config["catchment_slide_count"]
            catchment_slide_per_ha = catchment_slide_per_ha if isnan(config["catchment_slide_per_ha"]) else config["catchment_slide_per_ha"]

            config["physPropCSV"] = soil_properties_path
            config["outputPath"] = sim_path
            config["rainfall"] = config["perc_{0}".format(rp)]

            config["exportLandslides"] = 0

            config["demPath"] = config["demPath"].replace("..",plot_path)
            config["flowAccPath"] = config["flowAccPath"].replace("..",plot_path)
            config["soilTypeRaster"] =  config["soilTypeRaster"].replace("..",plot_path)  
            config["soilThicknessValueRaster"] = config["soilThicknessValueRaster"].replace("..",plot_path)

            if config["treeFilePath"] != None and ".." in config["treeFilePath"]:
                config["treeFilePath"] = config["treeFilePath"].replace("..",plot_path)

            if os.path.exists(prob_path): # Simulation already exists -> skip and continue
                continue

            with open(toml_output_path, mode="wt") as fp:
                tomlkit.dump(config, fp)
    
            # Run the simulation in SlideforMAP
            cmd = "{0} --config {1}".format(cfg["sfm_path"],toml_output_path)
            os.system(cmd)

            # Mask all results to perimeter and compress
            for rast in glob.glob(os.path.join(sim_path,"*.tif")):
                tmp_raster = rast.replace(".tif","_tmp.tif")
                with gdal_calc.Calc(calc="A*B",A=rast, B=perimeter_mask_path, outfile=tmp_raster, hideNoData=True, creation_options = ["COMPRESS=DEFLATE", "BIGTIFF=YES"] , overwrite=True) as ds:
                    ds=None
                delete_raster(rast)
                os.rename(tmp_raster,rast)

    return output_path

# Entry function that initiates and coordinates the process
def process_sfm_application(cfg):
    start_time = time.time()

    print("Start Simulations",datetime.now())

    # Preparing records for parallel processing describing the catchments and parameter sets
    print("Start Setup",datetime.now())
    sim_records = process_setup(cfg, num_processes = 1)

    # Running simulations
    num_processes=cfg["num_processes"]
    #process_records_linear(sim_records, process_simulation, cfg, num_processes = num_processes)
    process_records(sim_records, process_simulation, cfg, num_processes = num_processes)
    
    print("TOTAL PROCESSING TIME: %s (h:min:sec)" % str(timedelta(seconds=(time.time() - start_time))))


# Default entry point
if __name__ == "__main__":
    freeze_support()

    # Main configuration for processing
    cfg = {
        "result_base_path": "/mnt/data/ezg_inputs", # path where the simulations will be written
        "sfm_path": "/home/administrator/slideformap/slideformap_v2.2.2/SlideforMAP", # path to SlideforMAP executable
        "output_name_template": "opt_{0}", # prefix for simulation output folders
        # Running only one simulation leveraging the SlideforMAP internal parallelization is more efficient with areas as large as the catchments
        "num_processes": 1, # number of parallel simulations
        "num_cores": 2, # number of cores per instance
        "execute_simulations": True, # Whether to actually run the simulation or just generate the setup
        "remove_outputs": False, # whether to delete "unneeded" simulation outputs  
        "observation_period": 61, # observation period for actual slides
        "slide_buffer_size": 5, # Buffer to use on historical slides
        "slide_threshold": 0.0001, # observation period for actual slides
        "raster_resolution": 2, # resolution of the probability raster
    }

    # Start the overall process
    process_sfm_application(cfg)
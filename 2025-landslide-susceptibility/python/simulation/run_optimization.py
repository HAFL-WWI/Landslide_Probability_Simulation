######################################################################
# Copyright (C) 2025 BFH
#
# Script for running optimization procedure for the soil physical 
# properties using plots with and without observed landslides.
#
# Author: Christoph Schaller, BFH-HAFL, March 2025
######################################################################

#
# Imports
#

import os

from osgeo import gdal
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
import pandas as pd

from sklearn.metrics import roc_auc_score

from lmfit import Minimizer, Parameters, create_params, report_fit

import json
import pickle

ITER_COUNT = 0

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
def process_setup(soil_properties_path, cfg, num_processes = 1):
    
    result_base_path = cfg["result_base_path"]

    # Create Records for all plots
    records = []
    for plot_path in next(os.walk(result_base_path))[1]:
        plot_name = os.path.basename(plot_path)
        plot_number = int(plot_name.replace("plot_",""))
        output_name = cfg["output_name_template"].format(ITER_COUNT)

        rec = {
            "soil_properties_path": soil_properties_path,
            "step":ITER_COUNT,
            "output_name":output_name,
            "plot_path":os.path.join(result_base_path,plot_path),
            "plot_id": plot_number,
        }
        records.append(rec)

    return records

# Variant of the function for merging propabilities from multiple rainfall scenarios into one absolute landslide probability. 
# Assumes the number of slides per hectare in slide_count for calculating the correction factor.
def merge_probabilities(output_folder,sim_folder,return_periods,observation_period,slide_count,cfg):
    """
    This function processes susceptibility TIFF files for different return periods to obtain an absolute release probability map. 
    (1) Each return period susceptibility is first weighted with the occurence probability that an event for that return period will occur during the period T (observation_period).
    (2) Susceptibilities are weighted using the occurrence probability of the return period and then summed
    (3) Result from (2) is then scaled using the number of observed landslides in the observation period

    Author: Denis Cohen-Corticchiato, Christoph Schaller
    """
    global ITER_COUNT

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
        tiff_path = os.path.join(dir_path, 'probability_failure_clipped.tif')
        
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
    occ_prob_file_name = "RLOP_sim_tot_{0}".format(ITER_COUNT)
    occ_prob_file_path = os.path.join(base_dir, occ_prob_file_name + ".tif")
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
    ) as dst:
        dst.write(weighted_sum_occurence, 1)

    # Calculate the sum of the occurence probability in the plot
    cur_rlop_sum = np.nansum(weighted_sum_occurence)

    perimeter_cell_count = weighted_sum_occurence.shape[0]*weighted_sum_occurence.shape[1]
    perimeter_area_ha = perimeter_cell_count*2*2/10000

    if slide_count == 0:
        slide_count = 1/10000

    if cur_rlop_sum!=0:
        corr_factor = (slide_count*perimeter_area_ha) / cur_rlop_sum
        print("correction_factor",corr_factor)
    else:
        corr_factor = 0.00001
    weighted_sum_release *= corr_factor

    with rio.open(occ_prob_file_path) as src:
        # Read the data
        tiff_data = src.read(1)  # Read the first band

    # Save the release probability to TIFF files
    rel_prob_file_name = "AbsLSprob_{0}".format(ITER_COUNT)
    rel_prob_file_path = os.path.join(base_dir, rel_prob_file_name + ".tif")
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
    ) as dst:
        dst.write(weighted_sum_release, 1)

    return [rel_prob_file_path,occ_prob_file_path]


# Worker function to setup and run simulations
def process_simulation(sim_record, cfg):
    global ITER_COUNT
    # Get parameters for simulation from config
    step = sim_record["step"]
    output_name = sim_record["output_name"]  
    
    plot_path = sim_record["plot_path"]
    soil_properties_path = sim_record["soil_properties_path"]
    output_path = os.path.join(plot_path,output_name)

    # Create output directory
    ensure_dir(output_path)

    observation_period = cfg["observation_period"]
    return_periods = [2,10,30,100,300]
    plot_sensitivity = False

    # Execute a simulation for wach return period/rainfall scenario
    for rp in return_periods:
        sim_name = "rp_{0}".format(rp)
        sim_path = os.path.join(output_path,sim_name)
        ensure_dir(sim_path)

        # Write config file by modifying the "default_config.toml" template
        toml_template_path = os.path.join(plot_path,"default_config.toml")
        toml_output_path = os.path.join(sim_path,"artificial.toml")
        [x_min,y_min,x_max,y_max] = [0,0,0,0]
        ps_slide_count = 0
        catchment_slide_count = 0
        catchment_slide_per_ha = 0.00001
        with open(toml_template_path, mode="rt") as fp:
            config = tomlkit.load(fp)

            plot_sensitivity = config["plot_sensitivity"]

            perimeter_cell_count = config["perimeter_cell_count"] 
            perimeter_buffered_cell_count = config["perimeter_buffered_cell_count"] 
            [x_min,y_min,x_max,y_max] = [config["x_min"],config["y_min"],config["x_max"],config["y_max"]]
            ps_slide_count = ps_slide_count if isnan(config["ps_slide_count"]) else config["ps_slide_count"]
            catchment_slide_count = catchment_slide_count if isnan(config["catchment_slide_count"]) else config["catchment_slide_count"]
            catchment_slide_per_ha = catchment_slide_per_ha if isnan(config["catchment_slide_per_ha"]) else config["catchment_slide_per_ha"]

            config["physPropCSV"] = soil_properties_path
            config["outputPath"] = sim_path
            config["rainfall"] = config["perc_{0}".format(rp)] # set return period specific rainfall amount

            config["exportLandslides"] = 0

            config["demPath"] = config["demPath"].replace("..",plot_path)
            config["flowAccPath"] = config["flowAccPath"].replace("..",plot_path)
            config["soilTypeRaster"] =  config["soilTypeRaster"].replace("soil_type_constant.tif","soil_type_slide.tif").replace("..",plot_path)  
            config["soilThicknessValueRaster"] = config["soilThicknessValueRaster"].replace("..",plot_path)

            if config["treeFilePath"] != None and ".." in config["treeFilePath"]:
                config["treeFilePath"] = config["treeFilePath"].replace("..",plot_path)

            with open(toml_output_path, mode="wt") as fp:
                tomlkit.dump(config, fp)
    
            # Run the simulation in SlideforMAP
            cmd = "{0} --config {1}".format(cfg["sfm_path"],toml_output_path)
            os.system(cmd)

            # Clip failure probability to perimeter 
            pf_path = os.path.join(sim_path,"probability_failure.tif")
            pf_clip_path = os.path.join(sim_path,"probability_failure_clipped.tif")

            with gdal.Warp(destNameOrDestDS=pf_clip_path,srcDSOrSrcDSTab=pf_path, resampleAlg=gdal.GRA_NearestNeighbour, format="GTiff", creationOptions={"COMPRESS":"DEFLATE", "BIGTIFF":"YES"}, outputBounds=[x_min,y_min,x_max,y_max]) as ds:
                ds=None

    # Merge scenario probabilities
    slide_count = (perimeter_cell_count*2*2/10000)*catchment_slide_per_ha
    [AbsLSprob_path,occ_prob_file_path] = merge_probabilities(plot_path,output_name,return_periods,observation_period,slide_count,cfg)


    # Calculate AUC
    buffer_dist = cfg["slide_buffer_size"]
    slide_buff_mask_2m_path = os.path.join(plot_path,"slide_buff{0}_mask.tif".format(buffer_dist))
    slide_buff_mask_path = os.path.join(plot_path,"slide_buff{0}_mask_clipped.tif".format(buffer_dist))
    with gdal.Warp(destNameOrDestDS=slide_buff_mask_path,srcDSOrSrcDSTab=slide_buff_mask_2m_path, resampleAlg=gdal.GRA_NearestNeighbour, format="GTiff", creationOptions={"COMPRESS":"DEFLATE", "BIGTIFF":"YES"}, outputBounds=[x_min,y_min,x_max,y_max]) as ds:
        ds=None

    slide_actual = gdal_array.LoadFile(slide_buff_mask_path).flatten()
    slide_prediction = gdal_array.LoadFile(AbsLSprob_path).flatten()

    # Use (1-probability) as comparison in plots without observed slides
    if not plot_sensitivity:
        slide_prediction = 1-slide_prediction

    auc = roc_auc_score(slide_actual,slide_prediction)
    res = 1-auc

    # Write CSV with the most important parameters and stats
    comb_path = os.path.join(output_path,"optimization_{0}.csv".format(step))
    with open(comb_path, "w") as file:
        writer = csv.writer(file, delimiter=",")

        writer.writerow(["i","auc","res"] )    
        writer.writerow([step, auc, res])  
    
    # Optionally remove unneeded outputs to save space
    if cfg["remove_outputs"]:        
        delete_raster(AbsLSprob_path)
        delete_raster(occ_prob_file_path)

        os.system("rm -rf {0}".format(os.path.join(output_path,"rp_*")))
    return output_path


# Gathering results from the simulations for the specified parameter set
def collect_results(base_path,sim_folder,cfg):#
    global ITER_COUNT
    summary_out_path = os.path.join(base_path,"optimization_merged_{0}.csv".format(ITER_COUNT))

    i=0

    # Loop all plot folders
    plotfolders = [ f.path for f in os.scandir(base_path) if f.is_dir() ]
    combs = []
    for plot_path in plotfolders:

        plot_folder = os.path.basename(plot_path) 
        plot_fid = int(plot_folder.replace("plot_",""))

        # Read CSV with result summary for current simulation
        comb_path = os.path.join(plot_path,sim_folder,"optimization_{0}.csv".format(ITER_COUNT))
        if not os.path.exists(comb_path):
            continue
        df_comb = pd.read_csv(comb_path,sep=",")
        df_comb["plot_fid"] = plot_fid
        combs.append(df_comb)

        i+=1
        if i%1000==0:
            print(i)

    # Merge results from all plots into one Dataframe
    df_merged = pd.concat(combs)
    # Sort by plot ID so it always has the same order and matches the "X" input
    df_merged = df_merged.sort_values(by=["plot_fid"])

    # Save the merged result as CSV
    df_merged.to_csv(summary_out_path, sep=",")

    # Get the Result vecton (1-AUC) to return to the optimizer
    ret = df_merged["res"].to_list()
    return ret



def fit_soil_type_parameters(pars, x, cfg=None):
    global ITER_COUNT

    # Default soil parameters
    vals = dict(
        s_id_1 = 1, 		por_1 = 0.38,	fc_1 = 0.1,	r_1 = 36,	r_sd_1 = 0,	c_1 = 0,	c_sd_1 = 0,	ks_1 = 8600,	ks_sd_1 = 0,	alpha_1 = 0.001,	n_1 = 2,	d50_1 = 1,	 percentClaySilt_1 =5,	
        s_id_2 = 2, 		por_2 = 0.37,	fc_2 = 0.16,	r_2 = 34,	r_sd_2 = 0,	c_2 = 0,	c_sd_2 = 0,	ks_2 = 20,	ks_sd_2 = 0,	alpha_2 = 0.001,	n_2 = 2,	d50_2 = 1,	 percentClaySilt_2 =25,	
        s_id_3 = 3, 		por_3 = 0.37,	fc_3 = 0.21,	r_3 = 34,	r_sd_3 = 0,	c_3 = 0.5,	c_sd_3 = 0,	ks_3 = 20,	ks_sd_3 = 0,	alpha_3 = 0.001,	n_3 = 2,	d50_3 = 1,	 percentClaySilt_3 =50,	
        s_id_4 = 4, 		por_4 = 0.4,	fc_4 = 0.36,	r_4 = 32,	r_sd_4 = 0,	c_4 = 1.5,	c_sd_4 = 0,	ks_4 = 20,	ks_sd_4 = 0,	alpha_4 = 0.001,	n_4 = 2,	d50_4 = 1,	 percentClaySilt_4 =50,	
        s_id_5 = 5, 		por_5 = 0.4,	fc_5 = 0.32,	r_5 = 32,	r_sd_5 = 0,	c_5 = 0,	c_sd_5 = 0,	ks_5 = 20,	ks_sd_5 = 0,	alpha_5 = 0.001,	n_5 = 2,	d50_5 = 1,	 percentClaySilt_5 =50,	
        s_id_6 = 6, 		por_6 = 0.41,	fc_6 = 0.27,	r_6 = 27,	r_sd_6 = 0,	c_6 = 2,	c_sd_6 = 0,	ks_6 = 20,	ks_sd_6 = 0,	alpha_6 = 0.001,	n_6 = 2,	d50_6 = 1,	 percentClaySilt_6 =50,	
        s_id_7 = 7, 		por_7 = 0.47,	fc_7 = 0.3,	r_7 = 33,	r_sd_7 = 0,	c_7 = 0,	c_sd_7 = 0,	ks_7 = 20,	ks_sd_7 = 0,	alpha_7 = 0.001,	n_7 = 2,	d50_7 = 1,	 percentClaySilt_7 =70,	
        s_id_8 = 8, 		por_8 = 0.47,	fc_8 = 0.29,	r_8 = 33,	r_sd_8 = 0,	c_8 = 0,	c_sd_8 = 0,	ks_8 = 20,	ks_sd_8 = 0,	alpha_8 = 0.001,	n_8 = 2,	d50_8 = 1,	 percentClaySilt_8 =90,	
        s_id_9 = 9, 		por_9 = 0.41,	fc_9 = 0.29,	r_9 = 27,	r_sd_9 = 0,	c_9 = 2,	c_sd_9 = 0,	ks_9 = 20,	ks_sd_9 = 0,	alpha_9 = 0.001,	n_9 = 2,	d50_9 = 1,	 percentClaySilt_9 =70,	
        s_id_10 = 10, 		por_10 = 0.41,	fc_10 = 0.28,	r_10 = 27,	r_sd_10 = 0,	c_10 = 2,	c_sd_10 = 0,	ks_10 = 20,	ks_sd_10 = 0,	alpha_10 = 0.001,	n_10 = 2,	d50_10 = 1,	 percentClaySilt_10 =80,	
        s_id_11 = 11, 		por_11 = 0.56,	fc_11 = 0.22,	r_11 = 22,	r_sd_11 = 0,	c_11 = 2.5,	c_sd_11 = 0,	ks_11 = 20,	ks_sd_11 = 0,	alpha_11 = 0.001,	n_11 = 2,	d50_11 = 1,	 percentClaySilt_11 =70,	
        s_id_12 = 12, 		por_12 = 0.56,	fc_12 = 0.2,	r_12 = 22,	r_sd_12 = 0,	c_12 = 2.5,	c_sd_12 = 0,	ks_12 = 20,	ks_sd_12 = 0,	alpha_12 = 0.001,	n_12 = 2,	d50_12 = 1,	 percentClaySilt_12 =80,	
        s_id_13 = 13, 		por_13 = 0,	fc_13 = 0.17,	r_13 = 90,	r_sd_13 = 0,	c_13 = 999,	c_sd_13 = 0,	ks_13 = 20000,	ks_sd_13 = 0,	alpha_13 = 0.001,	n_13 = 2,	d50_13 = 1,	 percentClaySilt_13 =50,	
    )

    vals_pars = pars.valuesdict()

    vals.update(vals_pars)

    ITER_COUNT+=1

    # Generating soil property csv for this iteration based on parameters
    result_base_path = cfg["result_base_path"]
    comb_path = os.path.join(result_base_path,"combination_{0}.csv".format(ITER_COUNT))
    soil_properties_path = os.path.join(result_base_path,"SoilProperties_{0}.csv".format(ITER_COUNT))
    with open(soil_properties_path, "w") as file:
        writer = csv.writer(file, delimiter=",")

        writer.writerow(["Soil ID","Soil Name","Porosity","Field capacity","Friction angle mean (deg)","Friction angle std dev (deg)","Cohesion mean (kPa)","Cohesion std dev (kPa)","Ks mean (m/d)","Ks std dev (m/d)","alpha","n","D50 (mm)","percentClaySilt"] )    
        shifts = [0,20]
        soil_types = [1,2,3,4,5,6,7,8,9,10,11,12,13]
        for shift in shifts:
            for st in soil_types:
                s_id = int(vals["s_id_{0}".format(st)] + shift)
                por_mean = vals["por_{0}".format(st)]
                fc_mean = vals["fc_{0}".format(st)]
                r_mean = vals["r_{0}".format(st)]
                r_sd = vals["r_sd_{0}".format(st)]
                c_mean = vals["c_{0}".format(st)]
                c_sd = vals["c_sd_{0}".format(st)]
                ks_mean = vals["ks_{0}".format(st)]
                ks_sd = vals["ks_sd_{0}".format(st)]
                alpha = vals["alpha_{0}".format(st)]
                n = vals["n_{0}".format(st)]
                d50 = vals["d50_{0}".format(st)]
                percentClaySilt = vals["percentClaySilt_{0}".format(st)]

                ks_mean = ks_mean if shift==0 or st==13 else ks_mean+700
                writer.writerow([s_id,s_id,por_mean,fc_mean,r_mean,r_sd,c_mean,c_sd,ks_mean,ks_sd,alpha,n,d50,percentClaySilt])  

    # Preparing records for parallel processing
    print("Start Setup",datetime.now())
    sim_records = process_setup(soil_properties_path, cfg, num_processes = 1)
    print(len(sim_records))
    # print(sim_records)

    # Running simulations
    num_processes=cfg["num_processes"]

    #process_records_linear(sim_records, process_simulation, cfg, num_processes = num_processes)
    process_records(sim_records, process_simulation, cfg, num_processes = num_processes)

    # Gathering results from Simulation
    result = collect_results(result_base_path,cfg["output_name_template"].format(ITER_COUNT),cfg)
    # print("collected result",result)

    return result

# Entry function that initiates and coordinates the process
def process_sfm_optimization(cfg):
    start_time = time.time()

    print("Start Simulations",datetime.now())

    # determine the plot numbers as "X" input for the optimizer
    result_base_path = cfg["result_base_path"]
    x = []
    plot_folders = []
    for plot_path in next(os.walk(result_base_path))[1]:
        plot_name = os.path.basename(plot_path)
        plot_folders.append(plot_name)
        plot_number = int(plot_name.replace("plot_",""))
        x.append(plot_number)
    x.sort()

    # Define parameters for the optimization procedure including bounds
    params = create_params(

        # # # middle initial guess
        # # r_1 = dict(value=36, min=34, max=38, vary=True),	c_1 = dict(value=0.5, min=0, max=1, vary=True),	ks_1 = dict(value=4000, min=80, max=9000, vary=True),
        # r_3 = dict(value=31, min=29, max=34, vary=True),	c_3 = dict(value=2.5, min=0, max=5, vary=True),	ks_3 = dict(value=55, min=10, max=100, vary=True),
        # r_4 = dict(value=32, min=28, max=36, vary=True),	c_4 = dict(value=2.5, min=0, max=5, vary=True),	ks_4 = dict(value=40, min=0, max=80, vary=True),
        # # r_5 = dict(value=32, min=28, max=36, vary=True),	c_5 = dict(value=2.5, min=0, max=5, vary=True),	ks_5 = dict(value=40, min=0, max=80, vary=True),
        # r_7 = dict(value=31, min=28, max=34, vary=True),	c_7 = dict(value=2.5, min=0, max=5, vary=True),	ks_7 = dict(value=40, min=0, max=80, vary=True),
        # r_8 = dict(value=33, min=29, max=37, vary=True),	c_8 = dict(value=2.5, min=0, max=5, vary=True),	ks_8 = dict(value=40, min=0, max=80, vary=True),
        # # r_9 = dict(value=27, min=23, max=31, vary=True),	c_9 = dict(value=3, min=1, max=5, vary=True),	ks_9 = dict(value=40, min=0, max=80, vary=True),
        # # r_10 = dict(value=28, min=24, max=32, vary=True),	c_10 = dict(value=3, min=1, max=5, vary=True),	ks_10 = dict(value=40, min=0, max=80, vary=True),
        # # r_11 = dict(value=24, min=20, max=28, vary=True),	c_11 = dict(value=5.5, min=1, max=10, vary=True),	ks_11 = dict(value=40, min=0, max=80, vary=True),
        # # r_12 = dict(value=22, min=18, max=26, vary=True),	c_12 = dict(value=5.5, min=1, max=10, vary=True),	ks_12 = dict(value=40, min=0, max=80, vary=True),


        # # expert input initial guess
        # r_1 = dict(value=36, min=34, max=38, vary=True),	c_1 = dict(value=0, min=0, max=1, vary=True),	ks_1 = dict(value=8600, min=80, max=9000, vary=True),
        r_3 = dict(value=34, min=29, max=34, vary=True),	c_3 = dict(value=0.5, min=0, max=5, vary=True),	ks_3 = dict(value=20, min=10, max=100, vary=True),
        r_4 = dict(value=32, min=28, max=36, vary=True),	c_4 = dict(value=1.5, min=0, max=5, vary=True),	ks_4 = dict(value=20, min=0, max=80, vary=True),
        # r_5 = dict(value=32, min=28, max=36, vary=True),	c_5 = dict(value=0, min=0, max=5, vary=True),	ks_5 = dict(value=20, min=0, max=80, vary=True),
        r_7 = dict(value=33, min=28, max=34, vary=True),	c_7 = dict(value=0, min=0, max=5, vary=True),	ks_7 = dict(value=20, min=0, max=80, vary=True),
        r_8 = dict(value=33, min=29, max=37, vary=True),	c_8 = dict(value=0, min=0, max=5, vary=True),	ks_8 = dict(value=20, min=0, max=80, vary=True),
        # r_9 = dict(value=27, min=23, max=31, vary=True),	c_9 = dict(value=2, min=1, max=5, vary=True),	ks_9 = dict(value=20, min=0, max=80, vary=True),
        # r_10 = dict(value=27, min=24, max=32, vary=True),	c_10 = dict(value=2, min=1, max=5, vary=True),	ks_10 = dict(value=20, min=0, max=80, vary=True),
        # r_11 = dict(value=22, min=20, max=28, vary=True),	c_11 = dict(value=2.5, min=1, max=10, vary=True),	ks_11 = dict(value=20, min=0, max=80, vary=True),
        # r_12 = dict(value=22, min=18, max=26, vary=True),	c_12 = dict(value=2.5, min=1, max=10, vary=True),	ks_12 = dict(value=20, min=0, max=80, vary=True),


        # # minimum initial guess
        # # r_1 = dict(value=34, min=34, max=38, vary=True),	c_1 = dict(value=0, min=0, max=1, vary=True),	ks_1 = dict(value=80, min=80, max=9000, vary=True),
        # r_3 = dict(value=29, min=29, max=34, vary=True),	c_3 = dict(value=0, min=0, max=5, vary=True),	ks_3 = dict(value=0, min=10, max=100, vary=True),
        # r_4 = dict(value=28, min=28, max=36, vary=True),	c_4 = dict(value=0, min=0, max=5, vary=True),	ks_4 = dict(value=0, min=0, max=80, vary=True),
        # # r_5 = dict(value=28, min=28, max=36, vary=True),	c_5 = dict(value=0, min=0, max=5, vary=True),	ks_5 = dict(value=0, min=0, max=80, vary=True),
        # r_7 = dict(value=28, min=28, max=34, vary=True),	c_7 = dict(value=0, min=0, max=5, vary=True),	ks_7 = dict(value=0, min=0, max=80, vary=True),
        # r_8 = dict(value=29, min=29, max=37, vary=True),	c_8 = dict(value=0, min=0, max=5, vary=True),	ks_8 = dict(value=0, min=0, max=80, vary=True),
        # # r_9 = dict(value=23, min=23, max=31, vary=True),	c_9 = dict(value=1, min=1, max=5, vary=True),	ks_9 = dict(value=0, min=0, max=80, vary=True),
        # # r_10 = dict(value=24, min=24, max=32, vary=True),	c_10 = dict(value=1, min=1, max=5, vary=True),	ks_10 = dict(value=0, min=0, max=80, vary=True),
        # # r_11 = dict(value=20, min=20, max=28, vary=True),	c_11 = dict(value=1, min=1, max=10, vary=True),	ks_11 = dict(value=0, min=0, max=80, vary=True),
        # # r_12 = dict(value=18, min=18, max=26, vary=True),	c_12 = dict(value=1, min=1, max=10, vary=True),	ks_12 = dict(value=0, min=0, max=80, vary=True),

        # # maximum initial guess
        # # r_1 = dict(value=38, min=34, max=38, vary=True),	c_1 = dict(value=1, min=0, max=1, vary=True),	ks_1 = dict(value=9000, min=80, max=9000, vary=True),
        # r_3 = dict(value=34, min=29, max=34, vary=True),	c_3 = dict(value=5, min=0, max=5, vary=True),	ks_3 = dict(value=100, min=10, max=100, vary=True),
        # r_4 = dict(value=36, min=28, max=36, vary=True),	c_4 = dict(value=5, min=0, max=5, vary=True),	ks_4 = dict(value=80, min=0, max=80, vary=True),
        # # r_5 = dict(value=36, min=28, max=36, vary=True),	c_5 = dict(value=5, min=0, max=5, vary=True),	ks_5 = dict(value=80, min=0, max=80, vary=True),
        # r_7 = dict(value=34, min=28, max=34, vary=True),	c_7 = dict(value=5, min=0, max=5, vary=True),	ks_7 = dict(value=80, min=0, max=80, vary=True),
        # r_8 = dict(value=37, min=29, max=37, vary=True),	c_8 = dict(value=5, min=0, max=5, vary=True),	ks_8 = dict(value=80, min=0, max=80, vary=True),
        # # r_9 = dict(value=31, min=23, max=31, vary=True),	c_9 = dict(value=5, min=1, max=5, vary=True),	ks_9 = dict(value=80, min=0, max=80, vary=True),
        # # r_10 = dict(value=32, min=24, max=32, vary=True),	c_10 = dict(value=5, min=1, max=5, vary=True),	ks_10 = dict(value=80, min=0, max=80, vary=True),
        # # r_11 = dict(value=28, min=20, max=28, vary=True),	c_11 = dict(value=10, min=1, max=10, vary=True),	ks_11 = dict(value=80, min=0, max=80, vary=True),
        # # r_12 = dict(value=26, min=18, max=26, vary=True),	c_12 = dict(value=10, min=1, max=10, vary=True),	ks_12 = dict(value=80, min=0, max=80, vary=True),

    )
    
    #
    # Minimize with leastsq algorithm (Levenberg-Marquardt)
    #

    # Set the parameters for the minimizer
  
    # minner = Minimizer(fit_soil_type_parameters, params, max_nfev=10000, fcn_args=(x,cfg))
    minner = Minimizer(fit_soil_type_parameters, params, max_nfev=10000, fcn_args=(x,cfg),epsfcn=0.01)
    # minner = Minimizer(fit_soil_type_parameters, params, max_nfev=10000, fcn_args=(x,cfg),epsfcn=0.1)
    # minner = Minimizer(fit_soil_type_parameters, params, max_nfev=10000, fcn_args=(x,cfg),epsfcn=0.25)
    # minner = Minimizer(fit_soil_type_parameters, params, max_nfev=10000, fcn_args=(x,cfg),epsfcn=0.5)
    # minner = Minimizer(fit_soil_type_parameters, params, max_nfev=10000, fcn_args=(x,cfg),epsfcn=0.75)

    # Start the minimization procedure 
    result = minner.minimize()



    #
    # Minimize with Nelder-Mead algorithm using an initial simplex based on expert input, min, middle, 
    # and max start values as well as on results from Levenberg-Marquardt optimizations
    #

    # # Set initial simplex
    # init_simplex = [
    #     [34.0,0.0,20.0, 32.0,1.5,20.0, 33.0,0.0,20.0,  33.0,0.0,20.0],
    #     [29.0,0.0,10.0, 28.0,0.0,20.0, 28.0,0.0,20.0,  29.0,0.0,20.0],
    #     [34.0,5.0,100.0, 36.0,5.0,80.0, 34.0,5.0,80.0,  37.0,5.0,80.0],
    #     [31.5,2.5,60.0, 32.0,2.5,50.0, 31.0,2.5,50.0,  33.0,2.5,50.0],
    #     [32.0,4.0,40.00, 32.0,3.0,10.0, 31.0,1.0,10.0, 32.0,2.50,30.0],
    #     [32.54,1.71,40.82, 30.26,2.36,8.21,	28.00,3.77,47.48, 30.15,2.58,11.44],
    #     [30.61,4.35,29.97, 32.25,2.92,26.64, 28.54,1.13,0.45, 32.98,2.31,26.02],
    #     [33.03,3.28,10.16, 28.00,1.86,7.59, 33.99,0.22,1.06, 32.67,2.23,4.98],
    #     [33.80,4.86,99.85, 31.35,1.15,8.90, 34.0,5.0,80.0, 35.57,0.5,28.89],
    #     [31.55,4.14,88.00, 29.79,1.74,5.33, 31.30,1.30,15.06, 32.46,2.26,6.74],
    #     [30.61,4.35,29.97, 32.25,2.92,26.64, 28.54,1.13,0.45, 32.98,2.31,26.02],
    #     [29.89,1.77,20.59, 31.46,1.77,14.90, 28.80,0.06,0.44, 32.24,2.34,4.68],
    #     [33.99,0.52,20.00, 32.01,1.53,19.99, 32.99,0.00,19.99, 33.01,0.00,19.99],
    # ]
    #
    # # Set the parameters for the minimizer
    # minner = Minimizer(fit_soil_type_parameters, params, fcn_args=(x,cfg), options={"initial_simplex":init_simplex})
    #     
    # # Start the minimization procedure 
    # result = minner.minimize(method="nelder")


    # Write error report to console
    report_fit(result)

    # Save the optimization result by pickling it (may not work with all optimization algorithms)
    model_result_pickle_path = os.path.join(result_base_path,"model_fit_reslt.pkl")

    with open(model_result_pickle_path, 'wb') as f:
        pickle.dump(result, f)
    
    print("TOTAL PROCESSING TIME: %s (h:min:sec)" % str(timedelta(seconds=(time.time() - start_time))))


# Default entry point
if __name__ == "__main__":
    freeze_support()

    # Main configuration for processing
    cfg = {
        "result_base_path": "/mnt/data/optimization_inputs", # path where the simulations will be written
        "sfm_path": "/home/administrator/slideformap/slideformap_v2.2.2/SlideforMAP", # path to SlideforMAP executable
        "output_name_template": "opt_{0}", # prefix for simulation output folders
        "num_processes": 40, # number of parallel simulations
        "num_cores": 2, # number of cores per instance
        "execute_simulations": True, # Whether to actually run the simulation or just generate the setup
        "remove_outputs": True, # whether to delete "unneeded" simulation outputs  
        "observation_period": 61, # observation period for actual slides
        "slide_buffer_size": 5, # Buffer to use on historical slides
        "slide_threshold": 0.0001, # observation period for actual slides
    }

    # Start the overall process
    process_sfm_optimization(cfg)
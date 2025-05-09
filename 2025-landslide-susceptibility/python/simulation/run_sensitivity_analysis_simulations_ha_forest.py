######################################################################
# Copyright (C) 2025 BFH
#
# Script for running the simulations for the sensitivity analysis 
# for forest effect using the plots with observed landslides and  with 
# at least 50% forest cover.
#
# Author: Christoph Schaller, BFH-HAFL, March 2025
######################################################################

#
# Imports
#

import os
import glob

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
import geopandas as gpd

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
def process_setup(plot_path, cfg, num_processes = 1):
    
    # Read the parameter values
    result_base_path = cfg["result_base_path"]
    param_in_path = os.path.join(result_base_path,"sensitivity_param_values_with_forest.txt")
    param_values = np.loadtxt(param_in_path) 

    records = []

    # Check if plot is intended for sensitivity analysis
    plot_sensitivity = False
    [x_min,y_min,x_max,y_max] = [0,0,0,0]
    perimeter_cell_count = 0
    forest_cell_count = 0 
    toml_template_path = os.path.join(plot_path,"default_config.toml")
    with open(toml_template_path, mode="rt") as fp:
        config = tomlkit.load(fp)
        [x_min,y_min,x_max,y_max] = [config["x_min"],config["y_min"],config["x_max"],config["y_max"]]
        
        perimeter_cell_count = config["perimeter_cell_count"] 
        forest_cell_count = config["forest_cell_count"] 

        plot_sensitivity = config["plot_sensitivity"] 

    # Determine forest cover in percent and set indicator flag 
    plot_forest = False
    if forest_cell_count/perimeter_cell_count>=0.5:
        plot_forest = True

    combination_path = os.path.join(plot_path,"combination_forest.csv")
    if os.path.exists(combination_path):
        return records

    # Only process plots with observed slides and >50% forest cover
    if plot_sensitivity and plot_forest:
       
        # Generate "records" for all plots and parameter combinations that will be iterated over when executing the simulations
        for i in range(len(param_values)):
            [c_mean,fc_mean, ks_mean, por_mean, r_mean,rainfall,forest] = param_values[i]

            output_name = cfg["output_name_template"].format(i)
            rec = {
                "c_mean": c_mean,"fc_mean": fc_mean, "ks_mean": ks_mean, "por_mean": por_mean, "r_mean": r_mean, "rainfall": rainfall,"simulate_with_forest": forest==1,
                "step":i,
                "output_name":output_name,
                "plot_path":plot_path,
            }
            records.append(rec)

    return records

# Worker function to setup and run simulations
def process_simulation(sim_record, cfg):
    # Get parameters for simulation from config
    step = sim_record["step"]
    output_name = sim_record["output_name"]  
    
    plot_path = sim_record["plot_path"]
    output_path = os.path.join(plot_path,output_name)

    if os.path.exists(output_path):
        comb_path = os.path.join(output_path,cfg["combination_name"])
        if os.path.exists(comb_path):
            return output_path
    # Create output directory
    ensure_dir(output_path)

    c_mean = sim_record["c_mean"]
    fc_mean = sim_record["fc_mean"]
    ks_mean = sim_record["ks_mean"]
    por_mean = sim_record["por_mean"]
    r_mean = sim_record["r_mean"]
    rainfall = sim_record["rainfall"]



    # Generating soil property csv
    soil_properties_path = os.path.join(output_path,"SoilProperties_uniform.csv")
    with open(soil_properties_path, "w") as file:
        writer = csv.writer(file, delimiter=",")

        writer.writerow(["Soil ID","Soil Name","Porosity","Field capacity","Friction angle mean (deg)","Friction angle std dev (deg)","Cohesion mean (kPa)","Cohesion std dev (kPa)","Ks mean (m/d)","Ks std dev (m/d)","alpha","n","D50 (mm)","percentClaySilt"] )    
        writer.writerow([1,"Sensitivity",por_mean,fc_mean,r_mean,0,c_mean,0,ks_mean,0,0.001,2,1,50])  
    

    # Write config file by modifying the "default_config.toml" template
    toml_template_path = os.path.join(plot_path,"default_config.toml")
    toml_output_path = os.path.join(output_path,"artificial.toml")
    perimeter_cell_count = 0
    perimeter_buffered_cell_count = 0
    forest_cell_count = 0
    [x_min,y_min,x_max,y_max] = [0,0,0,0]
    with open(toml_template_path, mode="rt") as fp:
        config = tomlkit.load(fp)

        perimeter_cell_count = config["perimeter_cell_count"] 
        perimeter_buffered_cell_count = config["perimeter_buffered_cell_count"] 
        forest_cell_count = config["forest_cell_count"] 
        [x_min,y_min,x_max,y_max] = [config["x_min"],config["y_min"],config["x_max"],config["y_max"]]

        config["physPropCSV"] = soil_properties_path
        config["outputPath"] = output_path
        config["rainfall"] = rainfall


        config["demPath"] = config["demPath"].replace("..",plot_path)
        config["flowAccPath"] = config["flowAccPath"].replace("..",plot_path)
        config["soilTypeRaster"] =  config["soilTypeRaster"].replace("..",plot_path)  
        config["soilThicknessValueRaster"] = config["soilThicknessValueRaster"].replace("..",plot_path)

        if sim_record["simulate_with_forest"]:
            if config["treeFilePath"] != None and ".." in config["treeFilePath"]:
                config["treeFilePath"] = config["treeFilePath"].replace("..",plot_path)
        else:
            config.pop("treeFilePath")

        with open(toml_output_path, mode="wt") as fp:
            tomlkit.dump(config, fp)


    # Execute the actual simulation
    if cfg["execute_simulations"]: 
        #os.system("export OMP_NUM_THREADS={0}".format(cfg["num_cores"])) # Needs to be set outside the python instance
        
        # Run the simulation in SlideforMAP
        cmd = "{0} --config {1}".format(cfg["sfm_path"],toml_output_path)
        os.system(cmd)
    else: # Stop at setup and return
        return output_path
    
    # Convert simulation outputs to compressed GeoTIFF to save space
    for rast_in in glob.glob(os.path.join(output_path,"*.tif")):
        rast_out = rast_in.replace(".tif","_conf.tif")
        with gdal.Translate(destName=rast_out, srcDS=rast_in, format="GTiff", creationOptions = ["COMPRESS=DEFLATE"]) as ds:
            ds=None
        delete_raster(rast_in)
        os.rename(rast_out,rast_in)

    # Aggregate the exported landslides by 0.01 steps of the FOS
    ls_path = os.path.join(output_path,"ls_export.csv")
    ls_zip_path = os.path.join(output_path,"ls_export.tar.bz2")
    os.system("tar -cjf {1} {0}".format(ls_path,ls_zip_path))


    # Determine number of failed landslides in vore plot
    ls_df = pd.read_csv(ls_path,sep=",")
    ls_failed_count = 0
    ls_failed_all_count = 0
    ls_failed_forest = 0
    if len(ls_df)>0:
        ls_failed_all_count = len(ls_df)
        ls_failed_count = len(ls_df.loc[(ls_df["X (m)"]>x_min) & (ls_df["X (m)"]<=x_max) & (ls_df[" Y (m)"]>y_min) & (ls_df[" Y (m)"]<=y_max)])

        forest_mask_gpkg = os.path.join(plot_path, "forest_mask.gpkg")
        forest_mask_layer = "forest_mask"

        forest_gdf = gpd.GeoDataFrame.from_file(forest_mask_gpkg,layer=forest_mask_layer,engine="pyogrio",)
        ls_gdf = gpd.GeoDataFrame(ls_df, geometry=gpd.points_from_xy(ls_df["X (m)"], ls_df[" Y (m)"]), crs="EPSG:2056")
        ls_gdf_clip = gpd.clip(ls_gdf, forest_gdf)
        ls_failed_forest = len(ls_gdf_clip)

    # Clip failure probability to perimeter 
    pf_path = os.path.join(output_path,"probability_failure.tif")
    pf_clip_path = os.path.join(output_path,"probability_failure_clipped.tif")
    perimeter_mask_path = os.path.join(plot_path,"perimeter_mask.tif")
    
    # Mask to perimeter while multiplying by 10000 and converting to INT to save space 
    with gdal_calc.Calc(calc="A*10000*B",A=pf_path, B=perimeter_mask_path, type=GDT_Int16, outfile=pf_clip_path, creation_options = ["COMPRESS=DEFLATE"] , overwrite=True) as ds:
        ds=None
    
    # Determine Stats from the failure probability
    [pf_min, pf_max, pf_mean, pf_5,p_10,pf_q1, pf_median, pf_q2,pf_90,pf_95] =  [-1,-1,-1, -1,-1,-1,-1,-1,-1,-1]

    perimeter_gpkg = os.path.join(plot_path,"perimeter.gpkg")
    
    with rio.open(pf_clip_path, "r") as src:
        array = src.read().astype(np.float16)
        array_selection = array[array>=0]
        pf_min = np.min(array_selection) 
        pf_max = np.max(array_selection) 
        pf_mean = np.mean(array_selection) 
        [pf_5,p_10,pf_q1, pf_median, pf_q2,pf_90,pf_95] = np.nanpercentile(array_selection,[5,10,25,50,75,90,95]) 

    # Write CSV with the most important parameters and stats
    comb_path = os.path.join(output_path,cfg["combination_name"])
    with open(comb_path, "w") as file:
        writer = csv.writer(file, delimiter=",")

        writer.writerow(["i","c","fc","ks", "por", "r","rainfall","min","max","mean","p5","p10","q1","median","q2","p90","p95","with_forest","cells","cells_buffered","cells_forest","ls_failed","ls_failed_all","ls_failed_forest"] )    
        writer.writerow([step,c_mean,fc_mean, ks_mean, por_mean, r_mean,rainfall,pf_min, pf_max, pf_mean, pf_5,p_10,pf_q1, pf_median, pf_q2,pf_90,pf_95, 1 if sim_record["simulate_with_forest"] else 0,perimeter_cell_count, perimeter_buffered_cell_count,forest_cell_count, ls_failed_count, ls_failed_all_count,ls_failed_forest])  

    # Optionally remove unneeded outputs to save space
    if cfg["remove_outputs"]:
        for f in ["friction_angle_deg.tif","cohesion_kPa.tif","hydraul_conduct_ms-1.tif","probability_runout.tif","soil_thickness_m.tif","water_pressure_kPa.tif","root_basal_kPa.tif", "root_lateral_kNm-1.tif", "wetness_index.tif", "ls_export.csv","artificial.toml"]:
            f_path = os.path.join(output_path,f)
            if os.path.isfile(f_path):
                os.remove(f_path)

    return output_path

# Function that collects the result CSV from the individual simulations of a plot into a single CSV
def collect_results(plot_path,cfg):
    summary_out_path = os.path.join(plot_path,cfg["combination_name"].replace(".csv","_merged.csv"))

    plot_folder = os.path.basename(plot_path)
    print(plot_folder)

    # Open CSV for saving results
    i=0 # Counter for row number
    with open( summary_out_path, "w", newline="") as comb_file:
        comb_writer = csv.writer(comb_file, delimiter=",")

        # Write header
        comb_writer.writerow(["plot_fid","i","c","fc","ks", "por", "r","rainfall","min","max","mean","p5","p10","q1","median","q2","p90","p95","with_forest","cells","cells_buffered","cells_forest","ls_failed","ls_failed_all","ls_failed_forest"] )

        plot_fid = int(plot_folder.replace("plot_",""))

        # Loop all simulation folders
        simfolders = [ f.path for f in os.scandir(plot_path) if f.is_dir() ]
        for sim in simfolders:

            sim_folder = os.path.basename(sim)
            sim_nr = int(sim_folder.replace("sim_",""))

            # Read csv with simulation result summary for current simulation
            comb_path = os.path.join(sim,cfg["combination_name"])
            df_comb = pd.read_csv(comb_path,sep=',')

            # Write result to central CSV
            comb_writer.writerow([plot_fid,df_comb["i"].iloc[0],df_comb["c"].iloc[0],df_comb["fc"].iloc[0],df_comb["ks"].iloc[0], df_comb["por"].iloc[0], df_comb["r"].iloc[0],df_comb["rainfall"].iloc[0],
                                df_comb["min"].iloc[0],df_comb["max"].iloc[0],df_comb["mean"].iloc[0],df_comb["p5"].iloc[0],df_comb["p10"].iloc[0],df_comb["q1"].iloc[0],
                                df_comb["median"].iloc[0],df_comb["q2"].iloc[0],df_comb["p90"].iloc[0],df_comb["p95"].iloc[0],df_comb["with_forest"].iloc[0] ,df_comb["cells"].iloc[0],
                                df_comb["cells_buffered"].iloc[0],df_comb["cells_forest"].iloc[0],df_comb["ls_failed"].iloc[0],df_comb["ls_failed_all"].iloc[0],df_comb["ls_failed_forest"].iloc[0]])  

            # Update row counter and print progress occasionally 
            i+=1
            if i%1000==0:
                print(i)

    # Save the actual simulation outputs as a TAR archive (SLOW)
    if cfg["tar_outputs"]:
        plot_tar_path = os.path.join(cfg["result_base_path"],"{0}.tar".format(plot_folder))
        os.system("tar -cf {1} {0}".format(plot_path,plot_tar_path))

    # Delete the folder with the individual simulation results
    os.system("rm -rf {0}".format(os.path.join(plot_path,"sim*")))

# Entry function that initiates and coordinates the process
def process_sfm_simulations(cfg):
    start_time = time.time()

    print("Start Simulations",datetime.now())

    result_base_path = cfg["result_base_path"]

    # Get all folders for the plots and run simulations over all combinations per plot
    for folder_name in next(os.walk(result_base_path))[1]:
        plot_path = os.path.join(result_base_path,folder_name)
        print(plot_path)

        # Simulations have already been executed -> Skip plot
        if os.path.exists(os.path.join(plot_path,cfg["combination_name"].replace(".csv","_merged.csv"))):
            print("Exists")
            continue

        # Generate records with work units for all combinations
        print("Start Setup",datetime.now())
        sim_records = process_setup(plot_path,cfg)

        # No records -> skip processing
        print(len(sim_records))
        if len(sim_records)==0:
            continue

        # Execute the actual simulations
        num_processes=cfg["num_processes"]

        #process_records_linear(sim_records, process_simulation, cfg, num_processes = num_processes)
        process_records(sim_records, process_simulation, cfg, num_processes = num_processes)

        # Collect the result from all simulations to a single CSV for the plot
        collect_results(plot_path,cfg)

    print("TOTAL PROCESSING TIME: %s (h:min:sec)" % str(timedelta(seconds=(time.time() - start_time))))



# Default entry point
if __name__ == "__main__":
    freeze_support()

    # Main configuration for processing
    cfg = {
        "result_base_path": "/mnt/data/analysis_inputs_forest", # path where the simulations will be written
        "sfm_path": "/home/administrator/slideformap/slideformap_v2.2.2/SlideforMAP", # path to SlideforMAP executable
        "output_name_template": "sim_{0}", # prefix for simulation output folders
        "num_processes": 40, # number of parallel simulations
        "num_cores": 2, # number of cores per instance
        "combination_name": "combination_forest.csv", # filename for colected output
        "execute_simulations": True, # Whether to actually run the simulation or just generate the setup
        "remove_outputs": True, # whether to delete "unneeded" simulation outputs  
        "tar_outputs": False # whether to tar the simulation outputs before deleting the simulation folders 
    }

    # Start the overall process
    process_sfm_simulations(cfg)


######################################################################
# Copyright (C) 2025 BFH
#
# Script for generating a prediction raster for the failure depth
# based on the previously trained random forests model. The previously
# prepared covariate rasters from training are needed as input. The 
# prediction rasters are generated in tiles based on aggregated 
# catchments of about 150 km2 area and then merged.
#
# Author: Christoph Schaller, BFH-HAFL, March 2025
######################################################################

#
# Initial configuration
#

gdal_base_path <- "C:/OSGeo4W/bin"
Sys.setenv(GDAL_DATA=gdal_base_path)
Sys.setenv(GDAL_DRIVER_PATH=gdal_base_path)
# Sys.setenv(PROJ_PATH=file.path(sf_path,"proj"))
Sys.getenv("GDAL_DATA")

library(sf)
library(gdalUtilities)
library(terra)

library(RSAGA)
library(dplyr)
library(parallel)

#Paths, names and other configs
use_smooth_rasters <- FALSE

base_path <- "E:/GIS_Projekte/Paper_2"

watershed_path <- "../EZG_Gewaesser_aggregiert_150.gpkg"
watershed_layer <- "EZG_Gewaesser_aggregiert_150"
watershed_mask_path <- "EZG_TEZGNR150_mask.tif"

perimeter_path  <- "E:/GIS_Projekte/GHK/BE/BE_Perimeter.gpkg"
perimeter_layer <- "BE_Perimeter"

output_base_path <- "E:/GIS_Projekte/Paper_2/data/geomorph"

#Rounded Swiss extent plus a buffer of 4000m
extent <- c(2481300,1071200,2838000,1300000)

res_l1 <- 5
res_l2 <- 10
res_l3 <- 25

scale0 <- 20
scale1 <- 50
scale2 <- 200

perimeter_buffer <- 4000

epsg_lv95 <- 2056
crs_lv95 <- paste0("EPSG:",epsg_lv95)
crs_lv95 <- terra::crs(crs_lv95)

setwd(output_base_path)

tile_folder <- "ch_tiles"
tile_name_template <- paste0("dem_tile_",res_l1,"_ezg")
tile_name_template_l2 <- paste0("dem_tile_",res_l2,"_ezg")
tile_name_template_l3 <- paste0("dem_tile_",res_l3,"_ezg")

# Helper function to derive the path of a input raster based onthe tile folder, raster name and extension 
get_raster_path <- function(tile_output_folder,raster_name,raster_extension){
  rast_path <- file.path(".",tile_output_folder,paste0(raster_name,".",raster_extension))
  return(rast_path)
}

# Function returning the extend of a spatial object rounded to a specific resolution and modified by the specified buffer
get_rounded_extent <- function (sobj,resolution,buffer){
  extent_input <- terra::ext(sobj)
  
  xmin <- resolution*floor(xmin(extent_input)/resolution)-resolution-buffer
  ymin <- resolution*floor(ymin(extent_input)/resolution)-resolution-buffer
  xmax <- resolution*floor(xmax(extent_input)/resolution)+resolution+buffer
  ymax <- resolution*floor(ymax(extent_input)/resolution)+resolution+buffer
  
  # extent_output <- ext(xmin,xmax,ymin,ymax)
  extent_output <- terra::ext(c(xmin,xmax,ymin,ymax))
  
  return(extent_output)
}

# Load the saved models used for prediction
model_list <- list()
model_list["rf_model_hmdb"] <- "../../rf_model_hmdb.RData"
model_list["rf_model_be"] <- "../../rf_model_be.RData"

for(model_name in names(model_list)){
  model_path <- model_list[[model_name]]
  load(model_path)
}

# Prepare folder for prediction output
pred.out.dir <- file.path(tile_folder,perimeter_layer) 
dir.create(pred.out.dir)

# Determine perimeter and all catchment/tiles that overlap
perimeter <- st_read(perimeter_path, layer=perimeter_layer)
watersheds <- st_read(watershed_path, layer=watershed_layer)
watershed_mask <- rast(watershed_mask_path)

crs.lv95 <- st_crs("EPSG:2056")

perimeter.joined <- st_join(x = perimeter, y = watersheds)
tile_ids <- unique(perimeter.joined$TEZGNR150)
tile_ids

# Prepare inputs and predict raster for all tiles
for(tile_id in tile_ids) {

  # Preparing for preocessing
  print(tile_id)
  tile_name_l1 <- paste0(tile_name_template,tile_id)
  tile_name_l2 <- paste0(tile_name_template_l2,tile_id)
  tile_name_l3 <- paste0(tile_name_template_l3,tile_id)
  
  tile_output_folder_l1 <- file.path(tile_folder,tile_name_l1) 
  tile_output_folder_l2 <- file.path(tile_folder,tile_name_l2)   
  tile_output_folder_l3 <- file.path(tile_folder,tile_name_l3)   
  rasters <- list()
  
  smoothed_name <-""
  if(use_smooth_rasters) {
    smoothed_name <-"smoothed_"
  }
  
  pred_output_path_l1 <- file.path(pred.out.dir,paste0("pred_rf_model_hmdb_",smoothed_name,tile_name_template,tile_id,".tif"))
  if (file.exists(pred_output_path_l1)) {
    next
  }

  # Get input raster paths based on prepared covariate rasters. Not all are used for prediction.

  # Percipitation, 'prcp_60m_10a_0975',prcp_60m_30a_0975','prcp_60m_100a_0975','prcp_24h_30a_0975' ,'prcp_24h_100a_0975'
  rasters[paste0("prcp_60m_10a_0975")] <- get_raster_path(tile_output_folder_l3,paste0("prcp_60m_10a_0975",tile_name_l3),"tif")
  rasters[paste0("prcp_60m_30a_0975")] <- get_raster_path(tile_output_folder_l3,paste0("prcp_60m_30a_0975",tile_name_l3),"tif")
  rasters[paste0("prcp_60m_100a_0975")] <- get_raster_path(tile_output_folder_l3,paste0("prcp_60m_100a_0975",tile_name_l3),"tif")
  rasters[paste0("prcp_24h_30a_0975")] <- get_raster_path(tile_output_folder_l3,paste0("prcp_24h_30a_0975",tile_name_l3),"tif")
  rasters[paste0("prcp_24h_100a_0975")] <- get_raster_path(tile_output_folder_l3,paste0("prcp_24h_100a_0975",tile_name_l3),"tif")
  # Terrain Ruggedness Index (TRI),'tri_r5_10','tri_r5_5'
  rasters[paste0("tri_r5_",res_l1)] <- get_raster_path(tile_output_folder_l1,paste0("tri_r5_",tile_name_l1),"tif")
  rasters[paste0("tri_r5_",res_l2)] <- get_raster_path(tile_output_folder_l2,paste0("tri_r5_",tile_name_l2),"tif")
  #Topographic Position Indexes,'tpi_4km_25','tpi_500m_25'
  rasters[paste0("tpi_4km_",res_l3)] <- get_raster_path(tile_output_folder_l3,paste0("tpi_4km__",tile_name_l3),"tif")
  rasters[paste0("tpi_500m_",res_l3)] <- get_raster_path(tile_output_folder_l3,paste0("tpi_500m__",tile_name_l3),"tif")
  #VHM,'vhm_max_5'
  # rasters[paste0("vhm_max_",res_l3)] <- get_raster_path(tile_output_folder_l3,paste0("vhm_max_",tile_name_l3),"tif")
  # rasters[paste0("vhm_q3_",res_l3)] <- get_raster_path(tile_output_folder_l3,paste0("vhm_q3_",tile_name_l3),"tif")
  rasters[paste0("vhm_max_",res_l1)] <- get_raster_path(tile_output_folder_l1,paste0("vhm_max_",tile_name_l1),"tif")
  #Raw DEM height above sea level,'h_5'
  rasters[paste0("h_",res_l1)] <- get_raster_path(tile_output_folder_l1,tile_name_l1,"tif")
  # Vector Ruggedness Measure (VRM),'vrm_r5_25'
  rasters[paste0("vrm_r5_",res_l3)] <- get_raster_path(tile_output_folder_l3,paste0("vrm_r5_",tile_name_l3),"tif")
  # Multiresolution Index of the Ridge Top Flatness (MRRTF),'mrrtf_25'
  rasters[paste0("mrrtf_",res_l3)] <- get_raster_path(tile_output_folder_l3,paste0("mrrtf_",tile_name_l3),"tif")
  # Multiresolution Index of Valley Bottom Flatness (MRVBF),'mrvbf_10'
  rasters[paste0("mrvbf_",res_l2)] <- get_raster_path(tile_output_folder_l2,paste0("mrvbf_",tile_name_l2),"tif")
  rasters[paste0("mrvbf_",res_l3)] <- get_raster_path(tile_output_folder_l3,paste0("mrvbf_",tile_name_l3),"tif")
  # Slope Height,'slope_height_25'
  rasters[paste0("slope_height_",res_l3)] <- get_raster_path(tile_output_folder_l3,paste0("heigslope_",tile_name_l3),"tif")
  #Catchment slope,'catchment_slope_25'
  rasters[paste0("catchment_slope_",res_l3)] <- get_raster_path(tile_output_folder_l3,paste0("slope_twi_",tile_name_l3),"tif")
  # Valley Depth,'valley_depth_25'
  rasters[paste0("valley_depth_",res_l3)] <- get_raster_path(tile_output_folder_l3,paste0("heigval_",tile_name_l3),"tif")
  # Standardized Height,'standardized_height_25'
  rasters[paste0("standardized_height_",res_l3)] <- get_raster_path(tile_output_folder_l3,paste0("heigstand_",tile_name_l3),"tif")
  #Slope,'slope_10'
  rasters[paste0("slope_",res_l2)] <- get_raster_path(tile_output_folder_l2,paste0("slope_ind_",tile_name_l2),"tif")
 
  # Aspect Northness,'aspect_nness_10'
  rasters[paste0("aspect_nness_",res_l2)] <- get_raster_path(tile_output_folder_l2,paste0("asp_nness_",tile_name_l2),"tif")
  #Topographic Wetness Index SAGA,'twi_5','twi_25'
  rasters[paste0("twi_",res_l1)] <- get_raster_path(tile_output_folder_l1,paste0("TWIsaga_",tile_name_l1),"tif")
  #rasters[paste0("twi_",res_l3)] <- get_raster_path(tile_output_folder_l3,paste0("TWIsaga_",tile_name_l3),"tif")

  #General Curvature,'curvature_10'
  rasters[paste0("curvature_",res_l2)] <- get_raster_path(tile_output_folder_l2,paste0("curv_",tile_name_l2),"tif")
  #Profile Curvature,'curvature_profile_25'
  rasters[paste0("curvature_profile_",res_l3)] <- get_raster_path(tile_output_folder_l3,paste0("curvprov_",tile_name_l3),"tif")
  #Catchement Area in cells,'catchment_area_25'
  rasters[paste0("catchment_area_",res_l3)] <- get_raster_path(tile_output_folder_l3,paste0("catchment_",tile_name_l3),"tif")
  #Openness,'openness_pos_200_5'
  rasters[paste0("openness_pos_",scale2,"_",res_l1)] <- get_raster_path(tile_output_folder_l1,paste0("openness_pos_",tile_name_l1,"_",scale2),"tif")
  #Openness,'openness_neg_200_5'
  rasters[paste0("openness_neg_",scale2,"_",res_l1)] <- get_raster_path(tile_output_folder_l1,paste0("openness_neg_",tile_name_l1,"_",scale2),"tif")
  
  #Slide susceptibility deposits, 'susceptibility_unconsolidated_deposits_5f_0''susceptibility_unconsolidated_deposits_5f_1''susceptibility_unconsolidated_deposits_5f_2''susceptibility_unconsolidated_deposits_5f_3'
  rasters[paste0("susceptibility_unconsolidated_deposits_",res_l1,"f_0")] <- get_raster_path(tile_output_folder_l1,paste0("unconsolidated_deposits_susceptibility_",tile_name_l1,"f_0"),"tif")
  rasters[paste0("susceptibility_unconsolidated_deposits_",res_l1,"f_1")] <- get_raster_path(tile_output_folder_l1,paste0("unconsolidated_deposits_susceptibility_",tile_name_l1,"f_1"),"tif")
  rasters[paste0("susceptibility_unconsolidated_deposits_",res_l1,"f_2")] <- get_raster_path(tile_output_folder_l1,paste0("unconsolidated_deposits_susceptibility_",tile_name_l1,"f_1"),"tif")
  rasters[paste0("susceptibility_unconsolidated_deposits_",res_l1,"f_3")] <- get_raster_path(tile_output_folder_l1,paste0("unconsolidated_deposits_susceptibility_",tile_name_l1,"f_2"),"tif")
  
  #Slide susceptibility bedrock, 'susceptibility_gk500_5f_0''susceptibility_gk500_5f_1''susceptibility_gk500_5f_2''susceptibility_gk500_5f_3'
  
  rasters[paste0("susceptibility_gk500_",res_l1,"f_0")] <- get_raster_path(tile_output_folder_l1,paste0("gk500_susceptibility_",tile_name_l1,"f_0"),"tif")
  rasters[paste0("susceptibility_gk500_",res_l1,"f_1")] <- get_raster_path(tile_output_folder_l1,paste0("gk500_susceptibility_",tile_name_l1,"f_1"),"tif")
  rasters[paste0("susceptibility_gk500_",res_l1,"f_2")] <- get_raster_path(tile_output_folder_l1,paste0("gk500_susceptibility_",tile_name_l1,"f_1"),"tif")
  rasters[paste0("susceptibility_gk500_",res_l1,"f_3")] <- get_raster_path(tile_output_folder_l1,paste0("gk500_susceptibility_",tile_name_l1,"f_2"),"tif")
  
  #Gesteinsdichte,'rhob_m'
  rasters["rhob_m"] <- get_raster_path(tile_output_folder_l1,paste0("rhob_m__",tile_name_l1),"tiff")
  
  wtrshd<- watersheds[watersheds$TEZGNR150==tile_id,]

  tile_output_extent <- get_rounded_extent(wtrshd,5,5)
  tile_output_extent

  tmp.dir <- file.path(tile_output_folder_l1,"tmp_stack")
  dir.create(tmp.dir)

  # Create or set mask used for masking the tile to the catchment
  tile_mask_name <- "tile_mask.tif"
  tile_mask_path <- file.path(tmp.dir,tile_mask_name)
  
  if (!file.exists(tile_mask_path)){
    # Mask doesn't exist yet -> create mask 
    tile_mask <- watershed_mask==as.numeric(tile_id)
    tile_mask <- crop(tile_mask,tile_output_extent,filename=tile_mask_path,overwrite = TRUE, filetype = "GTiff", 
                      gdal=c("TILED=YES","COMPRESS=DEFLATE"),datatype="INT1U")
  } else {
    tile_mask <- rast(tile_mask_path)
  }

  # Greate raster stack with inputs for prediction based on covariate rasters. All inputs are resampled to 5m resolution.
  stack_rasters <- list()
  for(n in names(rasters)){
    raster_path <- rasters[[n]]
    t_ras <- rast(raster_path)
    t_res <- terra::res(t_ras)
    t_ext <- ext(t_ras)
    
    t_base_name <- basename(raster_path)
    print(t_base_name)
    
    layer_path <- file.path(tmp.dir,t_base_name)
    layer_path_smooth <- gsub(".tif","_smooth.tif",layer_path)
    # Create 5m resolution version of the raster if it doesn't exist yet
    if (!file.exists(layer_path)){
      gdalwarp(srcfile=normalizePath(raster_path),dstfile=normalizePath(layer_path),t_srs=crs_lv95, overwrite=TRUE,co=c("TILED=YES","COMPRESS=DEFLATE"), r = "near", tr=c(res_l1,res_l1), te = c(tile_output_extent[1],tile_output_extent[3],tile_output_extent[2],tile_output_extent[4])
               ,dryrun = FALSE)
    }
    
    # Create raster smoothed with a mean function and a 11 pixel window for inputs with >5m resolution. 
    if((startsWith(n,"prcp") | endsWith(n,"_25") | endsWith(n,"_10"))&(!file.exists(layer_path_smooth))) {
      w.width <- 11
      # print(layer_path_smooth)
      t_rast <- rast(layer_path)
      t_rast_smooth <- focal(t_rast, w = matrix(1, ncol = w.width, nrow = w.width), fun = mean, expand=TRUE,
                             overwrite = TRUE, filetype = "GTiff", 
                             gdal=c("TILED=YES","COMPRESS=DEFLATE"),
                             filename = layer_path_smooth)
    }
    
    # Switch path of input to smoothed version where applicable if the global option for smoothing is set.
    if (use_smooth_rasters) {
      if(startsWith(n,"prcp") | endsWith(n,"_25") | endsWith(n,"_10")) {
        layer_path <- layer_path_smooth
      }
    }
    
    # Add input path to input list for raster stack
    stack_rasters[n] <- layer_path
  }
  
  # Create actual raster stack
  pred_in_stk  <-  rast(unlist(stack_rasters, use.names = FALSE))

  names(pred_in_stk) <- c(names(stack_rasters))
  crs(pred_in_stk) <- crs_lv95

  # Fill NA values in inputs with default values
  pred_in_stk$vhm_max_5[is.na(pred_in_stk$vhm_max_5)]<-0
  pred_in_stk$prcp_60m_10a_0975[is.na(pred_in_stk$prcp_60m_10a_0975)]<-0
  pred_in_stk$prcp_60m_30a_0975[is.na(pred_in_stk$prcp_60m_30a_0975)]<-0
  pred_in_stk$prcp_60m_100a_0975[is.na(pred_in_stk$prcp_60m_100a_0975)]<-0
  pred_in_stk$prcp_24h_30a_0975[is.na(pred_in_stk$prcp_24h_30a_0975)]<-0
  pred_in_stk$prcp_24h_100a_0975[is.na(pred_in_stk$prcp_24h_100a_0975)]<-0
  pred_in_stk$aspect_nness_10[is.na(pred_in_stk$aspect_nness_10)]<-0
  #pred_in_stk$rhob_m[is.na(pred_in_stk$rhob_m)]<-0
  
  # Predict raster for model trained with HMDB data and mask to catchment
  pred <- terra::predict(object=pred_in_stk,model=rf_model_hmdb, na.action = "na.pass")
  pred <- tile_mask*pred

  pred_output_path_l1 <- file.path(pred.out.dir,paste0("pred_rf_model_hmdb_",smoothed_name,tile_name_template,tile_id,".tif"))
  print(pred_output_path_l1)
  writeRaster(pred, pred_output_path_l1, filetype = "GTiff", overwrite=TRUE,
              gdal=c("TILED=YES","COMPRESS=DEFLATE"))

  # Predict raster for model trained with KtBE data  and mask to catchment
  pred <- terra::predict(object=pred_in_stk,model=rf_model_be, na.action = "na.pass")
  pred <- tile_mask*pred

  pred_output_path_l1 <- file.path(pred.out.dir,paste0("pred_rf_model_be_",smoothed_name,tile_name_template,tile_id,".tif"))
  print(pred_output_path_l1)
  writeRaster(pred, pred_output_path_l1, filetype = "GTiff", overwrite=TRUE,
              gdal=c("TILED=YES","COMPRESS=DEFLATE"))
  
  # Remove temp directory with prediction inputs
  unlink(tmp.dir, recursive=TRUE)
  
}

# Merge the predicted tiles into one output raster
vrt_output_path_l1 <- file.path(pred.out.dir,"pred_rf_model_hmdb.vrt")
tif_output_path_l1 <- file.path(pred.out.dir,"pred_rf_model_hmdb.tif")
gdalUtilities::gdalbuildvrt(gdalfile = Sys.glob(file.path(getwd(),pred.out.dir, "pred_rf_model_hmdb_*.tif")), output.vrt = vrt_output_path_l1, srcnodata=0)
gdalUtilities::gdal_translate(src_dataset = vrt_output_path_l1,dst_dataset = tif_output_path_l1,of="GTiff",co=c("TILED=YES","COMPRESS=DEFLATE","BIGTIFF=YES"))

vrt_output_path_l1 <- file.path(pred.out.dir,"pred_rf_model_be.vrt")
tif_output_path_l1 <- file.path(pred.out.dir,"pred_rf_model_be.tif")
gdalUtilities::gdalbuildvrt(gdalfile = Sys.glob(file.path(getwd(),pred.out.dir, "pred_rf_model_be_*.tif")), output.vrt = vrt_output_path_l1, srcnodata=0)
gdalUtilities::gdal_translate(src_dataset = vrt_output_path_l1,dst_dataset = tif_output_path_l1,of="GTiff",co=c("TILED=YES","COMPRESS=DEFLATE","BIGTIFF=YES"))




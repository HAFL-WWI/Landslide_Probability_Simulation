#######################################################################
# Script for generating with helper functions.
#
# Author: Christoph Schaller, BFH-HAFL, April 2025
######################################################################

import os
import sys 
import glob

import rasterio
from rasterio import windows
import rasterio.mask
import rasterio.merge 

from osgeo import ogr, gdal

import numpy as np


import logging
from multiprocessing import current_process

# Function to ensure that a directory exists
def ensure_dir (path):
    if not os.path.isdir(path):
        os.mkdir(path)

# Function that returns a list of immediate subdirectories
def get_subdirs (path):
    dirs = []
    for f in os.listdir(path):
        full_path = os.path.join(path,f)
        if(os.path.isdir(full_path)):
          dirs.append([f,full_path])

    return dirs

# Function that deletes an existing raster using gdal
def delete_raster(raster):
    data = gdal.Open(raster, gdal.GA_ReadOnly)
    driver = data.GetDriver()
    data = None
    if os.path.exists(raster):
        driver.Delete(raster)

# Function that deletes all files with a specific extension in a path
def delete_files_by_extension (path,extension):
    for filepath in glob.iglob(os.path.join(path, '*.'+extension)):
        os.remove(filepath)

#Function to copy a GeoTIFF raster
def copy_raster_tiff(in_raster, out_raster):
    driver = gdal.GetDriverByName('GTiff')
    in_ds = gdal.Open(in_raster)
    out_ds = driver.CreateCopy(out_raster, in_ds, 0)
    in_ds = None
    out_ds = None

# Function for cropping a raster based on a set of vector geometries
def crop_image(input_path, output_path, mask, dtype=None):
    with rasterio.open(input_path,mode='r') as input_image:
        no_data = input_image.nodata
        out_image, out_transform = rasterio.mask.mask(input_image, mask,
                                                        all_touched=True, nodata=no_data, crop=True)


        out_meta = input_image.meta.copy()
        out_meta.update({#"driver": "GTiff",
                         "height": out_image.shape[1],
                         "width": out_image.shape[2],
                         "transform": out_transform,
                         "compress":"lzw"})
        if dtype:
            out_meta.update({"dtype":dtype})
        if no_data:
            out_meta.update({"nodata":no_data})

        with rasterio.open(output_path, "w", **out_meta) as dest:
            dest.write(out_image)


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


# Function that deletes an existing raster using gdal
def delete_raster(raster):
    if os.path.exists(raster):
        data = gdal.Open(raster, gdal.GA_ReadOnly)
        driver = data.GetDriver()
        data = None

        driver.Delete(raster)

# Function that returns the extent or a shapefile
def get_shp_extent (shapefile_path):
    inDriver = ogr.GetDriverByName("ESRI Shapefile")
    inDataSource = inDriver.Open(shapefile_path, 0)
    inLayer = inDataSource.GetLayer()
    extent = inLayer.GetExtent()
    inDataSource = None
    return extent

# Function that returns the extent of a geopackage vector layer
def get_gpkg_extent (gpkg_path,layer_name):
    inDriver = ogr.GetDriverByName("GPKG")
    inDataSource = inDriver.Open(gpkg_path, 0)
    inLayer = inDataSource.GetLayerByName(layer_name)
    extent = inLayer.GetExtent()
    inDataSource = None
    return extent

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

# Function that checks if a layer exists in a geopackage
def gpkg_layer_exists(gpkg_path, layer_name):
    if os.path.exists(gpkg_path):
        con = ogr.Open(gpkg_path)
        for l in con:
            if l.GetName()==layer_name:
                return True
        return False
    else:
        return False

# Function for adding two rasters of the same size 
def add_rasters (a_raster, b_raster, out_raster):
        # Open the datasets
        dsIn_a = gdal.Open(a_raster, gdal.GA_ReadOnly)
        b1_a = dsIn_a.GetRasterBand(1)
        vNA_a = b1_a.GetNoDataValue()
        dataIn_a = np.array(b1_a.ReadAsArray())

        dsIn_b = gdal.Open(b_raster, gdal.GA_ReadOnly)
        b1_b = dsIn_b.GetRasterBand(1)
        vNA_b = b1_b.GetNoDataValue()
        dataIn_b = np.array(b1_b.ReadAsArray())

        # Add rasters
        dataIn_b[dataIn_b == vNA_b] = 0
        dataOut = np.copy(dataIn_a)
        dataOut = dataOut+dataIn_b

        # Write the out file
        delete_raster(out_raster)
        driver = gdal.GetDriverByName("GTiff")
        dst_options = ['COMPRESS=LZW']
        dsOut = driver.Create(out_raster, dsIn_a.RasterXSize, dsIn_a.RasterYSize, 1, b1_a.DataType, dst_options)
        dsOut.GetRasterBand(1).SetNoDataValue(vNA_a)
        gdal.CopyDatasetInfo(dsIn_a, dsOut)
        bandOut = dsOut.GetRasterBand(1)
        gdal.BandWriteArray(bandOut, dataOut)


        dsIn_a = None
        dsIn_b = None
        dsOut = None

# Function that removes the no data value from all bands in a raster
def unset_nodata (raster):
    ds = gdal.Open(raster, gdal.GA_Update)
    for i in range(ds.RasterCount):
        ds.GetRasterBand(i + 1).DeleteNoDataValue()
    ds = None


# Function that sets the no data value on all bands in a raster toe the specified value
def set_nodata (raster,nodata):
    ds = gdal.Open(raster, gdal.GA_Update)
    for i in range(ds.RasterCount):
        ds.GetRasterBand(i + 1).SetNoDataValue(nodata)

# Function that first unsets the no data value in a raster and then sets a new value
def reset_nodata (raster,nodata):
    unset_nodata (raster)
    set_nodata (raster,nodata)


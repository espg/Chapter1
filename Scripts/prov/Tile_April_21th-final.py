from osgeo import gdal
import georasters as gr
import hashlib
import numpy as np
from pylab import r_
import osr
import os
import glob
from tqdm import tqdm

# Better version of this function exists...
# ...one that isn't hard coded to Polar Stereo North
def resample_tif(tif_file, pixel_spacing=0.5):
    source = gdal.Open(tif_file)
    source.GetRasterBand(1).SetNoDataValue(-32767)
    wgs84 = osr.SpatialReference() # slopppy 
    wgs84.ImportFromEPSG(3413)     # only for polar stereographic!
    # Get the Geotransform vector
    geo_t = source.GetGeoTransform ()
    x_size = source.RasterXSize # Raster xsize
    y_size = source.RasterYSize # Raster ysize
    # Work out the boundaries of the new dataset in the target projection
    ulx, uly = geo_t[0], geo_t[3]
    lrx = geo_t[0] + geo_t[1]*x_size
    lry = geo_t[3] + geo_t[5]*y_size
    # Now, we create an in-memory raster
    mem_drv = gdal.GetDriverByName('MEM')
    # The size of the raster is given the new projection and pixel spacing
    # Using the values we calculated above. Also, setting it to store one band
    # and to use Float32 data type.
    dest = mem_drv.Create('', int((lrx - ulx)/pixel_spacing),             int((uly - lry)/pixel_spacing), 1, gdal.GDT_Float32)
    # Calculate the new geotransform
    new_geo = ( ulx, pixel_spacing, geo_t[2],                 uly, geo_t[4], -pixel_spacing )
    # Set the geotransform
    dest.SetGeoTransform(new_geo)
    dest.SetProjection(source.GetProjection())
    # Perform the projection/resampling 
    res = gdal.ReprojectImage(source, dest, None, None,
                gdal.GRA_Bilinear)
    graster = gr.GeoRaster(dest.GetRasterBand(1).ReadAsArray(),
                           dest.GetGeoTransform(),0,
                           wgs84,gdal.GDT_Float32)
    return graster


def grid_write(graster):
    rows, cols = graster.shape
    for i in r_[0:rows:200]:
        for j in r_[0:cols:200]:
            temp = graster[i:i+200,j:j+200]
            if np.sum(temp.raster.mask) > 0:
                pass
            else:
                if temp.shape == (200,200):
                    temp = temp - temp.mean()
                    temp.to_tiff(hashlib.sha256(temp.raster).hexdigest())


input_files = "/Users/sgrigsby/Desktop/April_21st_DEMs/IODMS3*.tif"

os.chdir("/Users/sgrigsby/Desktop/April_21st_tile_output")

for filename in tqdm(glob.glob(input_files)):
    grid_write(resample_tif(filename))


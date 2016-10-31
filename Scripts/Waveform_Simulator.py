from numba import jit, autojit
from osgeo import gdal
import numpy as np
import pandas as pd
import os
from tqdm import tqdm


def waveformC(raster):
    # Define constants
    cumpowr=np.zeros((544,1))   # Signal length
    # Define spatial weighting 
    x, y = np.r_[-50:50:0.5], np.r_[-50:50:0.5]
    X, Y = np.meshgrid(x,y)
    w=np.exp(-(X**2 + Y**2)/2/17.5**2)
    # Conversion to nanoseconds... 'D' is nanoseconds
    Z = raster[:,:] - np.max(raster[:,:])
    D = (abs(Z) / 0.15) + 10    # Start all recordings after 10 ns
    for i in range(200):
        for j in range(200):
            ibin=np.int(D[i,j]) #finding ibin
            #add power within pixel (i,j)-which is w(i,j) to the returned
            #power in the bin number ibin
            if ibin < 544:
                cumpowr[ibin]=cumpowr[ibin]+w[i,j]
    return cumpowr

wave_numba = autojit(waveformC)


folder = "/Users/sgrigsby/Desktop/April_21st_tile_output/"


idx = []
for i, filename in enumerate(os.listdir(folder)):
    fn, ext = os.path.splitext(filename)
    # source = gdal.Open(folder + filename)
    # (source.RasterYSize & source.RasterXSize) == 200
    idx.append(fn)


#April_20th_wf_no_conv = pd.DataFrame(index=idx,columns=np.r_[0:544:1])
space = np.zeros((len(idx),544),dtype=float)

for i, filename in enumerate(tqdm(os.listdir(folder))):
    #fn, ext = os.path.splitext(filename)
    source = gdal.Open(folder + filename)
    if (source.RasterYSize & source.RasterXSize) == 200:
        space[i,:] = wave_numba(source.ReadAsArray()).ravel()
    del(source)

    #April_20th_wf_no_conv.iloc[i] = wave_numba(source.ReadAsArray()).ravel()


# In[ ]:
April_21st_wf_no_conv = pd.DataFrame(space,index=idx,columns=np.r_[0:544:1])
# store = pd.HDFStore('/Users/sgrigsby/Desktop/April_20th_wfs.h5','a')
store = pd.HDFStore('/Users/sgrigsby/Desktop/April_21st_wfs.h5','a')
store['raw'] = April_21st_wf_no_conv
store.close()



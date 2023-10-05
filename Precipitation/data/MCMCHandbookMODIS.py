#!/usr/bin/env python
# coding: utf-8

import os
import warnings

import matplotlib.pyplot as plt
import numpy.ma as ma
import xarray as xr
import rasterio as rio
import rioxarray as rxr
from shapely.geometry import mapping, box
import geopandas as gpd
import earthpy as et
import earthpy.spatial as es
import earthpy.plot as ep
# https://www.earthdatascience.org/courses/use-data-open-source-python/hierarchical-data-formats-hdf/open-MODIS-hdf4-files-python/


from pyhdf.SD  import *
fileName = 'MOD05_L2.A2023173.1320.061.2023174011800.hdf'
hdf = SD(fileName, SDC.READ)

# Read dataset.
DATAFIELD_NAME='Water_Vapor_Near_Infrared'
data3D = hdf.select(DATAFIELD_NAME)
data = data3D[:,:]

# Read dataset.
DATAFIELD_NAME='Water_Vapor_Infrared'
data3DCoarse = hdf.select(DATAFIELD_NAME)
dataCoarse = data3DCoarse[:,:]

# Read geolocation dataset.
lat = hdf.select('Latitude')
latitude = lat[:,:]
print(latitude.shape)

lon = hdf.select('Longitude')
longitude = lon[:,:]

from numpy import savetxt

savetxt('latitude.csv', latitude, delimiter=',')
savetxt('longitude.csv', longitude, delimiter=',')
savetxt('data.csv', data, delimiter=',')
savetxt('data5km.csv', dataCoarse, delimiter=',')



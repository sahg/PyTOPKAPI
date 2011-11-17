#!/usr/bin/python
import os
import sys

import numpy as np
import grass.script as grass

catch_shp = sys.argv[1]
catch_id = sys.argv[2]
buffer_range = sys.argv[3]

# make sure the region is set correctly
grass.run_command('g.region', rast='srtm_dem')

# create initial catchment mask
grass.run_command('v.buffer',
                  flags = '-o',
                  input = catch_shp,
                  output = '%s_buffer_py' % catch_id,
                  distance='%s' % buffer_range)

grass.run_command('v.to.rast',
                  flags = '-o',
                  input = '%s_buffer_py' % catch_id,
                  out = '%s_mask_py' % catch_id,
                  use = 'val')

# use mask to extract catchment from DEM
grass.run_command('r.mask', input='%s_mask_py' % catch_id)
grass.mapcalc('%s_dem=srtm_dem' % catch_id, overwrite=True)
grass.run_command('r.mask', flags='r')

# obtain parameters from DEM
grass.run_command('r.watershed',
                  flags = '-o',
                  elevation = '%s_dem' % catch_id,
                  drainage = '%s_dir' % catch_id,
                  stream = '%s_str' % catch_id,
                  threshold=250)

# Thin the stream raster and convert to vector for visualization
# purposes
grass.run_command('r.thin',
                  flags = '-o',
                  input = '%s_str' % catch_id,
                  out = '%s_str_thin' % catch_id)

grass.run_command('r.to.vect',
                  flags = '-o',
                  input = '%s_str_thin' % catch_id,
                  out = '%s_streams' % catch_id,
                  feature = 'line')

# Find outlet point and upstream cells
grass.run_command('v.patch',
                  flags = '-o',
                  input = '%s_streams,%s' % (catch_id, catch_shp),
                  out = '%s_streams_catch' % catch_id)

grass.run_command('v.clean',
                  flags = '-o',
                  input = '%s_streams_catch' % catch_id,
                  out = '%s_streams_catch_clean' % catch_id,
                  tool = 'break',
                  error = '%s_catch_outlet' % catch_id)

grass.run_command('v.out.ascii',
                  flags = '-o',
                  input = '%s_catch_outlet' % catch_id,
                  out = '%s_catch_outlet.txt' % catch_id,
                  format = 'point',
                  fs = 'space')

easting, northing = np.loadtxt('%s_catch_outlet.txt' % catch_id)
os.remove('%s_catch_outlet.txt' % catch_id)

grass.run_command('r.water.outlet',
                  flags = '-o',
                  drain = '%s_dir' % catch_id,
                  east = easting,
                  north = northing,
                  basin = '%s_catch_mask' % catch_id)

grass.run_command('r.null',
                  map = '%s_catch_mask' % catch_id,
                  setnull = '0')

# crop region to catchment mask
grass.run_command('g.region',
                  zoom = '%s_catch_mask' % catch_id,
                  align = '%s_catch_mask' % catch_id)

# use new mask to extract catchment information from regional rasters
# or catchment buffer
grass.run_command('r.mask', input='%s_catch_mask' % catch_id)
grass.mapcalc('%s_dem=%s_dem' % (catch_id, catch_id), overwrite=True)
grass.mapcalc('%s_slope=srtm_slope' % catch_id, overwrite=True)
grass.mapcalc('%s_n_overland=manning_overland' % catch_id, overwrite=True)
grass.mapcalc('%s_soil_depth=soil_depth' % catch_id, overwrite=True)
grass.mapcalc('%s_sat_moisture_content=sat_moisture_content' \
              % catch_id, overwrite=True)
grass.mapcalc('%s_residual_moisture_content=residual_moisture_content' \
              % catch_id, overwrite=True)
grass.mapcalc('%s_hydraulic_conductivity=hydraulic_conductivity' \
              % catch_id, overwrite=True)
grass.mapcalc('%s_pore_size=pore_size' % catch_id, overwrite=True)
grass.mapcalc('%s_bubbling_pressure=bubbling_pressure' \
              % catch_id, overwrite=True)
grass.mapcalc('%s_dir=%s_dir' % (catch_id, catch_id), overwrite=True)
grass.mapcalc('%s_str_thin=%s_str_thin' % (catch_id, catch_id), overwrite=True)
grass.run_command('r.mask', flags='r')

# output to GTiff for further processing
grass.run_command('r.out.gdal',
                  type='Int16',
                  input='%s_dir' % catch_id,
                  format='GTiff',
                  output='%s-flow-dir.tif' % catch_id)

grass.run_command('r.out.gdal',
                  input='%s_catch_mask' % catch_id,
                  format='GTiff',
                  output='%s-mask.tif' % catch_id)

grass.run_command('r.out.gdal',
                  input='%s_dem' % catch_id,
                  format='GTiff',
                  output='%s-dem.tif' % catch_id)

grass.run_command('r.out.gdal',
                  input='%s_str_thin' % catch_id,
                  format='GTiff',
                  output='%s-channel-network.tif' % catch_id)

grass.run_command('r.out.gdal',
                  input='%s_n_overland' % catch_id,
                  format='GTiff',
                  output='%s-manning-overland.tif' % catch_id)

grass.run_command('r.out.gdal',
                  input='%s_soil_depth' % catch_id,
                  format='GTiff',
                  output='%s-soil-depth.tif' % catch_id)

grass.run_command('r.out.gdal',
                  input='%s_sat_moisture_content' % catch_id,
                  format='GTiff',
                  output='%s-sat-moisture-content.tif' % catch_id)

grass.run_command('r.out.gdal',
                  input='%s_residual_moisture_content' % catch_id,
                  format='GTiff',
                  output='%s-residual-moisture-content.tif' % catch_id)

grass.run_command('r.out.gdal',
                  input='%s_hydraulic_conductivity' % catch_id,
                  format='GTiff',
                  output='%s-hydraulic-conductivity.tif' % catch_id)

grass.run_command('r.out.gdal',
                  input='%s_pore_size' % catch_id,
                  format='GTiff',
                  output='%s-pore-size.tif' % catch_id)

grass.run_command('r.out.gdal',
                  input='%s_bubbling_pressure' % catch_id,
                  format='GTiff',
                  output='%s-bubbling-pressure.tif' % catch_id)

grass.run_command('r.out.gdal',
                  input='%s_slope' % catch_id,
                  format='GTiff',
                  output='%s-slope.tif' % catch_id)

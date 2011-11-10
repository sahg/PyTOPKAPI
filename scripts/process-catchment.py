import numpy as np
import grass.script as grass

buffer_range = 5000
catch_id = 'B60'

# make sure the region is set correctly
grass.run_command('g.region', rast='srtm-dem')

# create catchment mask
grass.run_command('v.extract',
                  flags = '-o',
                  input = 'sa_quaternaries',
                  output = '%s_py' % catch_id,
                  where="TERTIARY = '%s'" % catch_id)

grass.run_command('v.dissolve',
                  flags = '-o',
                  input = '%s_py' % catch_id,
                  output = '%s_dissolve_py' % catch_id,
                  column='TERTIARY')

grass.run_command('v.buffer',
                  flags = '-o',
                  input = '%s_dissolve_py' % catch_id,
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
grass.run_command('r.slope.aspect',
                  flags = '-o',
                  elevation = '%s_dem' % catch_id,
                  slope = '%s_slope' % catch_id)

grass.run_command('r.watershed',
                  flags = '-o',
                  elevation = '%s_dem' % catch_id,
                  accumulation = '%s_accum' % catch_id,
                  drainage = '%s_dir' % catch_id,
                  stream = '%s_str' % catch_id,
                  threshold=250)

# the stream raster usually requires thinning
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
                  input = '%s_streams,%s_dissolve_py' % (catch_id, catch_id),
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

# use new mask to extract catchment information
grass.run_command('r.mask', input='%s_catch_mask' % catch_id)
grass.mapcalc('%s_slope=srtm_slopes' % catch_id, overwrite=True)
grass.mapcalc('%s_dir=%s_dir' % (catch_id, catch_id), overwrite=True)
grass.run_command('r.mask', flags='r')

# output to GTiff for further processing
grass.run_command('r.out.gdal',
                  type='Int16',
                  input='%s_dir' % catch_id,
                  format='GTiff',
                  output='lieb-flow-dir.tif')

grass.run_command('r.out.gdal',
                  input='%s_catch_mask' % catch_id,
                  format='GTiff',
                  output='lieb-mask.tif')

grass.run_command('r.out.gdal',
                  input='%s_slope' % catch_id,
                  format='GTiff',
                  output='lieb-slope.tif')

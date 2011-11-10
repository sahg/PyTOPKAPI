import numpy as np
import grass.script as grass

buffer_range = 5000
tertiary = 'B60'

# make sure the region is set correctly
grass.run_command('g.region', rast='srtm-dem')

# create catchment mask
grass.run_command('v.extract',
                  flags = '-o',
                  input = 'sa_quaternaries',
                  output = '%s_py' % tertiary,
                  where="TERTIARY = '%s'" % tertiary)

grass.run_command('v.dissolve',
                  flags = '-o',
                  input = '%s_py' % tertiary,
                  output = '%s_dissolve_py' % tertiary,
                  column='TERTIARY')

grass.run_command('v.buffer',
                  flags = '-o',
                  input = '%s_dissolve_py' % tertiary,
                  output = '%s_buffer_py' % tertiary,
                  distance='%s' % buffer_range)

grass.run_command('v.to.rast',
                  flags = '-o',
                  input = '%s_buffer_py' % tertiary,
                  out = '%s_mask_py' % tertiary,
                  use = 'val')

# use mask to extract catchment from DEM
grass.run_command('r.mask', input='%s_mask_py' % tertiary)
grass.mapcalc('%s_dem=srtm_dem' % tertiary, overwrite=True)
grass.run_command('r.mask', flags='r')

# obtain parameters from DEM
grass.run_command('r.slope.aspect',
                  flags = '-o',
                  elevation = '%s_dem' % tertiary,
                  slope = '%s_slope' % tertiary)

grass.run_command('r.watershed',
                  flags = '-o',
                  elevation = '%s_dem' % tertiary,
                  accumulation = '%s_accum' % tertiary,
                  drainage = '%s_dir' % tertiary,
                  stream = '%s_str' % tertiary,
                  threshold=250)

# the stream raster usually requires thinning
grass.run_command('r.thin',
                  flags = '-o',
                  input = '%s_str' % tertiary,
                  out = '%s_str_thin' % tertiary)

grass.run_command('r.to.vect',
                  flags = '-o',
                  input = '%s_str_thin' % tertiary,
                  out = '%s_streams' % tertiary,
                  feature = 'line')

# Find outlet point and upstream cells
grass.run_command('v.patch',
                  flags = '-o',
                  input = '%s_streams,%s_dissolve_py' % (tertiary, tertiary),
                  out = '%s_streams_catch' % tertiary)

grass.run_command('v.clean',
                  flags = '-o',
                  input = '%s_streams_catch' % tertiary,
                  out = '%s_streams_catch_clean' % tertiary,
                  tool = 'break',
                  error = '%s_catch_outlet' % tertiary)

grass.run_command('v.out.ascii',
                  flags = '-o',
                  input = '%s_catch_outlet' % tertiary,
                  out = '%s_catch_outlet.txt' % tertiary,
                  format = 'point',
                  fs = 'space')

easting, northing = np.loadtxt('%s_catch_outlet.txt' % tertiary)

grass.run_command('r.water.outlet',
                  flags = '-o',
                  drain = '%s_dir' % tertiary,
                  east = easting,
                  north = northing,
                  basin = '%s_catch_mask' % tertiary)

grass.run_command('r.null',
                  map = '%s_catch_mask' % tertiary,
                  setnull = '0')

# crop region to catchment mask
grass.run_command('g.region',
                  zoom = '%s_catch_mask' % tertiary,
                  align = '%s_catch_mask' % tertiary)

# use new mask to extract catchment information
grass.run_command('r.mask', input='%s_catch_mask' % tertiary)
grass.mapcalc('%s_slope=srtm_slopes' % tertiary, overwrite=True)
grass.mapcalc('%s_dir=%s_dir' % (tertiary, tertiary), overwrite=True)
grass.run_command('r.mask', flags='r')

# output to GTiff for further processing
grass.run_command('r.out.gdal',
                  type='Int16',
                  input='%s_dir' % tertiary,
                  format='GTiff',
                  output='lieb-flow-dir.tif')

grass.run_command('r.out.gdal',
                  input='%s_catch_mask' % tertiary,
                  format='GTiff',
                  output='lieb-mask.tif')

grass.run_command('r.out.gdal',
                  input='%s_slope' % tertiary,
                  format='GTiff',
                  output='lieb-slope.tif')

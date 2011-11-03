import grass.script as grass

accum_thresh = 30

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
                  distance='5000.0')

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
                  input = '%s_str_thin' % tertiary,
                  out = '%s_streams' % tertiary,
                  feature = 'line')

grass.mapcalc('%s_accum_thresh=abs(%s_accum) >= %s' % (tertiary, tertiary, accum_thresh),
              overwrite=True)

grass.mapcalc('%s_accum_abs=abs(%s_accum)' % (tertiary, tertiary),
              overwrite=True)

# crop region to mask
grass.run_command('g.region',
                  zoom='%s_mask_py' % tertiary,
                  align='%s_mask_py' % tertiary)

# output to GTiff for further processing
grass.run_command('r.out.gdal',
                  type='Int16',
                  input='%s_dir' % tertiary,
                  format='GTiff',
                  output='lieb-flow-dir.tif')

grass.run_command('r.out.gdal',
                  input='%s_mask_py' % tertiary,
                  format='GTiff',
                  output='lieb-mask.tif')

grass.run_command('r.out.gdal',
                  input='%s_slope' % tertiary,
                  format='GTiff',
                  output='lieb-slope.tif')

grass.run_command('r.out.gdal',
                  input='%s_accum_abs' % tertiary,
                  format='GTiff',
                  output='lieb-accum-abs.tif')

grass.run_command('r.out.gdal',
                  input='%s_accum_thresh' % tertiary,
                  format='GTiff',
                  output='lieb-accum-thresh.tif')

import grass.script as grass

accum_thresh = 30

# make sure the region is set correctly
grass.run_command('g.region', rast='srtm-dem')

# create catchment mask
grass.run_command('v.extract',
                  flags = '-o',
                  input = 'sa_quaternaries',
                  output = 'liebenbergsvlei_py',
                  where="TERTIARY = 'C83'")

grass.run_command('v.dissolve',
                  flags = '-o',
                  input = 'liebenbergsvlei_py',
                  output = 'lieb_dissolve_py',
                  column='TERTIARY')

grass.run_command('v.buffer',
                  flags = '-o',
                  input = 'lieb_dissolve_py',
                  output = 'lieb_buffer_py',
                  distance='5000.0')

grass.run_command('v.to.rast',
                  flags = '-o',
                  input = 'lieb_buffer_py',
                  out = 'lieb_mask_py',
                  use = 'val')

# use mask to extract catchment from DEM
grass.run_command('r.mask', input='lieb_mask_py')
grass.mapcalc('lieb_dem=srtm_dem', overwrite=True)
grass.run_command('r.mask', flags='r')

# obtain parameters from DEM
grass.run_command('r.slope.aspect',
                  flags = '-o',
                  elevation = 'lieb_dem',
                  slope = 'lieb_slope')

grass.run_command('r.watershed',
                  flags = '-o',
                  elevation = 'lieb_dem',
                  accumulation = 'lieb_accum',
                  drainage = 'lieb_dir',
                  threshold=250)

grass.mapcalc('lieb_accum_thresh=abs(lieb_accum) >= %s' % accum_thresh,
              overwrite=True)

grass.mapcalc('lieb_accum_abs=abs(lieb_accum)', overwrite=True)

# crop region to mask
grass.run_command('g.region', zoom='lieb_mask_py', align='lieb_mask_py')

# output to GTiff for further processing
grass.run_command('r.out.gdal',
                  type='Int16',
                  input='lieb_dir',
                  format='GTiff',
                  output='lieb-flow-dir.tif')

grass.run_command('r.out.gdal',
                  input='lieb_mask',
                  format='GTiff',
                  output='lieb-mask.tif')

grass.run_command('r.out.gdal',
                  input='lieb_slope',
                  format='GTiff',
                  output='lieb-slope.tif')

grass.run_command('r.out.gdal',
                  input='lieb_accum_abs',
                  format='GTiff',
                  output='lieb-accum-abs.tif')

grass.run_command('r.out.gdal',
                  input='lieb_accum_thresh',
                  format='GTiff',
                  output='lieb-accum-thresh.tif')

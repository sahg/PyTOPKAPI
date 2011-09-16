import grass.script as grass

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

grass.run_command('v.to.rast',
                  flags = '-o',
                  input = 'lieb_dissolve_py',
                  out = 'lieb_mask_py',
                  use = 'val')

# use mask to extract portion of DEM
grass.run_command('r.mask', input='lieb_mask_py')
grass.mapcalc('lieb_dem=srtm_dem', overwrite=True)
grass.run_command('r.mask', flags='r')

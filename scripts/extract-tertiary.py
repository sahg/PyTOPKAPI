import numpy as np
import grass.script as grass

catch_id = 'C83'

# create catchment mask
grass.run_command('v.extract',
                  flags = '-o',
                  input = 'wr2005_catchments',
                  output = '%s_py' % catch_id,
                  where="TERTIARY = '%s'" % catch_id)

grass.run_command('v.dissolve',
                  flags = '-o',
                  input = '%s_py' % catch_id,
                  output = '%s_dissolve_py' % catch_id,
                  column='TERTIARY')

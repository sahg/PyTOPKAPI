# This test works on a Windows XP system running:
# Python 2.5.1
# numpy.__version__ '1.0.3'
# scipy.__version__ '0.5.2'
# tables.__version__ '2.0'

import TOPKAPI # import the TOPKAPI package

# create a cell parameter file using the config in create_file.ini
from TOPKAPI.parameter_utils import create_file
create_file.run('create_the_parameter_files/create_file.ini')
# Modify the created file
from TOPKAPI.parameter_utils import modify_file
# create a new cell parameter file for a subcatchment
modify_file.subcatch('create_the_parameter_files/subcatch.ini')
# create a new cell parameter file with modified parameter values
modify_file.new_param('create_the_parameter_files/new_param.ini')
# create a new cell parameter file with external flows connected to the channel network
modify_file.connect_external_flow('create_the_parameter_files/connect_external_flow.ini')
# create a new cell parameter file with initial values of reservoirs extracted from a simulation file
modify_file.initial_pVs_Vo_Qc_from_simu('create_the_parameter_files/initial_pVs_Vo_Qc_from_simu.ini')

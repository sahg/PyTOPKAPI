import os

from ConfigParser import SafeConfigParser

import h5py
import numpy as np

import pytopkapi

old_settings = None # global variable for numpy error settings

def read_hdf5_array(fname, dset_string):
    h5file = h5py.File(fname)
    array = h5file[dset_string][...]
    h5file.close()

    return array

def compute_precip_volume(precip_fname, group_name, X):
    """Compute the volume of precipitation over the catchment.

    """
    dset_string = '/%s/rainfall' % group_name
    data = read_hdf5_array(precip_fname, dset_string)

    depth = data.sum()/1000.0 # convert mm --> m

    precip_vol = depth*(X**2)

    return precip_vol

def compute_evapot_volume(result_fname, X):
    dset_string = '/ET_out'
    data = read_hdf5_array(result_fname, dset_string)

    depth = data.sum()/1000.0 # convert mm --> m

    evapo_vol = depth*(X**2)

    return evapo_vol

def compute_evap_volume(result_fname, channel_indices):
    """Compute total channel evaporation.

    Compute the total open water evaporation from channel cells. This is
    reported as a volume in the TOPKAPI results file.

    """
    h5file = h5py.File(result_fname)
    data = h5file['/Channel/Ec_out'][:, channel_indices]
    h5file.close()

    return data.sum()

def compute_storage(result_fname):
    """Compute the initial and final storage.

    The total storage in the catchment is calculated at the start and end of
    the simulation period.

    """
    h5file = h5py.File(result_fname)

    # soil stores
    initial = h5file['/Soil/V_s'][0, :]
    final = h5file['/Soil/V_s'][-1, :]
    initial_storage = initial.sum()
    final_storage = final.sum()

    # overland stores
    initial = h5file['/Overland/V_o'][0, :]
    final = h5file['/Overland/V_o'][-1, :]
    initial_storage += initial.sum()
    final_storage += final.sum()

    # channel stores
    initial = h5file['/Channel/V_c'][0, :]
    final = h5file['/Channel/V_c'][-1, :]
    initial_storage += initial.sum()
    final_storage += final.sum()

    h5file.close()

    return initial_storage, final_storage

def compute_channel_runoff(result_fname, delta_t, cell_id):
    h5file = h5py.File(result_fname)

    dset = h5file['/Channel/Qc_out']
    flows = dset[1:, cell_id]
    runoff_vol = flows.sum() * delta_t

    h5file.close()

    return runoff_vol

def compute_overland_runoff(result_fname, delta_t, cell_id):
    h5file = h5py.File(result_fname)

    dset = h5file['/Overland/Qo_out']
    flows = dset[1:, cell_id]
    runoff_vol = flows.sum() * delta_t

    h5file.close()

    return runoff_vol

def compute_soil_drainage(result_fname, delta_t, cell_id):
    h5file = h5py.File(result_fname)

    flows = h5file['/Soil/Qs_out'][1:, cell_id]
    runoff_vol = flows.sum() * delta_t

    h5file.close()

    return runoff_vol

def compute_down_drainage(result_fname, delta_t, cell_id):
    h5file = h5py.File(result_fname)

    flows = h5file['/Q_down'][1:, cell_id]
    runoff_vol = flows.sum() * delta_t

    h5file.close()

    return runoff_vol

def continuity_error(ini_fname, delta_t, cell_id, X, channel_indices):
    # run the model using the supplied configuration
    pytopkapi.run(ini_fname)

    # parse the config file
    config = SafeConfigParser()
    config.read(ini_fname)

    precip_fname = config.get('input_files', 'file_rain')
    ET_fname = config.get('input_files', 'file_ET')
    group_name = config.get('groups', 'group_name')
    result_fname = config.get('output_files', 'file_out')

    # write model version
    print 'PyTOPKAPI version = %s' % pytopkapi.__version__

    # compute the terms in the continuity eqn.
    initial_storage, final_storage = compute_storage(result_fname)
    print 'Initial storage = ', initial_storage
    print 'Final storage = ', final_storage

    precip_vol = compute_precip_volume(precip_fname, group_name, X)
    print 'Precipitation = ', precip_vol

    evapo_vol = compute_evapot_volume(result_fname, X)
    print 'Evapotranspiration = ', evapo_vol

    open_water_evap_vol = compute_evap_volume(result_fname, channel_indices)
    print 'Channel evaporation = ', open_water_evap_vol

    channel_runoff_vol = compute_channel_runoff(result_fname, delta_t, cell_id)
    print 'Channel runoff (outlet) = ', channel_runoff_vol

    overland_runoff_vol = compute_overland_runoff(result_fname,
                                                  delta_t, cell_id)
    print 'Overland runoff (outlet) = ', overland_runoff_vol

    soil_drainage_vol = compute_soil_drainage(result_fname, delta_t, cell_id)
    print 'Soil drainage (outlet) = ', soil_drainage_vol

    down_drainage_vol = compute_down_drainage(result_fname, delta_t, cell_id)
    print 'Non-channel drainage (outlet) = ', down_drainage_vol

    input = precip_vol
    output = evapo_vol \
             + open_water_evap_vol \
             + channel_runoff_vol \
             + down_drainage_vol

    delta_storage = final_storage - initial_storage
    error = delta_storage - (input - output)

    if precip_vol > 0:
        precip_error = abs((error/precip_vol)*100.0)
    else:
        precip_error = None
    stor_error = abs((error/initial_storage)*100.0)

    print 'Continuity error = ', error
    print 'Error as % precip. = ', precip_error
    print 'Error as % initial storage = ', stor_error

    os.remove(result_fname)

    return error, precip_error, stor_error

# Continuity tests and test environment
def setup():
    "set up test fixtures"
    os.chdir(os.path.join(os.getcwd(), 'pytopkapi/tests/test_continuity'))

    old_settings = np.seterr(all='ignore')

def teardown():
    "tear down test fixtures"
    os.chdir('../../..')
    np.seterr(all=old_settings)

def test_4cell_continuity():
    """Test continuity on the 4 cell catchment (Parak, 2006).

    """
    ini_fname = '4cells.ini'
    delta_t = 3600.0
    cell_id = 3
    X = 1000.0
    channel_indices = [1, 2, 3]

    error, precip_error, stor_error = continuity_error(ini_fname,
                                                       delta_t,
                                                       cell_id, X,
                                                       channel_indices)

    assert precip_error < 2.8e-05
    assert stor_error < 2.1e-03

def test_d8_continuity():
    """Test continuity on a generic catchment.

    The cell connectivity is D8 and evaporative forcing is applied in addition
    to rainfall.

    """
    ini_fname = 'D8.ini'
    delta_t = 21600.0
    cell_id = 0
    X = 1000.0
    channel_indices = [0, 1, 4]

    error, precip_error, stor_error = continuity_error(ini_fname,
                                                       delta_t,
                                                       cell_id, X,
                                                       channel_indices)

    assert precip_error < 3.6e-04
    assert stor_error < 1.3e-03

def test_lieb_continuity():
    """Test continuity on a sub-catchment of Liebenbergsvlei.

    The cell connectivity is D8, rainfall and evaporative forcing are
    both zero.

    """
    ini_fname = 'lieb.ini'
    delta_t = 21600.0
    cell_id = 0
    X = 1000.0
    channel_indices = [0]

    error, precip_error, stor_error = continuity_error(ini_fname,
                                                       delta_t,
                                                       cell_id, X,
                                                       channel_indices)
    assert precip_error == None
    assert stor_error < 1.5e-05

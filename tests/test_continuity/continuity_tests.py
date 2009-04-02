import os

from ConfigParser import SafeConfigParser

import tables as h5
import numpy as np

import TOPKAPI.model as topkapi

def compute_precip_volume(precip_fname, group_name, X):
    """Compute the volume of precipitation over the catchment.

    """
    h5file = h5.openFile(precip_fname)
    node = h5file.getNode('/'+group_name, 'rainfall')
    data = node.read()
    h5file.close()

    depth = data.sum()/1000.0 # convert mm --> m

    precip_vol = depth*(X**2)

    return precip_vol

def compute_evapot_volume(result_fname, X):
    h5file = h5.openFile(result_fname)
    node = h5file.getNode('/', 'ET_out')
    data = node.read()
    h5file.close()

    depth = data.sum()/1000.0 # convert mm --> m

    evapo_vol = depth*(X**2)

    return evapo_vol

def compute_evap_volume(result_fname, channel_indices):
    """Compute total channel evaporation.

    Compute the total open water evaporation from channel cells. This is
    reported as a volume in the TOPKAPI results file.

    """
    h5file = h5.openFile(result_fname)
    node = h5file.getNode('/'+'Channel', 'Ec_out')
    data = node.read()[:, channel_indices]
#    data = node.read()
    h5file.close()

    return data.sum()

def compute_storage(result_fname):
    """Compute the initial and final storage.

    The total storage in the catchment is calculated at the start and end of
    the simulation period.

    """
    h5file = h5.openFile(result_fname)

    # soil stores
    node = h5file.getNode('/'+'Soil', 'V_s')
    initial = node.read()[0, :]
    final = node.read()[-1, :]
    initial_storage = initial.sum()
    final_storage = final.sum()

    # overland stores
    node = h5file.getNode('/'+'Overland', 'V_o')
    initial = node.read()[0, :]
    final = node.read()[-1, :]
    initial_storage += initial.sum()
    final_storage += final.sum()

    # channel stores
    node = h5file.getNode('/'+'Channel', 'V_c')
    initial = node.read()[0, :]
    final = node.read()[-1, :]
    initial_storage += initial.sum()
    final_storage += final.sum()

    h5file.close()

    return initial_storage, final_storage

def compute_channel_runoff(result_fname, delta_t, cell_id):
    h5file = h5.openFile(result_fname)

    node = h5file.getNode('/'+'Channel', 'Qc_out')
    flows = node.read()[1:, cell_id]
    runoff_vol = flows.sum() * delta_t

    h5file.close()

    return runoff_vol

def compute_overland_runoff(result_fname, delta_t, cell_id):
    h5file = h5.openFile(result_fname)

    node = h5file.getNode('/'+'Overland', 'Qo_out')
    flows = node.read()[1:, cell_id]
    runoff_vol = flows.sum() * delta_t

    h5file.close()

    return runoff_vol

def compute_soil_drainage(result_fname, delta_t, cell_id):
    h5file = h5.openFile(result_fname)

    node = h5file.getNode('/'+'Soil', 'Qs_out')
    flows = node.read()[1:, cell_id]
    runoff_vol = flows.sum() * delta_t

    h5file.close()

    return runoff_vol

def compute_down_drainage(result_fname, delta_t, cell_id):
    h5file = h5.openFile(result_fname)

    node = h5file.getNode('/Q_down')
    flows = node.read()[1:, cell_id]
    runoff_vol = flows.sum() * delta_t

    h5file.close()

    return runoff_vol

def continuity_error(ini_fname, delta_t, cell_id, X, channel_indices):
    # run the model using the supplied configuration
    topkapi.run(ini_fname)

    # parse the config file
    config = SafeConfigParser()
    config.read(ini_fname)

    precip_fname = config.get('input_files', 'file_rain')
    ET_fname = config.get('input_files', 'file_ET')
    group_name = config.get('groups', 'group_name')
    result_fname = config.get('output_files', 'file_out')

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

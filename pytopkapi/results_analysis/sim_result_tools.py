__all__ = ['extract_ssi', 'extract_ssi_to_file',
           'extract_Q_channel', 'extract_Q_down']

from datetime import timedelta
from ConfigParser import SafeConfigParser

import h5py
import numpy as np
import numpy.ma as ma

# gzip compression flag
comp = 6

def extract_Q_down(control_fname):
    """Extract combined soil and overland out flow rates.

    Read a PyTOPKAPI simulation file and return the combined overland
    andsoil store outflows in a Numpy array.

    Parameters
    ----------
    control_fname : string
        The file name of a PyTOPKAPI simulation control file. The name
        should contain the full path relative to the current
        directory.

    Returns
    -------
    Qdown : Numpy array
        A Numpy array containing the simulated outflow flow rates from
        the overland and soil store of each cell.

    """
    config = SafeConfigParser()
    config.read(control_fname)

    sim_fname = config.get('output_files', 'file_out')

    tkpi_file = h5py.File(sim_fname)
    Qdown = tkpi_file['/Q_down'][...]
    tkpi_file.close()

    return Qdown

def extract_Q_channel(control_fname):
    """Extract channel flow rates from a PyTOPKAPI simulation file.

    Read a PyTOPKAPI simulation file and return the simulated channel
    flows in a Numpy masked array.

    Parameters
    ----------
    control_fname : string
        The file name of a PyTOPKAPI simulation control file. The name
        should contain the full path relative to the current
        directory.

    Returns
    -------
    Qc : Numpy masked array
        A Numpy masked array containing the simulated flow rates for
        channel cells.

    """
    config = SafeConfigParser()
    config.read(control_fname)

    param_fname = config.get('input_files', 'file_cell_param')
    sim_fname = config.get('output_files', 'file_out')

    params = np.loadtxt(param_fname)

    tkpi_file = h5py.File(sim_fname)
    Qc = tkpi_file['/Channel/Qc_out'][...]
    tkpi_file.close()

    channel_mask = params[:, 3]
    cond = params[:, 3]*np.ones(Qc.shape, dtype=np.int) != 1

    Qc = np.ma.masked_where(cond, Qc)

    return Qc

def extract_ssi(control_fname):
    """Extract SSI from a PyTOPKAPI simulation file.

    Read a PyTOPKAPI simulation file and it's associated parameter
    file and compute the Soil Saturation Index (SSI) for each model
    cell and timestep. The results are returned as a Numpy array.

    Parameters
    ----------
    control_fname : string
        The file name of a PyTOPKAPI simulation control file. The name
        should contain the full path relative to the current
        directory.

    Returns
    -------
    ssi : Numpy ndarray
        A Numpy array containing the calculated SSI values.

    """
    config = SafeConfigParser()
    config.read(control_fname)

    global_param_fname = config.get('input_files', 'file_global_param')
    param_fname = config.get('input_files', 'file_cell_param')
    sim_fname = config.get('output_files', 'file_out')
    fac_L = config.getfloat('calib_params', 'fac_L')

    params = np.loadtxt(param_fname)
    glob_params = np.genfromtxt(global_param_fname, names=True)

    soil_depth = fac_L*params[:, 8]
    factor = params[:, 11] - params[:, 10]
    cell_area = glob_params['X']**2 # m^2

    soil_depth = ma.masked_values(soil_depth, 0.0)
    factor = ma.array(factor, mask=soil_depth.mask)
    div = factor*soil_depth*cell_area

    tkpi_file = h5py.File(sim_fname)
    soil_vol = tkpi_file['/Soil/V_s'][...]
    tkpi_file.close()

    # ssi = (Vs/cell_vol)*100
    # cell_vol = (theta_s - theta_r)*soil_depth*cell_area
    sv = ma.array(soil_vol, mask=soil_depth.mask)
    ssi = (sv/(div))*100.0

    return ssi

def extract_ssi_to_file(sim_fname, param_fname, result_fname, start_dt):
    """
    Read a TOPKAPI simulation file and it's associated parameter file
    and compute the SSI for each timestep. Store the results in a new
    HDF5 file, grouped by date and containing datasets of latitude,
    longitude and SSI value.

    """
    params = np.loadtxt(param_fname)
    lon = params[:, 1]
    lat = params[:, 2]
    
    soil_depth = params[:, 8]
    factor = params[:, 11] - params[:, 10]
    cell_area = 1000.0**2 # m^2

    soil_depth = ma.masked_values(soil_depth, 0.0)
    factor = ma.array(factor, mask=soil_depth.mask)
    lon = ma.array(lon, mask=soil_depth.mask).compressed()
    lat = ma.array(lat, mask=soil_depth.mask).compressed()

    div = factor*soil_depth*cell_area

    tkpi_file = h5py.File(sim_fname)
    result_file = h5py.File(result_fname, 'w')

    soil_vol = tkpi_file['/Soil/V_s'][...]
    tkpi_file.close()
    rows, cols = soil_vol.shape

    # lat
    dset = result_file.require_dataset('lat', shape=lat.shape,
                                       dtype=np.float32, compression=comp)
    dset[...] = lat

    dset.attrs['name'] = 'latitude'
    dset.attrs['units'] = 'Decimal degrees'

    # lon
    dset = result_file.require_dataset('lon', shape=lon.shape,
                                       dtype=np.float32, compression=comp)
    dset[...] = lon

    dset.attrs['name'] = 'longitude'
    dset.attrs['units'] = 'Decimal degrees'

    curr_dt = start_dt
    for k in range(rows):
        print curr_dt
        # ssi = (Vs/cell_vol)*100
        # cell_vol = (theta_s - theta_r)*soil_depth*cell_area
        sv = ma.array(soil_vol[k], mask=soil_depth.mask)
        ssi = (sv/(div))*100.0

        ssi = ssi.compressed()

        # ssi
        dset = result_file.require_dataset(curr_dt.strftime('%Y%m%d%H00'),
                                           shape=ssi.shape,
                                           dtype=np.float32, compression=comp)
        dset[...] = ssi

        dset.attrs['name'] = 'TOPKAPI soil saturation index'
        dset.attrs['units'] = '% saturation'

        curr_dt += timedelta(hours=3)
    
    result_file.close()

__all__ = ['extract_ssi', 'extract_ssi_to_file',
           'extract_eta', 'extract_eta_to_file',
           'extract_Q_channel', 'extract_Q_down',
           'extract_overland_volume', 'extract_overland_volume_to_file']

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

def extract_overland_volume(control_fname):
    """Extract the volumes in the overland stores.

    Read a PyTOPKAPI simulation file and return the combined overland
    and store volumes in a Numpy array.

    Parameters
    ----------
    control_fname : string
        The file name of a PyTOPKAPI simulation control file. The name
        should contain the full path relative to the current directory
        (or the root of the file system).

    Returns
    -------
    Vo : Numpy array
        A Numpy array containing the simulated storage volume in the
        overland store of each cell.

    """
    config = SafeConfigParser()
    config.read(control_fname)

    sim_fname = config.get('output_files', 'file_out')

    tkpi_file = h5py.File(sim_fname)
    Vo = tkpi_file['/Overland/V_o'][...]
    tkpi_file.close()

    return Vo

def extract_overland_volume_to_file(sim_fname, param_fname,
                                    result_fname, start_dt, timestep):
    """Extract the volumes in the overland stores to a file.

    Read a TOPKAPI simulation file and it's associated parameter file
    and extract the overland store volumes for each timestep. Store
    the results in a new HDF5 file, grouped by date and containing
    datasets of latitude, longitude and storage volume.

    Parameters
    ----------
    sim_fname : string
        The name of a PyTOPKAPI simulation file. This should include
        the full or relative path.
    param_fname : string
        The name of a parameter file describing the catchment. This
        should include the full or relative path.
    result_fname : string
        The name of an HDF5 file to store the output. This should
        include the full or relative path.
    start_dt : datetime.datetime
        The starting date and time of the simulated results in
        `sim_fname`.
    timestep : int
        The length of each model time-step in seconds.

    Returns
    -------
    Nothing

    """
    params = np.loadtxt(param_fname)
    x = params[:, 1]
    y = params[:, 2]
    soil_depth = params[:, 8]

    soil_depth = ma.masked_values(soil_depth, 0.0)
    x = ma.array(x, mask=soil_depth.mask).compressed()
    y = ma.array(y, mask=soil_depth.mask).compressed()

    tkpi_file = h5py.File(sim_fname)
    result_file = h5py.File(result_fname, 'w')

    overland_vol = tkpi_file['/Overland/V_o'][...]
    tkpi_file.close()
    rows, cols = overland_vol.shape

    # y
    dset = result_file.require_dataset('y', shape=y.shape,
                                       dtype=np.float32, compression=comp)
    dset[...] = y

    dset.attrs['name'] = 'y coordinate'
    dset.attrs['units'] = 'Projection dependent (Metres or Decimal degrees)'

    # x
    dset = result_file.require_dataset('x', shape=x.shape,
                                       dtype=np.float32, compression=comp)
    dset[...] = x

    dset.attrs['name'] = 'x coordinate'
    dset.attrs['units'] = 'Projection dependent (Metres or Decimal degrees)'

    curr_dt = start_dt
    for k in range(rows):
        print curr_dt

        ov = ma.array(overland_vol[k], mask=soil_depth.mask).compressed()

        dset = result_file.require_dataset(curr_dt.strftime('%Y%m%d%H00'),
                                           shape=ov.shape,
                                           dtype=np.float32, compression=comp)
        dset[...] = ov

        dset.attrs['name'] = 'TOPKAPI overland store volume'
        dset.attrs['units'] = 'm^3'

        curr_dt += timedelta(seconds=timestep)

    tkpi_file.close()
    result_file.close()

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

def extract_ssi_to_file(sim_fname, param_fname,
                        result_fname, start_dt, timestep):
    """Extract percentage saturation to a file

    Read a TOPKAPI simulation file and it's associated parameter file
    and compute the SSI for each timestep. Store the results in a new
    HDF5 file, grouped by date and containing datasets of latitude,
    longitude and SSI value.

    Parameters
    ----------
    sim_fname : string
        The name of a PyTOPKAPI simulation file. This should include
        the full or relative path.
    param_fname : string
        The name of a parameter file describing the catchment. This
        should include the full or relative path.
    result_fname : string
        The name of an HDF5 file to store the output. This should
        include the full or relative path.
    start_dt : datetime.datetime
        The starting date and time of the simulated results in
        `sim_fname`.
    timestep : int
        The length of each model time-step in seconds.

    Returns
    -------
    Nothing

    """
    params = np.loadtxt(param_fname)
    x = params[:, 1]
    y = params[:, 2]

    soil_depth = params[:, 8]
    factor = params[:, 11] - params[:, 10]
    cell_area = 1000.0**2 # m^2

    soil_depth = ma.masked_values(soil_depth, 0.0)
    factor = ma.array(factor, mask=soil_depth.mask)
    x = ma.array(x, mask=soil_depth.mask).compressed()
    y = ma.array(y, mask=soil_depth.mask).compressed()

    div = factor*soil_depth*cell_area

    tkpi_file = h5py.File(sim_fname)
    result_file = h5py.File(result_fname, 'w')

    soil_vol = tkpi_file['/Soil/V_s'][...]
    tkpi_file.close()
    rows, cols = soil_vol.shape

    # y
    dset = result_file.require_dataset('y', shape=y.shape,
                                       dtype=np.float32, compression=comp)
    dset[...] = y

    dset.attrs['name'] = 'y coordinate'
    dset.attrs['units'] = 'Projection dependent (Metres or Decimal degrees)'

    # x
    dset = result_file.require_dataset('x', shape=x.shape,
                                       dtype=np.float32, compression=comp)
    dset[...] = x

    dset.attrs['name'] = 'x coordinate'
    dset.attrs['units'] = 'Projection dependent (Metres or Decimal degrees)'

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

        curr_dt += timedelta(seconds=timestep)

    tkpi_file.close()
    result_file.close()

def extract_eta(control_fname):
    """Extract ETa from a PyTOPKAPI simulation file.

    Read a PyTOPKAPI simulation file and it's associated parameter file
    and extract the actual evapotranspiration for each model cell and
    timestep. The results are returned as a Numpy array.

    Parameters
    ----------
    control_fname : string
        The file name of a PyTOPKAPI simulation control file. The name
        should contain the full path relative to the current
        directory.

    Returns
    -------
    eta : Numpy ndarray
        A Numpy array containing the calculated ETa values.

    """
    config = SafeConfigParser()
    config.read(control_fname)

    param_fname = config.get('input_files', 'file_cell_param')
    sim_fname = config.get('output_files', 'file_out')

    params = np.loadtxt(param_fname)

    soil_depth = params[:, 8]
    soil_depth = ma.masked_values(soil_depth, 0.0)

    tkpi_file = h5py.File(sim_fname)
    eta = tkpi_file['/ET_out'][...]
    tkpi_file.close()

    eta = ma.array(eta, mask=soil_depth.mask)

    return eta

def extract_eta_to_file(sim_fname, param_fname,
                        result_fname, start_dt, timestep):
    """Extract actual evapotranspiration to a file

    Read a PyTOPKAPI simulation file and it's associated parameter
    file and extract the actual evapotranspiration for each
    timestep. Store the results in a new HDF5 file, grouped by date
    and containing datasets of latitude, longitude and ETa value.

    Parameters
    ----------
    sim_fname : string
        The name of a PyTOPKAPI simulation file. This should include
        the full or relative path.
    param_fname : string
        The name of a parameter file describing the catchment. This
        should include the full or relative path.
    result_fname : string
        The name of an HDF5 file to store the output. This should
        include the full or relative path.
    start_dt : datetime.datetime
        The starting date and time of the simulated results in
        `sim_fname`.
    timestep : int
        The length of each model time-step in seconds.

    Returns
    -------
    Nothing

    """
    params = np.loadtxt(param_fname)
    x = params[:, 1]
    y = params[:, 2]

    soil_depth = params[:, 8]
    soil_depth = ma.masked_values(soil_depth, 0.0)

    x = ma.array(x, mask=soil_depth.mask).compressed()
    y = ma.array(y, mask=soil_depth.mask).compressed()

    tkpi_file = h5py.File(sim_fname)
    result_file = h5py.File(result_fname, 'w')

    eta = tkpi_file['/ET_out'][...]
    tkpi_file.close()

    rows, cols = eta.shape

    # y
    dset = result_file.require_dataset('y', shape=y.shape,
                                       dtype=np.float32, compression=comp)
    dset[...] = y

    dset.attrs['name'] = 'y coordinate'
    dset.attrs['units'] = 'Projection dependent (Metres or Decimal degrees)'

    # x
    dset = result_file.require_dataset('x', shape=x.shape,
                                       dtype=np.float32, compression=comp)
    dset[...] = x

    dset.attrs['name'] = 'x coordinate'
    dset.attrs['units'] = 'Projection dependent (Metres or Decimal degrees)'

    curr_dt = start_dt
    for k in range(rows):
        print curr_dt

        et = ma.array(eta[k], mask=soil_depth.mask)
        et = et.compressed()

        # ETa
        dset = result_file.require_dataset(curr_dt.strftime('%Y%m%d%H00'),
                                           shape=et.shape,
                                           dtype=np.float32, compression=comp)
        dset[...] = et

        dset.attrs['name'] = 'PyTOPKAPI actual ET'
        dset.attrs['units'] = 'mm'

        curr_dt += timedelta(seconds=timestep)

    result_file.close()

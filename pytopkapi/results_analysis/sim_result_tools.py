__all__ = ['extract_ssi']

from datetime import datetime, timedelta

import h5py
import numpy as np
import numpy.ma as ma

# gzip compression flag
comp = 6

def extract_ssi(sim_fname, param_fname, result_fname, start_dt):
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

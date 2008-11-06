import numpy as np
import rpy
from rpy import r

from TOPKAPI import arcfltgrid

def krige_to_grid(grid_fname, obs_x, obs_y, obs_data, vgm_par):
    """Interpolate point data onto a grid using Kriging.
    
    Interpolate point data onto a regular rectangular grid of square cells using
    Kriging with a predefined semi-variogram.  The observed data locations must
    be specified in the same projection and coordinate system as the grid, which
    is defined in an ArcGIS raster file.
    
    Parameters
    ----------
    grid_fname : string
        Filename of an ArcGIS float grid raster defining the required grid to
        Krige onto.  All cells are included regardless of their value.
    obs_x : array_like
        The x coordinates of the observation locations.
    obs_y : array_like
        The y coordinates of the observation locations.
    obs_data : array_like
        The data values at the observation locations.
    vgm : dict
        A dictionary describing the semi-variogram model.  Required keys are:
        'model' can be one of {'Lin', 'Exp', 'Sph', 'Gau'}
        'nugget' must be a scalar
        'range' must be a scalar
        'sill' must be a scalar
    
    Returns
    -------
    kriged_est : 2darray
        A 2D array containing the Kriged estimates at each point on the 
        specified rectangular grid.
    
    Notes
    -----
    This function requires that R, RPy and the R gstat library are correctly
    installed.
    
    """
    grid, headers = arcfltgrid.read(grid_fname)
    cols = headers[0]
    rows = headers[1]
    x0 = headers[2]
    y0 = headers[3]
    cell_size = headers[4]
    
    # define the grid (pixel centre's)
    xt, yt = np.meshgrid(np.linspace(x0, x0 + (cols-1)*cell_size, num=cols),
                         np.linspace(y0, y0 + (rows-1)*cell_size, num=rows))
    
    xt = xt.flatten()
    yt = yt.flatten()
    
    # Krige using gstat via RPy
    r.library('gstat')
    rpy.set_default_mode(rpy.NO_CONVERSION)
    
    obs_frame = r.data_frame(x=obs_x, y=obs_y, data=obs_data)
    target_grid = r.data_frame(x=xt, y=yt)
    
    v = r.vgm(vgm_par['sill'], vgm_par['model'], 
              vgm_par['range'], vgm_par['nugget'])
    
    result = r.krige(r('data ~ 1'), r('~ x + y'), 
                     obs_frame, target_grid, model=v)
    
    rpy.set_default_mode(rpy.BASIC_CONVERSION)
    
    result = result.as_py()
    
    kriged_est = np.array(result['var1.pred'])
    kriged_est = kriged_est.reshape(rows, cols)
    kriged_est = np.flipud(kriged_est)
    
    return kriged_est

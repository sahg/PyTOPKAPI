"""
*** OBJECTIVE

1. creating the TOPKAPI parameter file from GIS binary files (binary
grid format).  -main routine 'creat_param_file'

*** COMMENT

In the present program, to write the parameter file,
generally the grid are tranformed in an array by rearranging the cell
order from the West to East and North to South, as in the following
example:

The grid format is like:

GIS bingrid file
-9999  -9999  -9999  -9999  -9999  -9999  -9999
-9999      0      1  -9999  -9999  -9999  -9999
-9999      2      3     4   -9999  -9999  -9999
    5      6      7     8       9  -9999  -9999
-9999     10     11    12      13     14  -9999
-9999  -9999     15    16      17  -9999  -9999
-9999  -9999  -9999  -9999     18  -9999  -9999
-9999  -9999  -9999 -9999   -9999  -9999  -9999

The corresponding array extracted and ordered from West to East, North
to South is: 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18

"""
import sys
from warnings import warn
from ConfigParser import SafeConfigParser

#Python modules
import numpy as np
import pylab as pl
from numpy import ma
from osgeo import gdal

#External modules from PyTOPKAPI
#Utilities
from pytopkapi import utils as ut
#pretreatment: used for subroutines to read the column type files.
from pytopkapi import pretreatment as pm

def compute_cell_coordinates(mask_fname):
    dset = gdal.Open(mask_fname)
    mask = dset.ReadAsArray()

    x0, dx, fy, y0, fx, dy = dset.GetGeoTransform()

    # Adjust x0 and y0 to give pixel centres. GDAL specifies (x0, y0)
    # as the top left corner of the top left pixel. PyTOPKAPI expects
    # pixel centres. At the centre of the first pixel the following
    # holds: (Xpixel, Yline) == (0.5, 0.5)
    x0 = x0 + dx/2.0 + fy/2.0
    y0 = y0 + fx/2.0 + dy/2.0

    Yline, Xpixel = np.nonzero(mask == 1)

    Xgeo = x0 + Xpixel*dx + Yline*fy
    Ygeo = y0 + Xpixel*fx + Yline*dy

    return Xgeo, Ygeo

def read_raster(rast_fname, file_format='GTiff'):
    """Read the data in a raster file

    Parameters
    ----------
    rast_fname : string
        The path to the raster file.
    file_format : string
        The file format of the raster file, currently this can only be
        GeoTIFF ('GTiff' - default).

    Returns
    -------
    data : Numpy ndarray
        The 2D raster data in a Numpy array. The dtype of the returned
        array is the same as the data type stored in the raster file.

    """
    if file_format != 'GTiff':
        err_str = 'Reading %s files not implemented.' % file_format
        raise NotImplementedError(err_str)
    else:
        dset = gdal.Open(rast_fname)
        data = dset.ReadAsArray()

    return data

def generate_param_file(ini_fname, isolated_cells=False):
    """Create a PyTOPKAPI parameter file

    Generates a parameter file from the catchment data provided in a
    set of georeferenced raster files. The input files may be GeoTIFF
    or 32-bit raster files with ArcGIS style headers.

    Parameters
    ----------
    ini_fname : string
        The path to an ini-style config file specifying the input data
        file locations and format.

    Returns
    -------
    Nothing

    Notes
    -----

    * Input:
      +Binary grid files (from GIS)

      - file_bin_streamnet: grid of the streamnet (code: 1 for channel
        cell, 0 for hillslope cell)

      - file_bin_beta: grid of cell slopes in degress.

      - file_bin_flowdir: grid of flow directions

      - file_bin_GLCC: grid file containing the GLCC land use/cover
        codes

      - file_bin_SIRI: grid of soil properties (code given by SIRI)

      - file_bin_WRC90: binary grid file containing the WRC90 soil
        property codes (Here only three codes are considered
        (cf. comments): 3 for Loamy Sand, 2 for Sandy Loam, 1 for
        Clay)

      - file_bin_strahler: grid of strahler order of the channel
        cells.

      +Tables from litterature for transformaing grids into TOPKAPI
      physical parameters

      - file_table_GLCC_manning: an ASCII file containing a table of
                                correspondance between the GLCC land
                                use/cover codes and the values of
                                manning's coef proposed in different
                                references (Chow et al., 1998;
                                Maidment,1993 - give a range of values
                                and the MUSIC report)

      - file_table_SIRI_soil: an ASCII file containing a table of
                              correspondance between the SIRI codes
                              and the values of a selection of the
                              soil parameters proposed by SIRI
                              (Land-type, depth A, depth B, WP A, WP
                              B, FC A, FC B, Por A, Por B)

      - file_table_WRC90_soil: an ASCII file containing a table of
                               correspondance between the WRC90 codes
                               and the values of Ks (permeability) and
                               theta_r (residual soil moisture) from
                               Maidment(1993)

      - file_table_strahler_manning: an ASCII file containing a table
                                     of correspondance between the
                                     strahler order and the values of
                                     manning strickler proposed by Liu
                                     and Todini (2002)

      +Constant parameters

      - Vs_t0: Constant value for the initial soil saturation of each
        cell.

      - Vc_t0: Constant value for the initial channel water

      - kc: Crop coefficient

    * Output:

      - file_out: parameter file (ASCII column format) containing:
       label X Y lambda Xc dam tan_beta L Ks Theta_r Theta_s n_o n_c
       cell_down pVs_t0 Vo_t0 Qc_t0 kc

    * Comment:

     1. !!!!VERY IMPORTANT!!!! The routine refers to several
     subroutines listed below that must be carefully read before
     running the programm. This programm is helpfull for creating the
     parameter file but it is not automated.  Some valuable
     informations are required inside these subroutines especially for
     data that are assigned using Tables.  The Tables are indeed
     simple ASCII files that were created for the special case of the
     Liebenbergsvlei catchment.  The application to another catchment
     might require the modification of the Tables and thus the
     subroutines.  One must refer to the headers of each subroutine
     for the detailed information.

     2. Note that here, the initial soil moisture value, as well as
     the initial channel saturation are constant.  Assigning a
     spatially variable initial soil moisture can be done through the
     routine 'from_param_to_new_param_catchVsi' in this file.

    """
    config = SafeConfigParser()
    config.read(ini_fname)

    mask_fname = config.get('raster_files', 'mask_fname')
    soil_depth_fname = config.get('raster_files', 'soil_depth_fname')
    conductivity_fname = config.get('raster_files', 'conductivity_fname')
    hillslope_fname = config.get('raster_files', 'hillslope_fname')
    theta_sat_fname = config.get('raster_files', 'sat_moisture_content_fname')
    theta_r_fname = config.get('raster_files', 'resid_moisture_content_fname')
    psi_b_fname = config.get('raster_files', 'bubbling_pressure_fname')
    lamda_fname = config.get('raster_files', 'pore_size_dist_fname')
    n_o_fname = config.get('raster_files', 'overland_manning_fname')
    network_fname = config.get('raster_files', 'channel_network_fname')
    flowdir_fname = config.get('raster_files', 'flowdir_fname')
    fdir_source = config.get('raster_files', 'flowdir_source')

    param_fname = config.get('output', 'param_fname')

    # Read the input rasters
    mask = read_raster(mask_fname)
    hillslope = read_raster(hillslope_fname)
    depth = read_raster(soil_depth_fname)
    theta_sat = read_raster(theta_sat_fname)
    theta_r = read_raster(theta_r_fname)
    conductivity = read_raster(conductivity_fname)
    psi_b = read_raster(psi_b_fname)
    lamda = read_raster(lamda_fname)
    n_o = read_raster(n_o_fname)
    channel_network = read_raster(network_fname)
    flowdir = read_raster(flowdir_fname)

    # Calculate parameters
    ncells = mask[mask == 1].size
    nparams = 21
    cell_labels = np.arange(ncells)
    tan_beta = np.tan((np.pi/180.0)*hillslope)
    X, Y = compute_cell_coordinates(mask_fname)

    channel_network[channel_network < 255] = 1
    channel_network[channel_network == 255] = 0

    if isolated_cells == True:
        cell_down = -999
        channel_length = 0
    else:
        # Calculate the network connections and channel lengths.
        cell_down = cell_connectivity(flowdir, mask, fdir_source)
        channel_length = compute_Xchannel(cell_labels,
                                          channel_network[mask == 1],
                                          X, Y, cell_down)

    # Write parameter file
    param_table = np.zeros((ncells, nparams))
    param_table[:,0] = cell_labels
    param_table[:,1] = X
    param_table[:,2] = Y
    param_table[:,3] = channel_network[mask == 1]
    param_table[:,4] = channel_length
    ## param_table[:,5] = dam_locations
    param_table[:,6] = tan_beta[mask == 1]
    ## param_table[:,7] = tan_beta_channel
    param_table[:,8] = depth[mask == 1]
    param_table[:,9] = conductivity[mask == 1]
    param_table[:,10] = theta_r[mask == 1]
    param_table[:,11] = theta_sat[mask == 1]
    param_table[:,12] = n_o[mask == 1]
    ## param_table[:,13] = n_c[mask == 1]
    param_table[:,14] = cell_down
    ## param_table[:,15] = pVs_t0
    ## param_table[:,16] = Vo_t0
    ## param_table[:,17] = Qc_t0
    param_table[:,18] = 1 # Kc
    param_table[:,19] = psi_b[mask == 1]
    param_table[:,20] = lamda[mask == 1]

    format = '%d %f %f %d %f %d %f %f %f %f %f %f %f %f %d %f %f %f %f %f %f'
    np.savetxt(param_fname, param_table, fmt=format)

##############################################
###  SUBROUTINE USED IN "creat_param_file" ###
##############################################

####### READING THE BINARY FILES ###########

def read_bin_data(fname):
    """Read data from a float32 binary file into an array.

    """

    f = open(fname, "rb")
    raw = f.read()
    f.close()
    data = np.fromstring(raw, 'f')
    if sys.byteorder == 'big':
        data = data.byteswap()

    return data

def read_arc_bin(bingrid_name):
    """Read data from an ArcGIS float32 binary file into an array.

    """

    li_headers=read_headers_arc_bin(bingrid_name)

    rows = li_headers[1] # fixed for UM grid
    cols = li_headers[0] # fixed for UM grid

    bin_name = bingrid_name + '.flt'

    a = read_bin_data(bin_name)

    a = a.reshape(rows, cols)

    return a

def read_headers_arc_bin(bingrid_name):
    """Read the ascii headers of the binary grid file:
    ncols         62
    nrows         121
    xllcorner     -288595.47161281
    yllcorner     -3158065.5722693
    cellsize      1000
    NODATA_value  -9999
    byteorder     LSBFIRST
    """

    hdr_name = bingrid_name + '.hdr'
    f=open(hdr_name,'r')
    tab_read=f.readlines()
    f.close()

    li_headers=[]
    i=-1
    for line in tab_read:
        i=i+1
        donnees=line.split()
        if i<6:
            li_headers.append(float(donnees[1]))
        else:
            li_headers.append(donnees[1])

    return li_headers

def arc_bin_plot(bin_name, fig_name, title='Plot'):
    """Create a plot of the data in a binary file."""

    a = read_arc_bin(bin_name)

    a_mask = ma.masked_where(a < 0, a)
    pl.imshow(a_mask, interpolation='nearest')
    pl.colorbar()
    pl.title(title)
    pl.savefig(fig_name)
    pl.close()

############# OTHER SUBROUTINES ###################
def from_grid_to_param(file_bin_grid):
    """
    * Objective

    Directly extract the values at the cell locations from the grid

    * Input:

    GIS binary grid

    * Output

    A 1D array with nb_cell values. Cells are ordered from West to
    East, North to South.

    """
    tab=read_arc_bin(file_bin_grid)
    nrows=np.shape(tab)[0]
    ncols=np.shape(tab)[1]
    tab=np.reshape(tab,ncols*nrows)
    ind=np.where(tab>-99)
    ar_param=tab[ind]

    return ar_param

def from_GLCC_to_manning(file_bin_GLCC,file_table_GLCC_manning):
    """
    * Objective:

      Estimation of the overland manning's coefficient for each
      catchment cell from the GLCC land cover/use map

    * Input

      - file_bin_GLCC is the binary grid file containing the GLCC land
        use/cover codes

      - file_table_GLCC_manning is an ASCII file containing a table of
      correspondance between the GLCC land use/cover codes and the
      values of manning's coef proposed in different references (Chow
      et al., 1998; Maidment,1993 - give a range of values and the
      MUSIC report)

    * Ouput

      This routine returns a 1D array (ar_n_o) containing the values
      of the Manning's coefficient for each cell. Cells are ordered
      from West to East, North to South.

    * Comment:

      This routine has to be used carefully since the code is
      dependant to:

       1. the format of the Table (file_table_GLCC_manning)

       2. the choice of the user to use a given reference

    """
    #Read the binary grid file of GLCC land use type
    tab=read_arc_bin(file_bin_GLCC)
    nrows=np.shape(tab)[0]
    ncols=np.shape(tab)[1]
    tab=np.reshape(tab,ncols*nrows)
    ind=np.where(tab>-99)
    ar_GLCC=tab[ind]

    #Read the Table file within a header line
    tab=pm.read_column_input(file_table_GLCC_manning,5)
    ar_code=tab[:,0]
    ar_n_chow=tab[:,1]
    ar_n_maidment_inf=tab[:,2]
    ar_n_maidment_sup=tab[:,3]
    ar_n_music=tab[:,4]

    ar_n_o=np.array(ar_GLCC)
    for i in ar_code:
        ind=np.where(ar_GLCC==i)
        #!!!! TO BE CHANGED ACCORDING TO USER CHOICE !!!#
        if i>1:
            ar_n_o[ind]=ar_n_chow[np.where(ar_code==i)][0]
        if i==1:
            ar_n_o[ind]=ar_n_music[np.where(ar_code==i)][0]

    return ar_n_o

def from_SIRI_to_soil_properties(file_bin_SIRI,file_table_SIRI_soil):
    """
    * Objective:

      Extraction of the parameters L (soil depth) and theta_s
      (porosity or humidity at saturation) for each catchment cell
      from the SIRI map

    * Input

      - file_bin_SIRI is the binary grid file containing the SIRI soil
        property codes

      - file_table_SIRI_soil is an ASCII file containing a table of
      correspondance between the SIRI codes and the values of a
      selection of the soil parameters proposed by SIRI (Land-type,
      depth A, depth B, WP A, WP B, FC A, FC B, Por A, Por B)

    * Ouput

      This routine returns two 1D array (ar_L, ar_theta_s) containing
      respectively the values of L and theta_s for each cell. Cells
      are ordered from West to East, North to South.

    * Comment:

      This routine has to be used carefully since the code is
      dependant to:

       1. the format of the Table (file_table_SIRI_soil)

       2. the choice of the user to use only one soil layer A, or
       combining the two soil layer A+B characteristics provided by
       SIRI. Here L=LA+LB and Theta_S=average(theta_sA,theta_sB)

    """
    #Read the binary grid file of GLCC land use type
    tab=read_arc_bin(file_bin_SIRI)
    nrows=np.shape(tab)[0]
    ncols=np.shape(tab)[1]
    tab=np.reshape(tab,ncols*nrows)
    ind=np.where(tab>-99)
    ar_SIRI=tab[ind]

    #Read the Table file within a header line
    tab=pm.read_column_input(file_table_SIRI_soil,9)
    ar_code=tab[:,0]
    ar_depthA=tab[:,1]
    ar_depthB=tab[:,2]
    ar_porA=tab[:,7]
    ar_porB=tab[:,8]

    ar_L=np.array(ar_SIRI)
    ar_theta_s=np.array(ar_SIRI)
    for i in ar_code:
        ind=np.where(ar_L==i)
        #!!!! TO BE CHANGED ACCORDING TO USER CHOICE !!!#
        ar_L[ind]=ar_depthA[np.where(ar_code==i)][0]+ar_depthB[np.where(ar_code==i)][0]
        ar_theta_s[ind]=0.5*(ar_porA[np.where(ar_code==i)][0]+ar_porB[np.where(ar_code==i)][0])

    return ar_L, ar_theta_s

def from_WRC90_to_soil_properties(file_bin_WRC90, file_table_WRC90_soil):
    """
    * Objective:

      Extraction of the parameters L (soil depth) and theta_s
      (porosity or humidity at saturation) for each catchment cell
      from the SIRI map

    * Input

      - file_bin_WRC90 is the binary grid file containing the WRC90
        soil property codes (Here only three 3 for Loamy Sand, 2 for
        Sandy Loam, 1 for Clay)

      - file_table_WRC90_soil is an ASCII file containing a table of
        correspondance between the WRC90 codes and the values of Ks
        (permeability) and theta_r (residual soil moisture)

    * Ouput

      This routine returns two 1D array (ar_theta_r, ar_theta_s)
      containing respectively the values of Ks and theta_r for each
      cell. Cells are ordered from West to East, North to South.

    """
    #Read the binary grid file of GLCC land use type
    tab=read_arc_bin(file_bin_WRC90)
    nrows=np.shape(tab)[0]
    ncols=np.shape(tab)[1]
    tab=np.reshape(tab,ncols*nrows)
    ind=np.where(tab>-99)
    ar_WRC90=tab[ind]

    #Read the Table file within a header line
    tab=pm.read_column_input(file_table_WRC90_soil,9)
    ar_code=tab[:,0]
    ar_theta_r_moy=tab[:,3]
    ar_theta_r_ect=tab[:,4]
    ar_conduct=tab[:,8]

    ar_theta_r=np.array(ar_WRC90)
    ar_Ks=np.array(ar_WRC90)
    for i in ar_code:
        ind=np.where(ar_WRC90==i)
        #!!!! TO BE CHANGED FOR PARAMETER ADJUSTMENT !!!#
        ar_theta_r[ind]=ar_theta_r_moy[np.where(ar_code==i)][0]
        ar_Ks[ind]=ar_conduct[np.where(ar_code==i)][0]

    return ar_theta_r, ar_Ks

def from_Strahler_to_channel_manning(file_bin_strahler,
                                     file_table_strahler_manning, ar_lambda):
    """
    * Objective:

      Extraction of the channel mannings for each channel cell within
      the catchment from the strahler order map

    * Input

      - file_bin_strahler is the binary grid file containing the
        strahler order of the channel cells

      - file_table_strahler_manning is an ASCII file containing a
      table of correspondance between the strahler order and the
      values of manning strickler proposed by Liu and Todini (2002)

      - ar_lambda is an array (dimension equal to the number of cell)
        with value 1 for channel cells, 0 otherwise. Cells being
        ordered from West to East, North to South.

    * Ouput

      This routine returns a 1D array (ar_n_c), ar_theta_s) containing
      the values of the manning coefficient for each channel
      cell. Cells are ordered from West to East, North to South.

    """

    #Read the binary grid file of GLCC land use type
    tab=read_arc_bin(file_bin_strahler)
    nrows=np.shape(tab)[0]
    ncols=np.shape(tab)[1]
    tab=np.reshape(tab,ncols*nrows)
    ind=np.where(tab>-99)
    ar_strahler=tab[ind]
    ar_lambda[ar_lambda==1]=ar_strahler

    #Read the Table file within a header line
    tab=pm.read_column_input(file_table_strahler_manning,2)
    ar_code=tab[:,0]
    ar_manning=tab[:,1]

    ar_n_c=np.array(ar_lambda)
    for i in ar_code:
        ind=np.where(ar_n_c==i)
        #!!!! TO BE CHANGED FOR PARAMETER ADJUSTMENT !!!#
        ar_n_c[ind]=ar_manning[np.where(ar_code==i)][0]

    return ar_n_c


def from_flowdir_to_celldown_4D(file_bin_flowdir):
    """
    * objective:

      Associate to each cell the label of the celldown (outcell)
      according to the 4D flow directions

    * Input

      - file_bin_flowdir: 4D flow direction binary file (from GIS)
        with codes: 1 right, 16 left, 64 up, 4 down

    * Output

      - ar_cell_down: an array (dimension equal to the number of cell)
        containing the label of the outcells.

    * Commment

     Some errors can occur in the direction files that should be
     detected by the routine. A manual correction is required if it
     does happen.

    """

    tab_flowdir=read_arc_bin(file_bin_flowdir)
    tab_label=from_bingrid_to_label(file_bin_flowdir)
    nrows=np.shape(tab_label)[0]
    ncols=np.shape(tab_label)[1]
    ar_cell_down=np.zeros(np.shape(np.where(tab_label>=0))[1])-99
    n=0
    num=0
    outlet=0
    for i in range(nrows):
        for j in range(ncols):
            flow=tab_flowdir[i,j]
            OK=0
            if tab_label[i,j] >= 0:
                num=num+1
                if flow==1:
                    x=i
                    y=j+1
                elif flow==16:
                    x=i
                    y=j-1
                elif flow==64:
                    x=i-1
                    y=j
                elif flow==4:
                    x=i+1
                    y=j
                else:
                    print 'ERROR'
                if tab_label[x,y]>=0:
                    ar_cell_down[tab_label[i,j]]=tab_label[x,y]
                else:
                    if outlet==0:
                        ar_cell_down[tab_label[i,j]]=tab_label[x,y]
                        outlet=1
                    else:
                        n=n+1
                        #Print the data to detect where the errors are
                        #(manual correction to be done...)
                        print n,num, flow,i,j,tab_label[i,j],x,y,tab_label[x,y]
    return ar_cell_down

def cell_connectivity(flowdir, mask, source='GRASS'):
    """Compute the connectivity between cells in the catchment

    Associate each cell in the catchment with the label of it's
    downstream neighbour. This defines the directed network of
    connections between the cells which make up a catchment in the
    model.

    Parameters
    ----------
    flowdir : Numpy ndarray
        8D flow direction codes as defined by your favourite GIS
        toolbox.
    mask : Numpy ndarray
        A 2D array with the same shape as `flowdir`. Cells which
        comprise the catchment should be marked by a value of 1.
    source : string
        A string describing the source of the flow direction
        codes. Current options are 'ArcGIS' or 'GRASS' (default)

    Returns
    -------
    cell_down : Numpy ndarray
        An ordered array containing the label of the immediate
        downstream cell for each cell in the catchment network. A 1D
        array with length equal to the number of cells. The catchment
        outlet is indiacted by a value of -999.

    """
    # Specify flow direction code from GRASS GIS r.watershed or ArcGIS
    # Hydrology toolbox flow-direction tool
    if source == 'GRASS':
        ddict = {1 : (-1,  1),
                 2 : (-1,  0),
                 3 : (-1, -1),
                 4 : ( 0, -1),
                 5 : ( 1, -1),
                 6 : ( 1,  0),
                 7 : ( 1,  1),
                 8 : ( 0,  1)}
    elif source == 'ArcGIS':
        ddict = {128 : (-1,  1),
                 64  : (-1,  0),
                 32  : (-1, -1),
                 16  : ( 0, -1),
                 8   : ( 1, -1),
                 4   : ( 1,  0),
                 2   : ( 1,  1),
                 1   : ( 0,  1)}
    else:
        raise ValueError('Unknown flow direction scheme: %s' % source)

    ncells = mask[mask == 1].size
    int_min = np.iinfo(np.int).min
    cell_id = np.ones(mask.shape, dtype=np.int)*int_min
    cell_id[mask == 1] = np.arange(ncells)

    cell_down = np.ones(ncells, dtype=np.int)*int_min

    nrows, ncols = mask.shape
    for i in range(nrows):
        for j in range(ncols):
            fdir = flowdir[i, j]

            if fdir in ddict.keys():
                r, c = ddict[fdir]
                m, n = i+r, j+c

                # To-do: Handle case where (m, n) is outside the array
                # bounds
                cell_down[cell_id[i, j]] = cell_id[m, n]

    if cell_down[cell_down == int_min].size > 1:
        warn_txt = """There are %d catchment cells without a downstream link.
Check the validity of the flow-direction raster."""  \
        % cell_down[cell_down == int_min].size

        warn(warn_txt)

    cell_down[cell_down == int_min] = -999

    return cell_down

def from_bingrid_to_label(file_bin_grid,write_file=False):
    """
    * Objective

      Replace values of the grid by the label of the cells.  The label
      are assigned from 0 to nb_cell, from West to East, North to
      South.

    * Input

      file_bin_grid: A binary grid file (whatever it is).

    * Output

      tab: a 1D array with nb_cell components.

    """
    tab=read_arc_bin(file_bin_grid)
    nrows=np.shape(tab)[0]
    ncols=np.shape(tab)[1]
    tab=np.reshape(tab,ncols*nrows)
    ind=np.where(tab>-99)
    tab_label=np.arange(len(tab[ind]))
    for i in tab_label:
        tab[ind[0][i]]=i
    tab=np.reshape(tab,(nrows,ncols))
    tab=tab.astype('int32')

    if write_file:
        f = file('c:/Theo/liebenbergsvlei/topkapi_model/parameters/label.dat', 'w')
        np.savetxt(f, tab)
        f.close()

    return tab

def from_bingrid_to_coordinate(file_bin_grid):
    """
    * Objective

    Compute the coordinates of the cells relatively to the low-left
    corner from any binary grid file

    * Input

      file_bin_grid: A binary grid file (whatever it is).

    * Output

      ar_coorX,ar_coorY: two 1D array with nb_cell components.

    """
    li_headers=read_headers_arc_bin(file_bin_grid)
    ncols = li_headers[0] # fixed for UM grid
    nrows = li_headers[1] # fixed for UM grid
    xmin = li_headers[2]
    ymin = li_headers[3]
    cellsize= li_headers[4]

    ar_line_coorX=xmin+np.arange(ncols)*cellsize
    ar_line_coorY=ymin+np.arange(nrows)*cellsize
    ar_line_coorY=ar_line_coorY[::-1]

    mat_coorX=np.zeros((nrows,ncols))
    mat_coorY=np.zeros((nrows,ncols))
    for i in range(int(nrows)):
        mat_coorX[i,:]=ar_line_coorX
    for j in range(int(ncols)):
        mat_coorY[:,j]=ar_line_coorY

    tab=read_arc_bin(file_bin_grid)
    ind=np.where(tab>-99.)

    ar_coorX=mat_coorX[ind]
    ar_coorY=mat_coorY[ind]

    return ar_coorX,ar_coorY


def compute_Xchannel(cell_labels, channel_network, X, Y, cell_down):
    """Compute the channel length in channel cells

    Cells draining diagonally have a different channel length from
    cells draining North, South, East or West. This function computes
    the channel length as a function of the drainage direction (based
    on the catchment connectivity network).

    Parameters
    ----------
    cell_labels : 1D Numpy ndarray
        An array of the labels associated with each cell in the
        catchment
    channel_network : 1D Numpy ndarray
        An ordered array with each channel cell indicated by a value
        of one.
    X : 1D Numpy ndarray
        An ordered array of the X coordinate of the centre of each
        cell.
    Y : 1D Numpy ndarray
        An ordered array of the Y coordinate of the centre of each
        cell.
    cell_down : 1D Numpy ndarray
        An ordered array giving the label of the downstream cell in
        the catchment network. The outlet of the catchment is
        indicated by a negative number.

    Returns
    -------
    Xc : 1D Numpy ndarray
        An array containing the length of the channel in each channel
        cell, zero otherwise.

    """
    ar_Xc = np.zeros(cell_labels.shape)

    for i in cell_labels[channel_network == 1]:
        if cell_down[i]>=0:
            indx = cell_down[i]
            Xcell = X[i]
            Ycell = Y[i]

            Xcell_down = X[cell_labels == indx]
            Ycell_down = Y[cell_labels == indx]

            ar_Xc[i] = ut.distance(Xcell, Ycell, Xcell_down, Ycell_down)

    ind_outlet=np.where(cell_down<0)
    ar_Xc[ind_outlet]=min(ar_Xc)

    return ar_Xc

def matrix_plot(matrix, fig_name, title='GRID Plot'):
    """Create a plot of the data in a GRIB1 file."""

    a=matrix

    a_mask = ma.masked_where(a < 0, a)
    pl.imshow(a_mask, interpolation='nearest')
    pl.colorbar()
    pl.title(title)
    pl.savefig(fig_name)
    pl.close()

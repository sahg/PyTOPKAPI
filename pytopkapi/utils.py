import os
import os.path
from subprocess import Popen, PIPE

import h5py
import numpy as np

import pytopkapi

def show_banner(ini_file, nb_cell, nb_time_step):
    """Show an ASCII banner at run time describing the model.

    Parameters
    ----------
    ini_file : str
        The name of the PyTOPKAPI initialization file passed to the
        run function.
    nb_cell : int
        The number of model cells to be processed.
    nb_time_step : int
        The number of model time-steps to be processed.

    """

    print('===============================================================\n',
          '\n PyTOPKAPI\n',
          'A Python implementation of the TOPKAPI Hydrological model\n\n',
          'Version {}\n'.format(pytopkapi.__version__),
          'Number of model cells: {:d}\n'.format(nb_cell),
          'Number of model time-steps: {:d}\n'.format(nb_time_step),
          'Running simulation from file: {}\n'.format(ini_file),
          '\r===============================================================\n')

def _create_dataset(h5file, grp_name, dset_name, shape, units):
    """Create HDF5 dataset if it doesn't exist.

    """
    if '{}/{}'.format(grp_name, dset_name) not in h5file:
        dset = h5file.create_dataset('{}/{}'.format(grp_name, dset_name),
                               shape, maxshape=(None, None), compression='gzip')

        dset.attrs['units'] = units

def open_simulation_file(file_out, fmode, Vs0, Vo0, Vc0, no_data,
                         nb_cell, nb_time_step, append_output, first_run):
    """Open simulation file and return handles to it's content.

    """
    h5file = h5py.File(file_out, fmode)

    dset_shape = (nb_time_step+1, nb_cell)

    h5file.attrs['title'] = 'PyTOPKAPI simulation'
    h5file.attrs['pytopkapi_version'] = pytopkapi.__version__
    h5file.attrs['pytopkapi_git_revision'] = pytopkapi.__git_revision__

    # create file structure as necessary
    grp_name = 'Soil'

    dset_name = 'Qs_out'
    _create_dataset(h5file, grp_name, dset_name, dset_shape, 'm3/s')
    dset_Qs_out = h5file['{}/{}'.format(grp_name, dset_name)]

    dset_name = 'V_s'
    _create_dataset(h5file, grp_name, dset_name, dset_shape, 'm3')
    dset_Vs = h5file['{}/{}'.format(grp_name, dset_name)]

    grp_name = 'Overland'

    dset_name = 'Qo_out'
    _create_dataset(h5file, grp_name, dset_name, dset_shape, 'm3/s')
    dset_Qo_out = h5file['{}/{}'.format(grp_name, dset_name)]

    dset_name = 'V_o'
    _create_dataset(h5file, grp_name, dset_name, dset_shape, 'm3')
    dset_Vo = h5file['{}/{}'.format(grp_name, dset_name)]

    grp_name = 'Channel'

    dset_name = 'Qc_out'
    _create_dataset(h5file, grp_name, dset_name, dset_shape, 'm3/s')
    dset_Qc_out = h5file['{}/{}'.format(grp_name, dset_name)]

    dset_name = 'V_c'
    _create_dataset(h5file, grp_name, dset_name, dset_shape, 'm3')
    dset_Vc = h5file['{}/{}'.format(grp_name, dset_name)]

    dset_name = 'Ec_out'
    _create_dataset(h5file, grp_name, dset_name, dset_shape, 'm3')
    dset_Ec_out = h5file['{}/{}'.format(grp_name, dset_name)]

    grp_name = ''

    dset_name = 'ET_out'
    _create_dataset(h5file, grp_name, dset_name, dset_shape, 'mm')
    dset_ET_out = h5file['{}/{}'.format(grp_name, dset_name)]

    dset_name = 'Q_down'
    _create_dataset(h5file, grp_name, dset_name, dset_shape, 'm3/s')
    dset_Q_down = h5file['{}/{}'.format(grp_name, dset_name)]

    if append_output is False or first_run is True:
        #Write the initial values into the output file
        dset_Vs[0] = Vs0
        dset_Vo[0] = Vo0
        dset_Vc[0] = Vc0

        dset_Qs_out[0] = np.ones(nb_cell)*no_data
        dset_Qo_out[0] = np.ones(nb_cell)*no_data
        dset_Qc_out[0] = np.zeros(nb_cell)

        dset_Q_down[0] = np.ones(nb_cell)*no_data

        dset_ET_out[0] = np.zeros(nb_cell)
        dset_Ec_out[0] = np.zeros(nb_cell)

    return h5file, dset_Vs, dset_Vo, dset_Vc, dset_Qs_out, \
           dset_Qo_out, dset_Qc_out, dset_Q_down, dset_ET_out, dset_Ec_out

# System utility functions
def exec_command(cmd_args):
    """Execute a shell command in a subprocess

    Convenience wrapper around subprocess to execute a shell command
    and pass back stdout, stderr, and the return code. This function
    waits for the subprocess to complete, before returning.

    Usage example:
    >>> stdout, stderr, retcode = exec_command(['ls', '-lhot'])

    Parameters
    ----------
    cmd_args : list of strings
        The args to pass to subprocess. The first arg is the program
        name.

    Returns
    -------
    stdout : string
        The contents of stdout produced by the shell command
    stderr : string
        The contents of stderr produced by the shell command
    retcode : int
        The return code produced by the shell command

    """
    proc = Popen(cmd_args, stdout=PIPE, stderr=PIPE)
    stdout, stderr = proc.communicate()
    proc.wait()

    return stdout, stderr, proc.returncode

########################
##   For graphics     ##
########################
def CRange(ar_x):
    '''
    Returns the range of an array
    '''
    lim=np.array([min(ar_x),max(ar_x)])
    return lim

def f_axe(p,xc):
    '''
    Returns the value in the array xc
    associated to a relative value inside [0,1]
    '''
    xc=np.sort(xc)
    pos=xc[0]+p*(xc[-1]-xc[0])
    return pos

def string(integer,len_str_out):
    """
    From a given integer, return an string of length len_str_out completed by zero
    Example:
    ut.string(1,3)-->'001'
    """
    str_zero='0'
    str_int=str(integer)
    len_int=len(str_int)
    if len_str_out-len_int<0:
        print('****ERROR: length of string too short')
        str_out=''
    else:
        str_out=(len_str_out-len_int)*str_zero+str_int
    return str_out

def from_float_array_to_string_array(ar_float,unique=False):
    if unique:
        a=str(np.unique(ar_float)).split()[1:]
    else:
        a=str(ar_float).split()[1:]
    a[-1]=a[-1].replace(']','')
    ar_string=a
    return ar_string

##############################
##   For file management    ##
##############################

def check_file_exist(filename):
    path_name, file_name = os.path.split(filename)
    if not os.path.exists(path_name) and path_name != '':
        os.makedirs(path_name)

def check_folder_exist(folder_name):
    if not os.path.exists(folder_name):
        os.mkdir(folder_name)

def read_one_array_hdf(file_h5, group, name):
    """Read a single array from a PyTOPKAPI simulation file.

    """
    h5file_in = h5py.File(file_h5)

    dset_string = '/%s/%s' % (group, name)
    array = h5file_in[dset_string][...]

    h5file_in.close()

    return array

##############################
##        Statistics        ##
##############################

def mov_avg(ar_float,period):
    '''
    period is a multiple of 2
    '''

    nb_ind=len(ar_float)-period
    ar_out=np.zeros(nb_ind)
    for i in range(nb_ind):
        n=period/2
        ind_mid=i+n
        ar_out[i]=np.average(ar_float[ind_mid-n:ind_mid+n])

    return ar_out

##~~~   Comparison of 2 vectors   ~~~~##
# Functions defining useful criteria comparing two vectors
## REFERENCE is Y
def R(ar_x,ar_y):
    R=np.corrcoef(ar_x,ar_y)
    return R[0,1]

def R2(ar_x,ar_y):
    R=np.corrcoef(ar_x,ar_y)
    return R[0,1]**2

def Nash(ar_x,ar_y):
    eff=1-sum((ar_y-ar_x)**2)/sum((ar_y-np.mean(ar_y))**2)
    return eff

def RMSE(ar_x,ar_y):
    rmserr=(np.mean((ar_y-ar_x)**2))**0.5
    return rmserr

def RMSE_norm(ar_x,ar_y):
    rmserr=(np.mean((ar_y-ar_x)**2))**0.5
    rmsenorm=rmserr/np.mean(ar_y)
    return rmsenorm

def Bias_cumul(ar_x,ar_y):
    b=sum(ar_x)/sum(ar_y)
    return b

def Diff_cumul(ar_x,ar_y):
    diff=sum(ar_x)-sum(ar_y)
    return diff

def Abs_cumul(ar_x,ar_y):
    abs_diff=abs(sum(ar_x)-sum(ar_y))
    return abs_diff

def Err_cumul(ar_x,ar_y):
    err_rel=abs(sum(ar_x)-sum(ar_y))/sum(ar_y)
    return err_rel


##############################
##        HDF5 files        ##
##############################

####HOW to remove a group
#h5file.removeNode('/', groupname)

##############################
##      Works on vectors    ##
##############################

def find_dist_max(ar_coorx,ar_coory):
    """
    Compute the maximum distance between several points defined by their coordinates ar_coorx and ar_coory
    """
    nb_cell=len(ar_coorx)
    max_dist=0.
    for i in range(nb_cell):
        for j in range(nb_cell):
            max_dist=max(max_dist,distance(ar_coorx[i],ar_coory[i],ar_coorx[j],ar_coory[j]))
    return max_dist

def distance(x1,y1,x2,y2):
    """
    Compute the distance between two points
    """
    dist=((x1-x2)**2+(y1-y2)**2)**0.5
    return dist

def find_cell_coordinates(ar_cell_label, Xoutlet, Youtlet,
                          ar_coorx, ar_coory, ar_lambda, channel=True):
    """Find the label of the cell closest to (Xoutlet, Youtlet).

    Find the label of the model cell containing the specified location. The
    co-ordinates of the location must be given in the same co-ordinate system
    as that specifying the model catchment.

    Parameters
    ----------
    ar_cell_label : (N,) int array
        Numbers labelling each cell.
    Xoutlet : float
        The x co-ordinate of a point. This is the Longitude expressed in
        metres using the same projection as `ar_coorx`.
    Youtlet : float
        The y co-ordinate of a point. This is the Longitude expressed in
        metres using the same projection as `ar_coory`.
    ar_coorx : (N,) float array
        The x co-ordinate of the centre of each cell (m). This is the Longitude
        expressed in metres using a Transverse Mercator projection, but any
        appropriate projection can be used.
    ar_coory : (N,) float array
        The y co-ordinate of the centre of each cell (m). This is the Latitude
        expressed in metres using a Transverse Mercator projection, but any
        appropriate projection can be used.
    ar_lambda : (N,) int array
        Switch indicating whether the current cell contains a channel. A value
        of `1` indicates a channel cell, `0` indicates no channel.
    channel : boolean (default=True)
        Allows cells with or without channels to be chosen.

    Returns
    -------
    cell_outlet : int
        The label for the cell closest to the defined location.

    """
    tab_x=np.unique(ar_coorx);X=abs(tab_x[0]-tab_x[1])
    dist_max=3*X
    dist_min=dist_max
    nb_cell=len(ar_cell_label)
    cell_outlet=-999.9
    for i in range(nb_cell):
        dist=distance(Xoutlet,Youtlet,ar_coorx[i],ar_coory[i])
        if channel:
            if dist < dist_min and ar_lambda[i]==1.:
                dist_min=dist
                cell_outlet=ar_cell_label[i]
        else:
            if dist<dist_min:
                dist_min=dist
                cell_outlet=ar_cell_label[i]


    if cell_outlet<0:
        print("Wrong coordinates")
        stop
    return cell_outlet

import os
import os.path
import tables as h5
import numpy as np

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
        print '****ERROR: length of string too short'
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
        print path_name, 'has been created'
        os.makedirs(path_name)

def check_folder_exist(folder_name):
    if not os.path.exists(folder_name):
        print folder_name, 'has been created'
        os.mkdir(folder_name)

def read_one_array_hdf(file_h5,group,name):
    h5file_in=h5.openFile(file_h5,mode='r')
    node = h5file_in.getNode(group+name)
    array=node.read()
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
        print "Wrong coordinates"
        stop
    return cell_outlet

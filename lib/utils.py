import os
import os.path
import tables as h5
import numpy as np

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

def check_file_exist(filename):
    folder_name=os.path.split(filename)[0]
    if not os.path.exists(folder_name):
        os.mkdir(folder_name)

def check_folder_exist(folder_name):
    if not os.path.exists(folder_name):
        os.mkdir(folder_name)
        
def read_one_array_hdf(file_h5,group,name):
    h5file_in=h5.openFile(file_h5,mode='r')
    node = h5file_in.getNode(group+name)
    array=node.read()
    h5file_in.close()
    return array

def from_float_array_to_string_array(ar_float,unique=False):
    if unique:
        a=str(np.unique(ar_float)).split()[1:]
    else:
        a=str(ar_float).split()[1:]
    a[-1]=a[-1].replace(']','')
    ar_string=a
    return ar_string

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
####HOW to remove a group
#h5file.removeNode('/', groupname)

def find_cell_coordinates(ar_cell_label,Xoutlet,Youtlet,ar_coorx,ar_coory,ar_lambda,channel=True):
    """
    Find the label of the closest cell from (Xoutlet, Youtlet)
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

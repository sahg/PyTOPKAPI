""" pretreatment_module.py

Functions required to compute the intrinsic TOPKAPI parameters from the
physical parameters
"""

__author__ = "Theo Vischel"
__version__ = "$Revision: 1.0 $"
__date__ = "$Date: 09/10/2006 $"

import scipy as sp
import numpy as np
import scipy.io as io


#```````````````````````````````````````````
def read_global_parameters(file_name):
    """ read_global_parameters
        Read the file containing the general parameters of the model
        similar for all cells
        Return:
           X,Dt,alpha_s,alpha_o,alpha_c,nb_cell
    """
    file_read=open(file_name,'r')
    ## reading lines
    tab_read=file_read.readlines()
    for line in tab_read[1:]:
        ##split the blanks
        donnees=line.split()
        ##append the values to variable lists
        X=float(donnees[0])
        Dt=float(donnees[1])
        alpha_s=float(donnees[2])
        alpha_o=float(donnees[3])
        alpha_c=float(donnees[4])
        A_thres=float(donnees[5])
        W_min=float(donnees[6])
        W_max=float(donnees[7])
        ##end of iteration the file is closed
        file_read.close()
        
    return X,Dt,alpha_s,alpha_o,alpha_c,A_thres,W_min,W_max
            
#`````````````````````````````````````````` 
def read_cell_parameters(file_name):
    """ read_cell_parameters
        Read the file containing the phsical parameters of each cell
    """
    tab_read=io.read_array(file_name)

    ar_cell_label=np.array(tab_read[:,0],int)
    ar_coorx=tab_read[:,1]
    ar_coory=tab_read[:,2]
    ar_lambda=np.array(tab_read[:,3],int)
    ar_Xc=np.array(tab_read[:,4])
    ar_dam=np.array(tab_read[:,5],int)
    ar_tan_beta=tab_read[:,6]
    ar_tan_beta_channel=tab_read[:,7]
    ar_L=tab_read[:,8]
    ar_Ks=tab_read[:,9]
    ar_theta_r=tab_read[:,10]
    ar_theta_s=tab_read[:,11]
    ar_n_o=tab_read[:,12]
    ar_n_c=tab_read[:,13]
    ar_cell_down=np.array(tab_read[:,14],int)
    ar_pVs_t0=tab_read[:,15]
    ar_Vo_t0=tab_read[:,16]
    ar_Qc_t0=tab_read[:,17]
    ar_kc=tab_read[:,18]

    return ar_cell_label,ar_coorx,ar_coory,ar_lambda,ar_Xc,ar_dam,ar_tan_beta,ar_tan_beta_channel,ar_L,ar_Ks,\
           ar_theta_r,ar_theta_s,ar_n_o,ar_n_c,\
           ar_cell_down,ar_pVs_t0,ar_Vo_t0,ar_Qc_t0,ar_kc

#``````````````````````````````````````````   
def sort_cell(ar_cell_label,ar_cell_down):
    """ sort_cell
        Sort the cell in a computation order required by the distance to outlet
        Return an array containing the sorted cell labels
    """
    ar_label_sort=sp.ones(len(ar_cell_label))*-99.9
    ar_dist_2_outlet=sp.ones(len(ar_cell_label))*-99.9
    ar_drained_area=sp.ones(len(ar_cell_label))*-99.9
    nb_cell=len(ar_cell_label)
    for cell in range(nb_cell):
        cell_down=ar_cell_down[cell]
        dist=0
        while cell_down > -99:
            cell_down=ar_cell_down[cell_down]
            dist=dist+1
        ar_dist_2_outlet[cell]=dist
    a=np.argsort(ar_dist_2_outlet)
    ar_label_sort=a[::-1]
    ar_label_sort=np.array(ar_label_sort,int)
    
    return ar_label_sort

#``````````````````````````````````````````   
def direct_up_cell(ar_cell_label,ar_cell_down,ar_label_sort):
    """ cell_up
        Give the direct upcells for each cell
    """
    nb_cell=len(ar_label_sort)
    li_cell_up=[]
    for ncell in range(len(ar_label_sort)):
        cell_label=ar_cell_label[ncell]
        up_cell=ar_cell_label[np.where(ar_cell_down==cell_label)]
        li_cell_up.append(up_cell)
    return li_cell_up

#``````````````````````````````````````````   
def drained_area(ar_label_sort,li_cell_up,X):
    """ drained_area
        Compute the drained area for each cell
    """
    A_cell=X**2
    nb_cell=len(ar_label_sort)
    ar_A_drained=sp.ones(nb_cell)*-99.9
    for cell in ar_label_sort:
        up_cell=li_cell_up[cell]
        ar_A_drained[cell]=sum(ar_A_drained[up_cell])+A_cell
    return ar_A_drained

#```````````````````````````````````````````````````````````
def compute_cell_param(X,ar_Xc,Dt,alpha_s,alpha_o,alpha_c,nb_cell,\
                       A_thres,W_max,W_min,\
                       ar_lambda,ar_tan_beta,ar_tan_beta_channel,ar_L,\
                       ar_Ks,ar_theta_r,ar_theta_s,ar_n_o,ar_n_c,\
                       ar_A_drained):
    """compute_cell_param
        Compute: For all cells
                 1. Soil cell parameters
                 - the saturated moisture volume   -->Vsm
                 - the constant term of the non differential equation
                   dV_s/dt=a_s-b_s*V_s**alpha_s    -->b_s
                 2. Overland cell parameters
                 - the constant term of the non differential equation
                   dV_o/dt=a_o-b_o*V_o**alpha_o    -->b_o
                 3. Channel cell parameters
                 - The channel width W             -->W
                 - The constant term of the non differential equation
                   dV_c/dt=a_c-b_c*V_c**alpha_c    -->b_c
        Return: ar_Vsm, ar_b_s, ar_b_o, ar_W, ar_b_c
        """
    ##Soil parameters    
    ar_Vsm=(ar_theta_s-ar_theta_r)*(X**2)*ar_L
    ar_Cs=(ar_L*ar_Ks*ar_tan_beta)/(((ar_theta_s-ar_theta_r)*ar_L)**alpha_s)
    ar_b_s=ar_Cs*X/(X**(2*alpha_s))
    ##Overland parameters
    ar_Co=(1/ar_n_o)*(ar_tan_beta)**0.5
    ar_b_o=ar_Co*X/(X**(2*alpha_o))
    ##Channel parameters
    A_total=nb_cell*X**2
    ar_W=W_max+((W_max-W_min)/(A_total**0.5-A_thres**0.5))*(ar_A_drained**0.5-A_total**0.5)
    ar_Cc=(1/ar_n_c)*(ar_tan_beta_channel)**0.5
    ar_b_c=ar_Cc*ar_W/((ar_Xc*ar_W)**(alpha_c))
    ar_W[ar_lambda==0]=-99.9
    ar_b_c[ar_lambda==0]=-99.9
    
    return ar_Vsm, ar_b_s, ar_b_o, ar_W, ar_b_c

#``````````````````````````````````````````   
def read_column_input(file_name,nb_cell):
    """ read_column_input
        Read the file containing data in column format:
        Cell1  Cell2  Cell3 ...
        1.3    4.3     5.2  ...
        2.3    5.6     4.2  ...

        Return a matrix mat_out(nrow, ncol) 
    """
    file_read=open(file_name,'r')
    tab_read=file_read.readlines()

    nb_time_step=len(tab_read)-1
    a=sp.zeros(nb_cell*nb_time_step)
    mat_out=a.reshape(nb_time_step,nb_cell)

    i=-1
    for line in tab_read[1:]:
        i=i+1
        donnees=line.split()
        mat_out[i,]=[float(elem) for elem in donnees]
    ##end of iteration the file is closed
    file_read.close()
    
    return mat_out

    

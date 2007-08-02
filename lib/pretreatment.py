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
def read_global_parameters_old(file_name):
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
        nb_cell=int(donnees[5])
        A_total=float(donnees[6])
        A_thres=float(donnees[7])
        W_min=float(donnees[8])
        W_max=float(donnees[9])
        ##end of iteration the file is closed
        file_read.close()
        
    return X,Dt,alpha_s,alpha_o,alpha_c,nb_cell,A_total,A_thres,W_min,W_max

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
def read_cell_parameters(file_name,begin=0):
    """ read_cell_parameters
        Read the file containing the phsical parameters of each cell
        Return:
           ar_num_cell,ar_lambda,ar_tan_beta,ar_L,ar_Ks,
           ar_theta_r,ar_theta_s,ar_n_o,ar_n_c,ar_A_drained,
           ar_cell_down,ar_pVs_t0,ar_pVc_t0
    """
    file_read=open(file_name,'r')
    ## variable list declaration
    li_num_cell=[]
    li_coorX=[]
    li_coorY=[]
    li_lambda=[]
    li_tan_beta=[]
    li_L=[]
    li_Ks=[]
    li_theta_r=[]
    li_theta_s=[]
    li_n_o=[]
    li_n_c=[]
    li_cell_down=[]
    li_pVs_t0=[]
    li_pVc_t0=[]
    li_kc=[]
    ## reading lines
    tab_read=file_read.readlines()
    for line in tab_read[begin:]:
        ##split the blanks
        donnees=line.split()
        ##append the values to variable lists
        li_num_cell.append(int(float(donnees[0])))
        li_coorX.append(float(donnees[1]))
        li_coorY.append(float(donnees[2]))
        li_lambda.append(int(float(donnees[3])))
        li_tan_beta.append(float(donnees[4]))
        li_L.append(float(donnees[5]))
        li_Ks.append(float(donnees[6]))
        li_theta_r.append(float(donnees[7]))
        li_theta_s.append(float(donnees[8]))
        li_n_o.append(float(donnees[9]))
        li_n_c.append(float(donnees[10]))
        li_cell_down.append(int(float(donnees[11])))
        li_pVs_t0.append(float(donnees[12]))
        li_pVc_t0.append(float(donnees[13]))
        li_kc.append(float(donnees[14]))
    ##end of iteration the file is closed
    file_read.close()

    ar_num_cell=sp.array(li_num_cell)
    ar_coorx=sp.array(li_coorX)
    ar_coory=sp.array(li_coorY)
    ar_lambda=sp.array(li_lambda)
    ar_tan_beta=sp.array(li_tan_beta)
    ar_L=sp.array(li_L)
    ar_Ks=sp.array(li_Ks)
    ar_theta_r=sp.array(li_theta_r)
    ar_theta_s=sp.array(li_theta_s)
    ar_n_o=sp.array(li_n_o)
    ar_n_c=sp.array(li_n_c)
    ar_cell_down=sp.array(li_cell_down)
    ar_pVs_t0=sp.array(li_pVs_t0)
    ar_pVc_t0=sp.array(li_pVc_t0)
    ar_kc=sp.array(li_kc)
    return ar_num_cell,ar_coorx,ar_coory,ar_lambda,ar_tan_beta,ar_L,ar_Ks,\
           ar_theta_r,ar_theta_s,ar_n_o,ar_n_c,\
           ar_cell_down,ar_pVs_t0,ar_pVc_t0,ar_kc

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
def compute_cell_param(X,Dt,alpha_s,alpha_o,alpha_c,nb_cell,\
                       A_thres,W_max,W_min,\
                       ar_lambda,ar_tan_beta,ar_L,\
                       ar_Ks,ar_theta_r,ar_theta_s,ar_n_o,ar_n_c,\
                       ar_A_drained,ar_pVs_t0,ar_pVc_t0):
    """compute_cell_param
        Compute: For all cells
                 1. Soil cell parameters
                 - the saturated moisture volume   -->Vsm
                 - the constant term of the non differential equation
                   dV_s/dt=a_s-b_s*V_s**alpha_s    -->b_s
                 - the initial soil moisture       -->Vs_t0
                 2. Overland cell parameters
                 - the constant term of the non differential equation
                   dV_o/dt=a_o-b_o*V_o**alpha_o    -->b_o
                 3. Channel cell parameters
                 - The channel width W             -->W
                 - The constant term of the non differential equation
                   dV_c/dt=a_c-b_c*V_c**alpha_c    -->b_c
                 - The initial stored moisture     -->Vc_t0
                  N.B= if not a channel cell (i.e. if lambda=0) W=-99.9 b_c=-99.9
        Return: ar_Vsm, ar_b_s, ar_b_o, ar_W, ar_b_c
        """
    ##Soil parameters    
    ar_Vsm=(ar_theta_s-ar_theta_r)*(X**2)*ar_L
    ar_Cs=(ar_L*ar_Ks*ar_tan_beta)/(((ar_theta_s-ar_theta_r)*ar_L)**alpha_s)
    ar_b_s=ar_Cs*X/(X**(2*alpha_s))
    ar_Vs_t0=ar_pVs_t0/100.*ar_Vsm
    ##Overland parameters
    ar_Co=(1/ar_n_o)*(ar_tan_beta)**0.5
    ar_b_o=ar_Co*X/(X**(2*alpha_o))
    ##Channel parameters
    A_total=nb_cell*X**2
    ar_W=W_max+((W_max-W_min)/(A_total**0.5-A_thres**0.5))*(ar_A_drained**0.5-A_total**0.5)
    ar_Cc=(1/ar_n_c)*(ar_tan_beta)**0.5
    ar_b_c=ar_Cc*ar_W/((X*ar_W)**(alpha_c))
    ar_W[ar_lambda==0]=-99.9
    ar_b_c[ar_lambda==0]=-99.9
    ar_Vc_t0=ar_pVc_t0/100.*(X*ar_W**2)
    ar_Vc_t0[ar_lambda==0]=0.
    
    return ar_Vsm, ar_b_s, ar_Vs_t0, ar_b_o, ar_W, ar_b_c, ar_Vc_t0


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

#``````````````````````````````````````````   
def generate_file_rain_ETr_ETo(nb_cell,\
                               ar_hyeto=np.array([5.,10.,20.,50.,50.,30.,20.,10.,5.]),nb_step_zero_rain=90,\
                               ETr_value=5.,\
                               ETo_value=5.):
    """
    Generate simple rainfall and ETr file for testing the model
    """
    ar_rain=np.concatenate((ar_hyeto,np.zeros(nb_step_zero_rain)))
    ndar_rain=np.zeros((len(ar_rain),nb_cell))
    for i in range(nb_cell):
        ndar_rain[:,i]=ar_rain
    ndar_ETr=np.zeros((len(ar_rain),nb_cell))+ETr_value
    ndar_ETo=np.zeros((len(ar_rain),nb_cell))+ETo_value
    
    return ndar_rain,ndar_ETr,ndar_ETo
    

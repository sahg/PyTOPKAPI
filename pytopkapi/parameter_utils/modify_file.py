"""
*** OBJECTIVE

Modifying the created parameter files in order to
   a. extract the parameter for a subcatchment
   b. changing the value of one of the parameters (for instance for calibration or sensitivity analyses)
   - part "MODIFYING THE PARAMETER FILES"
   - main routines

*** COMMENT

"""
import sys

#Internal modules
import numpy as np
import scipy as sp
import pylab as pl
from numpy import ma
import tables as h5
import os
import os.path
from ConfigParser import SafeConfigParser
config = SafeConfigParser()

#External modules
#Utilities
import pytopkapi.utils as ut
#pretreatment: used for subroutines to read the parameter files
import pytopkapi.pretreatment as pm
####################################
### MODIFYING THE PARAMETER FILE ###
####################################

def zero_slope_management(ini_file='zero_slope_management.ini'):

    ##================================##
    ##  Read the input file (*.ini)   ##
    ##================================##
    config.read(ini_file)
    print 'Read the file ',ini_file
    ##~~~~~~~~~~~ input files ~~~~~~~~~~~##
    #Param
    file_cell_param=config.get('input_files','file_cell_param')
    file_cell_param_out=config.get('output_files','file_cell_param_out')
    ##~~~~~~ numerical_values ~~~~~~##
    nb_param=config.getfloat('numerical_values','nb_param')
    X=config.getfloat('numerical_values','X')

    ##============================##
    ## Pretreatment of input data ##
    ##============================##
    print 'Pretreatment of input data'
    #~~~~Read Cell parameters file
    ar_cell_label,ar_coorx,ar_coory,ar_lambda,ar_Xc,ar_dam,ar_tan_beta,ar_tan_beta_channel,ar_L,ar_Ks,\
    ar_theta_r,ar_theta_s,ar_n_o,ar_n_c,\
    ar_cell_down,ar_pVs_t0,ar_Vo_t0,ar_Qc_t0,ar_kc\
        =pm.read_cell_parameters(file_cell_param)

    #~~~~Number of cell in the catchment
    nb_cell=len(ar_cell_label)
    #~~~~Computation of cell order
    ar_label_sort=pm.sort_cell(ar_cell_label,ar_cell_down)

    print 'Treatment of cell slopes'
    #Adjust the outlet slope if equal to zero
    if ar_tan_beta[ar_label_sort[-1]]==0:
        ar_tan_beta[ar_label_sort[-1]]=1./X

    for cell1 in ar_label_sort:
        print cell1
        cell=np.where(ar_cell_label==cell1)[0][0]

        li_label_path=[cell]
        li_slope_path=[ar_tan_beta[cell]]

        down_cell=ar_cell_down[cell]

        #~~~~ Extract the arrays of (i) all cells in the path (ar_label_path) (ii) all corresponding slopes (ar_slope_path)
        while down_cell>-1:
            li_label_path.append(down_cell)
            li_slope_path.append(ar_tan_beta[down_cell])
            down_cell=ar_cell_down[down_cell]
        ar_label_path=np.array(li_label_path)
        ar_slope_path=np.array(li_slope_path)

        #~~~~ Go into the slope vector to detect the series of zero slopes
        #The method first consists of creating an array containing the index of begining and end of zero series
        #Example: ar_slope_path=np.array([1, 2, 0, 0, 1, 0, 1, 0, 0, 0, 0])
        #--> ar_ind_start=array([2,5,7])
        #--> ar_ind_end=array([4,6,-1])
        li_ind_start=[]
        li_ind_end=[]
        #print len(ar_slope_path), np.where(ar_slope_path==0.)
        for i in np.arange(len(ar_slope_path)):
            if i==0:
                if ar_slope_path[i]==0.:
                    li_ind_start.append(i)
            elif i==np.arange(len(ar_slope_path))[-1]:
                if ar_slope_path[i]==0. and ar_slope_path[i-1]==0.:
                    li_ind_end.append(-1)
                if ar_slope_path[i]!=0. and ar_slope_path[i-1]==0.:
                    li_ind_end.append(i)
            else:
                if ar_slope_path[i]==0. and ar_slope_path[i-1]!=0:
                    li_ind_start.append(i)
                if ar_slope_path[i]!=0. and ar_slope_path[i-1]==0.:
                    li_ind_end.append(i)
        ar_ind_start=np.array(li_ind_start,int)
        ar_ind_end=np.array(li_ind_end,int)
        if len(ar_ind_start)-len(ar_ind_end)!=0:
            print 'problem length'
            print 'ar_ind_start:',len(ar_ind_start),'ar_ind_end:',len(ar_ind_end)
            print ar_ind_start
            print ar_ind_end
            stop
        a=ar_ind_end-ar_ind_start
        if len(np.where(a<0)[0])!=0:
            print 'problem index'
            stop

        #Then the slope are changed according to the defined index arrays
        if len(ar_ind_start)>0:
            print 'Number of sections with zero slopes:',len(ar_ind_start)
        for i in np.arange(len(ar_ind_start)):
            #Compute the length of the zero path
            if ar_ind_end[i]!=-1:
                length_path=ar_ind_end[i]-ar_ind_start[i]
            else:
                length_path=len(ar_slope_path)-ar_ind_start[i]
            #Compute the corresponding slope
            slope=1./(length_path*X)
            #Replace the values of slope in the initial ar_tan_beta
            #select the cell labels
            ar_label_cell_zero=ar_label_path[ar_ind_start[i]:ar_ind_end[i]]
            #select the cell index (if different from the label)
            ar_ind_cell_zero=np.array(ar_label_cell_zero,int)
            n=-1
            for label in ar_label_cell_zero:
                n=n+1
                ind=np.where(ar_cell_label==label)[0][0]
                ar_ind_cell_zero[n]=ind

            #Change the values
            ar_tan_beta[ar_ind_cell_zero]=slope

    print 'Treatment of channel slopes'
    #Adjust the outlet slope if equal to zero
    if ar_tan_beta_channel[ar_label_sort[-1]]==0:
        ar_tan_beta_channel[ar_label_sort[-1]]=1./X

    for cell1 in ar_label_sort:
        print cell1
        cell=np.where(ar_cell_label==cell1)[0][0]

        li_label_path=[cell]
        li_slope_path=[ar_tan_beta_channel[cell]]

        down_cell=ar_cell_down[cell]

        #~~~~ Extract the arrays of (i) all cells in the path (ar_label_path) (ii) all corresponding slopes (ar_slope_path)
        while down_cell>-1:
            li_label_path.append(down_cell)
            li_slope_path.append(ar_tan_beta_channel[down_cell])
            down_cell=ar_cell_down[down_cell]
        ar_label_path=np.array(li_label_path)
        ar_slope_path=np.array(li_slope_path)

        #~~~~ Go into the slope vector to detect the series of zero slopes
        #The method first consists of creating an array containing the index of begining and end of zero series
        #Example: ar_slope_path=np.array([1, 2, 0, 0, 1, 0, 1, 0, 0, 0, 0])
        #--> ar_ind_start=array([2,5,7])
        #--> ar_ind_end=array([4,6,-1])
        li_ind_start=[]
        li_ind_end=[]
        #print len(ar_slope_path), np.where(ar_slope_path==0.)
        for i in np.arange(len(ar_slope_path)):
            if i==0:
                if ar_slope_path[i]==0.:
                    li_ind_start.append(i)
            elif i==np.arange(len(ar_slope_path))[-1]:
                if ar_slope_path[i]==0. and ar_slope_path[i-1]==0.:
                    li_ind_end.append(-1)
                if ar_slope_path[i]!=0. and ar_slope_path[i-1]==0.:
                    li_ind_end.append(i)
            else:
                if ar_slope_path[i]==0. and ar_slope_path[i-1]!=0:
                    li_ind_start.append(i)
                if ar_slope_path[i]!=0. and ar_slope_path[i-1]==0.:
                    li_ind_end.append(i)
        ar_ind_start=np.array(li_ind_start,int)
        ar_ind_end=np.array(li_ind_end,int)
        if len(ar_ind_start)-len(ar_ind_end)!=0:
            print 'problem length'
            print 'ar_ind_start:',len(ar_ind_start),'ar_ind_end:',len(ar_ind_end)
            print ar_ind_start
            print ar_ind_end
            stop
        a=ar_ind_end-ar_ind_start
        if len(np.where(a<0)[0])!=0:
            print 'problem index'
            stop

        #Then the slope are changed according to the defined index arrays
        if len(ar_ind_start)>0:
            print 'Number of sections with zero slopes:',len(ar_ind_start)
        for i in np.arange(len(ar_ind_start)):
            #Compute the length of the zero path
            if ar_ind_end[i]!=-1:
                length_path=ar_ind_end[i]-ar_ind_start[i]
            else:
                length_path=len(ar_slope_path)-ar_ind_start[i]
            #Compute the corresponding slope
            slope=1./(length_path*X)
            #Replace the values of slope in the initial ar_tan_beta_channel
            #select the cell labels
            ar_label_cell_zero=ar_label_path[ar_ind_start[i]:ar_ind_end[i]]
            #select the cell index (if different from the label)
            ar_ind_cell_zero=np.array(ar_label_cell_zero,int)
            n=-1
            for label in ar_label_cell_zero:
                n=n+1
                ind=np.where(ar_cell_label==label)[0][0]
                ar_ind_cell_zero[n]=ind

            #Change the values
            ar_tan_beta_channel[ar_ind_cell_zero]=slope



    #~~~~~~Write parameter file~~~~~~#
    tab_param=np.zeros((len(ar_cell_label),nb_param))
    tab_param[:,0]=ar_cell_label
    tab_param[:,1]=ar_coorx
    tab_param[:,2]=ar_coory
    tab_param[:,3]=ar_lambda
    tab_param[:,4]=ar_Xc
    tab_param[:,5]=ar_dam
    tab_param[:,6]=ar_tan_beta
    tab_param[:,7]=ar_tan_beta_channel
    tab_param[:,8]=ar_L
    tab_param[:,9]=ar_Ks
    tab_param[:,10]=ar_theta_r
    tab_param[:,11]=ar_theta_s
    tab_param[:,12]=ar_n_o
    tab_param[:,13]=ar_n_c
    tab_param[:,14]=ar_cell_down
    tab_param[:,15]=ar_pVs_t0
    tab_param[:,16]=ar_Vo_t0
    tab_param[:,17]=ar_Qc_t0
    tab_param[:,18]=ar_kc

    np.savetxt(file_cell_param_out, tab_param)


def subcatch(ini_file='subcatch.ini'):
    """
    Create a subcatchment parameter file.

    Extract a subcatchment parameter file from an existing parameter
    file by defining the location of the sub-catchment outlet. The
    closest channel cell to the given coordinates is selected by
    default, based on a nearest neighbour search.

    A new parameter file for the subcatchment is written and a plot of
    the full catchment showing the location of the derived
    subcatchment is also produced.

    Parameters
    ----------
    ini_file : string
        The name of the ini file that specifies the input and output files,
        coordinates of the subcatchment outlet, the number of parameters and
        the grid dimension.

    Returns
    -------
    Nothing

    Notes
    -----
    This function creates a parameter file where the cell numbering doesn't
    follow the conventional North to South and West to East ordering. The
    cells are labelled in order from the outlet to the most distant upstream
    cell.

    """
    config.read(ini_file)
    print 'Read the file ', ini_file

    file_in = config.get('file_in', 'file_in')

    file_out = config.get('file_out', 'file_out')

    picture_out = config.get('picture_out', 'picture_out')

    Xoutlet = config.getfloat('coord_outlet', 'Xoutlet')
    Youtlet = config.getfloat('coord_outlet', 'Youtlet')

    nb_param = config.getfloat('flags', 'nb_param')
    X = config.getfloat('flags', 'X')

    #Reading of parameter file
    print 'Reading parameter file'
    ar_cell_label, ar_coorx, ar_coory, ar_lambda, ar_Xc, ar_dam, ar_tan_beta, \
    ar_tan_beta_channel, ar_L, ar_Ks, ar_theta_r, ar_theta_s, ar_n_o, ar_n_c, \
    ar_cell_down, ar_pVs_t0, ar_Vo_t0, ar_Qc_t0, ar_kc \
    = pm.read_cell_parameters(file_in)

    #Search for the cell close to the coordinates
    print 'Search for the outlet cell'
    cell_outlet = find_cell_coordinates(ar_cell_label, Xoutlet,
                                        Youtlet, ar_coorx, ar_coory, ar_lambda)

    #Search for the catchment cells
    print 'Search for the catchment cells'
    subcatch_label = all_up_cell(cell_outlet, ar_cell_down, ar_cell_label)

    #Select the subcatchmnent parameters
    print 'Select the subcatchmnent parameters'
    tab_param = np.zeros((len(subcatch_label),nb_param))
    new_label = np.arange(len(subcatch_label))

    tab_param[:,0] = new_label#ar_cell_label[subcatch_label]
    tab_param[:,1] = ar_coorx[subcatch_label]
    tab_param[:,2] = ar_coory[subcatch_label]
    tab_param[:,3] = ar_lambda[subcatch_label]
    tab_param[:,4] = ar_Xc[subcatch_label]
    tab_param[:,5] = ar_dam[subcatch_label]
    tab_param[:,6] = ar_tan_beta[subcatch_label]
    tab_param[:,7] = ar_tan_beta_channel[subcatch_label]
    tab_param[:,8] = ar_L[subcatch_label]
    tab_param[:,9] = ar_Ks[subcatch_label]
    tab_param[:,10] = ar_theta_r[subcatch_label]
    tab_param[:,11] = ar_theta_s[subcatch_label]
    tab_param[:,12] = ar_n_o[subcatch_label]
    tab_param[:,13] = ar_n_c[subcatch_label]
    for i in range(len(subcatch_label)):
        if i == 0:
            tab_param[i,14] = -9999.0
        else:
            ind = np.where(ar_cell_label[subcatch_label]
                           == ar_cell_down[subcatch_label][i])

            tab_param[i,14] = new_label[ind]

    tab_param[:,15]=ar_pVs_t0[subcatch_label]
    tab_param[:,16]=ar_Vo_t0[subcatch_label]
    tab_param[:,17]=ar_Qc_t0[subcatch_label]
    tab_param[:,18]=ar_kc[subcatch_label]

    #~~~~~~Write parameter file~~~~~~#
    np.savetxt(file_out, tab_param)

    ar_image=ar_cell_label*0.
    ar_image[subcatch_label]=1.
    ar_image[ar_lambda==1.]=10.
    ar_image[cell_outlet]=5.
    field_map(ar_image, ar_coorx, ar_coory, X, picture_out, 'Subcatchment')

def new_param(ini_file='new_param.ini'):
    """
    * Objective
        Modifies the param file by multiplying some variables (L, Ks, n_o,n_c)
        or replacing by a given new value for the initial level of reservoir in percent
    * Input
      - file_in: parameter file to be modified
      - the other parameters are assigned by default and have to be changed according to user's choice.
    * Output
      - file_out: new parameter file
    """
    ### READ THE INI FILE ###
    config.read(ini_file)
    print 'Read the file ',ini_file
    ##~~~~~~ file_in ~~~~~~##
    file_in=config.get('file_in','file_in')
    ##~~~~~~ file_out ~~~~~~##
    file_out=config.get('file_out','file_out')

    ##~~~~~~ factor_values ~~~~~##
    fac_L=config.getfloat('factor_values','fac_L')
    fac_Ks=config.getfloat('factor_values','fac_KS')
    fac_n_o=config.getfloat('factor_values','fac_n_o')
    fac_n_c=config.getfloat('factor_values','fac_n_c')

    ##~~~~~~ new_initial_values ~~~~~##
    new_pVs_t0=config.getfloat('new_initial_values','new_pVs_t0')
    new_Vo_t0=config.getfloat('new_initial_values','new_Vo_t0')
    new_Qc_t0=config.getfloat('new_initial_values','new_Qc_t0')

    ##~~~~~~ flags ~~~~~~##
    nb_param=config.getfloat('flags','nb_param')

    #Reading of parameter file
    print 'Reading parameter file'
    ar_cell_label,ar_coorx,ar_coory,ar_lambda,ar_Xc,ar_dam,ar_tan_beta,ar_tan_beta_channel,ar_L,ar_Ks,\
    ar_theta_r,ar_theta_s,ar_n_o,ar_n_c,\
    ar_cell_down,ar_pVs_t0,ar_Vo_t0,ar_Qc_t0,ar_kc\
    =pm.read_cell_parameters(file_in)

    #~~~~~~Change in parameters~~~~~~#
    #Multiplying factors for L, Ks, n_o and n_c
    if fac_L!=1.:
        print 'Change L'
        ar_L=ar_L*fac_L
    if fac_Ks!=1.:
        print 'Change Ks'
        ar_Ks=ar_Ks*fac_Ks
    if fac_n_o!=1.:
        print 'Change n_o'
        ar_n_o=ar_n_o*fac_n_o
    if fac_n_c!=1.:
        print 'Change n_c'
        ar_n_c=ar_n_c*fac_n_c
    #New values for pVs_t0, Vo_t0 and Qc_t0
    if new_pVs_t0!=ar_pVs_t0[0]:
        print 'Change pVs_t0'
        ar_pVs_t0=ar_pVs_t0*0.+new_pVs_t0
    if new_Vo_t0!=ar_Vo_t0[0]:
        print 'Change pVs_t0'
        ar_Vo_t0=ar_Vo_t0*0.+new_Vo_t0
    if new_Qc_t0!=ar_Qc_t0[0]:
        print 'Change pVc_t0'
        ar_Qc_t0=ar_Qc_t0*0.+new_Qc_t0

    #~~~~~~Write parameter file~~~~~~#
    tab_param=np.zeros((len(ar_cell_label),nb_param))
    tab_param[:,0]=ar_cell_label
    tab_param[:,1]=ar_coorx
    tab_param[:,2]=ar_coory
    tab_param[:,3]=ar_lambda
    tab_param[:,4]=ar_Xc
    tab_param[:,5]=ar_dam
    tab_param[:,6]=ar_tan_beta
    tab_param[:,7]=ar_tan_beta_channel
    tab_param[:,8]=ar_L
    tab_param[:,9]=ar_Ks
    tab_param[:,10]=ar_theta_r
    tab_param[:,11]=ar_theta_s
    tab_param[:,12]=ar_n_o
    tab_param[:,13]=ar_n_c
    tab_param[:,14]=ar_cell_down
    tab_param[:,15]=ar_pVs_t0
    tab_param[:,16]=ar_Vo_t0
    tab_param[:,17]=ar_Qc_t0
    tab_param[:,18]=ar_kc

    np.savetxt(file_out, tab_param)


def connect_external_flow(ini_file='connect_external_flow.ini'):
    """
    * Objective
      to link the location of external flows to the river network by transforming the initially non-channel cells into channel cells
      ar_lambda is modified (value from 0 to 1) ar_n_c is modified the value of the new channel cells are taken equal to the nearest existing channel cell
    """
    ### READ THE INI FILE ###
    config.read(ini_file)
    print 'Read the file ',ini_file
    ##~~~~~~ file_in ~~~~~~##
    file_in=config.get('file_in','file_in')
    ##~~~~~~ file_out ~~~~~~##
    file_out=config.get('file_out','file_out')

    ##~~~~~~ external_flow ~~~~~~##
    Xext_flow=config.getfloat('external_flow','Xext_flow')
    Yext_flow=config.getfloat('external_flow','Yext_flow')

    ##~~~~~~ flags ~~~~~~##
    nb_param=config.getfloat('flags','nb_param')

    #Reading of parameter file
    print 'Reading parameter file'
    ar_cell_label,ar_coorx,ar_coory,ar_lambda,ar_Xc,ar_dam,ar_tan_beta,ar_tan_beta_channel,ar_L,ar_Ks,\
    ar_theta_r,ar_theta_s,ar_n_o,ar_n_c,\
    ar_cell_down,ar_pVs_t0,ar_Vo_t0,ar_Qc_t0,ar_kc\
    =pm.read_cell_parameters(file_in)

    print 'Connect external flows to the network'
    ar_lambda,ar_n_c=link_channel_cell(ar_cell_label,ar_coorx,ar_coory,ar_lambda,ar_cell_down,ar_n_c,Xext_flow,Yext_flow)

    #~~~~~~Write parameter file~~~~~~#
    tab_param=np.zeros((len(ar_cell_label),nb_param))
    tab_param[:,0]=ar_cell_label
    tab_param[:,1]=ar_coorx
    tab_param[:,2]=ar_coory
    tab_param[:,3]=ar_lambda
    tab_param[:,4]=ar_Xc
    tab_param[:,5]=ar_dam
    tab_param[:,6]=ar_tan_beta
    tab_param[:,7]=ar_tan_beta_channel
    tab_param[:,8]=ar_L
    tab_param[:,9]=ar_Ks
    tab_param[:,10]=ar_theta_r
    tab_param[:,11]=ar_theta_s
    tab_param[:,12]=ar_n_o
    tab_param[:,13]=ar_n_c
    tab_param[:,14]=ar_cell_down
    tab_param[:,15]=ar_pVs_t0
    tab_param[:,16]=ar_Vo_t0
    tab_param[:,17]=ar_Qc_t0
    tab_param[:,18]=ar_kc

    np.savetxt(file_out, tab_param)


def initial_pVs_Vo_Qc_from_simu(ini_file='initial_pVs_Vo_Qc_from_simu.ini'):
    """
    * Objective

     """
    ### READ THE INI FILE ###
    config.read(ini_file)
    print 'Read the file ',ini_file

    ##~~~~~~ file_in ~~~~~~##
    file_in=config.get('file_in','file_in')
    file_in_global=config.get('file_in','file_in_global')
    file_h5=config.get('file_in','file_h5')

    ##~~~~~~ file_out ~~~~~~##
    file_out=config.get('file_out','file_out')

    ##~~~~~~ variables ~~~~~~##
    time_step=config.getint('variables','time_step')
    fac_L_simu=config.getfloat('variables','fac_L_simu')
    fac_Ks_simu=config.getfloat('variables','fac_Ks_simu')
    fac_n_o_simu=config.getfloat('variables','fac_n_o_simu')
    fac_n_c_simu=config.getfloat('variables','fac_n_c_simu')

    ##~~~~~~ flags ~~~~~~##
    nb_param=config.getfloat('flags','nb_param')

    #--Read and compute the parameters to have the values of parameter and ar_Vsm
    #~~~~Read Global parameters file
    #~~~~Read Global parameters file
    X,Dt,alpha_s,alpha_o,alpha_c,A_thres,W_min,W_max\
      =pm.read_global_parameters(file_in_global)
    #~~~~Read Cell parameters file
    ar_cell_label,ar_coorx,ar_coory,ar_lambda,ar_Xc,ar_dam,ar_tan_beta,ar_tan_beta_channel,ar_L,ar_Ks,\
    ar_theta_r,ar_theta_s,ar_n_o,ar_n_c,\
    ar_cell_down,ar_pVs_t0,ar_Vo_t0,ar_Qc_t0,ar_kc\
        =pm.read_cell_parameters(file_in)
    #~~~~Number of cell in the catchment
    nb_cell=len(ar_cell_label)
    #~~~~Computation of cell order
    ar_label_sort=pm.sort_cell(ar_cell_label,ar_cell_down)
    #~~~~Computation of upcells
    li_cell_up=pm.direct_up_cell(ar_cell_label,ar_cell_down,ar_label_sort)
    #~~~~Computation of drained area
    ar_A_drained=pm.drained_area(ar_label_sort,li_cell_up,X)
    #~~~~Modifies the values of the parameters
    ar_L1=ar_L*fac_L_simu
    ar_Ks1=ar_Ks*fac_Ks_simu
    ar_n_o1=ar_n_o*fac_n_o_simu
    ar_n_c1=ar_n_c*fac_n_c_simu
    #~~~~Computation of model parameters from physical parameters
    ar_Vsm, ar_b_s, ar_b_o, ar_W, ar_b_c\
      =pm.compute_cell_param(X,ar_Xc,Dt,alpha_s,alpha_o,alpha_c,nb_cell,\
                              A_thres,W_max,W_min,\
                              ar_lambda,ar_tan_beta,ar_tan_beta_channel,ar_L1,\
                              ar_Ks1,ar_theta_r,ar_theta_s,ar_n_o1,ar_n_c1,\
                              ar_A_drained)

    #Read the soil volume file
    ndar_Vs=np.array(ut.read_one_array_hdf(file_h5,'/Soil/','V_s'))
    #Read the overland volume file
    ndar_Vo=np.array(ut.read_one_array_hdf(file_h5,'/Overland/','V_o'))
    #Read the channel dischargefile
    ndar_Qc=np.array(ut.read_one_array_hdf(file_h5,'/Channel/','Qc_out'))

    ar_Vs=ndar_Vs[time_step,:]
    ar_pVs_t0=ar_Vs/ar_Vsm*100.
    ar_Vo_t0=ndar_Vo[time_step,:]
    ar_Qc_t0=ndar_Qc[time_step,:]


    #~~~~~~Write parameter file~~~~~~#
    tab_param=np.zeros((len(ar_cell_label),nb_param))

    tab_param[:,0]=ar_cell_label
    tab_param[:,1]=ar_coorx
    tab_param[:,2]=ar_coory
    tab_param[:,3]=ar_lambda
    tab_param[:,4]=ar_Xc
    tab_param[:,5]=ar_dam
    tab_param[:,6]=ar_tan_beta
    tab_param[:,7]=ar_tan_beta_channel
    tab_param[:,8]=ar_L
    tab_param[:,9]=ar_Ks
    tab_param[:,10]=ar_theta_r
    tab_param[:,11]=ar_theta_s
    tab_param[:,12]=ar_n_o
    tab_param[:,13]=ar_n_c
    tab_param[:,14]=ar_cell_down
    tab_param[:,15]=ar_pVs_t0
    tab_param[:,16]=ar_Vo_t0
    tab_param[:,17]=ar_Qc_t0
    tab_param[:,18]=ar_kc

    np.savetxt(file_out, tab_param)

def mean_simuVsi(ini_file='mean_simuVsi.ini'):
    """
    * Objective

    """

    ### READ THE INI FILE ###
    config.read(ini_file)
    print 'Read the file ',ini_file

    ##~~~~~~ file_in ~~~~~~##
    file_in=config.get('file_in','file_in')
    file_in_global=config.get('file_in','file_in_global')
    file_h5=config.get('file_in','file_h5')

    ##~~~~~~ file_out ~~~~~~##
    file_out=config.get('file_out','file_out')

    ##~~~~~~ variables ~~~~~~##
    mean_pVs_t0=config.getfloat('variables','mean_pVs_t0')
    fac_L_simu=config.getfloat('variables','fac_L_simu')
    fac_Ks_simu=config.getfloat('variables','fac_Ks_simu')
    fac_n_o_simu=config.getfloat('variables','fac_n_o_simu')
    fac_n_c_simu=config.getfloat('variables','fac_n_c_simu')

    ##~~~~~~ flags ~~~~~~##
    nb_param=config.getfloat('flags','nb_param')

    #--Read and compute the parameters to have the values of parameter and ar_Vsm
    #~~~~Read Global parameters file
    print 'Pretreatment of input data'
    #~~~~Read Global parameters file
    X,Dt,alpha_s,alpha_o,alpha_c,A_thres,W_min,W_max\
      =pm.read_global_parameters(file_in_global)
    #~~~~Read Cell parameters file
    ar_cell_label,ar_coorx,ar_coory,ar_lambda,ar_Xc,ar_dam,ar_tan_beta,ar_tan_beta_channel,ar_L,ar_Ks,\
    ar_theta_r,ar_theta_s,ar_n_o,ar_n_c,\
    ar_cell_down,ar_pVs_t0,ar_Vo_t0,ar_Qc_t0,ar_kc\
        =pm.read_cell_parameters(file_in)
    #~~~~Number of cell in the catchment
    nb_cell=len(ar_cell_label)
    #~~~~Computation of cell order
    ar_label_sort=pm.sort_cell(ar_cell_label,ar_cell_down)
    #~~~~Computation of upcells
    li_cell_up=pm.direct_up_cell(ar_cell_label,ar_cell_down,ar_label_sort)
    #~~~~Computation of drained area
    ar_A_drained=pm.drained_area(ar_label_sort,li_cell_up,X)
    #~~~~Modifies the values of the parameters
    ar_L1=ar_L*fac_L_simu
    ar_Ks1=ar_Ks*fac_Ks_simu
    ar_n_o1=ar_n_o*fac_n_o_simu
    ar_n_c1=ar_n_c*fac_n_c_simu
    #~~~~Computation of model parameters from physical parameters
    ar_Vsm, ar_b_s, ar_b_o, ar_W, ar_b_c\
      =pm.compute_cell_param(X,ar_Xc,Dt,alpha_s,alpha_o,alpha_c,nb_cell,\
                              A_thres,W_max,W_min,\
                              ar_lambda,ar_tan_beta,ar_tan_beta_channel,ar_L1,\
                              ar_Ks1,ar_theta_r,ar_theta_s,ar_n_o1,ar_n_c1,\
                              ar_A_drained)

    #Read the soil volume file
    ndar_Vs=np.array(ut.read_one_array_hdf(file_h5,'/Soil/','V_s'))

    #Read the file of catchment saturation rates
    ndar_Vs_sat=ndar_Vs/ar_Vsm*100.
    tab_rate=np.average(ndar_Vs_sat,axis=1)
    tab_rate_sort=np.sort(tab_rate)
    indice_sort=np.argsort(tab_rate)

    #Look for the rate closest to the expected mean value
    if mean_pVs_t0<tab_rate_sort[0]:
        ind=0
        print 'mean_pVs_t0 expected:', mean_pVs_t0,'effective:',tab_rate_sort[0]
    elif mean_pVs_t0>tab_rate_sort[-1]:
        ind=-1
        print 'mean_pVs_t0 expected:', mean_pVs_t0,' effective',tab_rate_sort[-1]
    else:
        loop=True
        i=-1
        while loop:
            i=i+1
            if mean_pVs_t0>=tab_rate_sort[i] and mean_pVs_t0<tab_rate_sort[i+1]:
                ind=i
                loop=False
                print 'mean_pVs_t0 expected:', mean_pVs_t0,' effective:',tab_rate_sort[i]

    ind_end=indice_sort[ind]
    print ind,ind_end
    ar_Vs=ndar_Vs[ind_end,:]
    ar_Vsi=ar_Vs/ar_Vsm*100.
    print ar_Vsi
    ar_pVs_t0=ar_Vsi

    #~~~~~~Write parameter file~~~~~~#
    tab_param=np.zeros((len(ar_cell_label),nb_param))

    tab_param[:,0]=ar_cell_label
    tab_param[:,1]=ar_coorx
    tab_param[:,2]=ar_coory
    tab_param[:,3]=ar_lambda
    tab_param[:,4]=ar_Xc
    tab_param[:,5]=ar_dam
    tab_param[:,6]=ar_tan_beta
    tab_param[:,7]=ar_tan_beta_channel
    tab_param[:,8]=ar_L
    tab_param[:,9]=ar_Ks
    tab_param[:,10]=ar_theta_r
    tab_param[:,11]=ar_theta_s
    tab_param[:,12]=ar_n_o
    tab_param[:,13]=ar_n_c
    tab_param[:,14]=ar_cell_down
    tab_param[:,15]=ar_pVs_t0
    tab_param[:,16]=ar_Vo_t0
    tab_param[:,17]=ar_Qc_t0
    tab_param[:,18]=ar_kc

    np.savetxt(file_out, tab_param)

########################################################################################
######## Surbroutines
########################################################################################

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

def distance(x1, y1, x2, y2):
    """
    Compute the distance between two points

    """
    dist = ((x1-x2)**2 + (y1-y2)**2)**0.5
    return dist

def find_cell_coordinates(ar_cell_label, Xtarget, Ytarget,
                          ar_coorx, ar_coory, ar_lambda, channel=True):
    """
    Find the closest cell to (Xtarget, Ytarget)

    """
    tab_x = np.unique(ar_coorx)
    X = abs(tab_x[0] - tab_x[1])
    dist_max = 3*X
    dist_min = dist_max
    nb_cell = len(ar_cell_label)
    cell_outlet = -999.9

    for i in range(nb_cell):
        dist = distance(Xtarget, Ytarget, ar_coorx[i], ar_coory[i])
        if channel:
            if dist < dist_min and ar_lambda[i] == 1.:
                dist_min = dist
                cell_outlet = ar_cell_label[i]
        else:
            if dist < dist_min:
                dist_min = dist
                cell_outlet = ar_cell_label[i]


    if cell_outlet < 0:
        err = 'Invalid outlet coordinates specified: X=%f, Y=%f' \
              % (Xtarget, Ytarget)
        raise ValueError, err

    return cell_outlet

def link_channel_cell(ar_cell_label,ar_coorx,ar_coory,ar_lambda,ar_cell_down,ar_n_c,Xext_flow,Yext_flow):
    """
    Change the ar_lambda to connnect the tunnel to the channel. Zero values are changed in 1 values.
    The values of the manning coefficient for the changed cells is taken as the manning value of the closest channel cell.
    """
    cell=find_cell_coordinates(ar_cell_label,Xext_flow,Yext_flow,ar_coorx,ar_coory,ar_lambda,channel=False)
    hillslope=True
    li_ind=[]
    cc=0.
    while hillslope:
        ind=np.where(ar_cell_label==cell)
        if ar_lambda[ind]==1.:
            hillslope=False
            last_ind=ind
        else:
            cc=cc+1
            print 'Cell',cell,'has been conected to the channel network via cell',ar_cell_down[ind]
            li_ind.append(ind)
            ar_lambda[ind]=1.
            cell=ar_cell_down[ind]
    for i in li_ind:
        ar_n_c[i]=ar_n_c[last_ind]
    if cc==0.:
        print 'External flows already connected'
    return ar_lambda,ar_n_c

#####################################
## SUBROUTINES                      #
#####################################

##~~ Compute all the upper cells ~~##
def direct_up_cell(ar_cell, ar_cell_down, ar_cell_label):
    """
    Calculate cells immediately uphill of target cells.

    Parameters
    ----------
    ar_cell : array_like
        Array of cell labels specifying the target cells to search from.
    ar_cell_down : array_like
        An array of cell labels, identifying the immediate downstream
        neighbour of each cell in the catchment.
    ar_cell_label : array_like
        The cell labels of the catchment. Ordered to match the information
        in `ar_cell_down`.

    Returns
    -------
    ar_out : ndarray
        An ordered array of immediate upstream cell labels.

    """
    ar_out = np.array([],int)

    for i in range(len(ar_cell)):
        ind = np.where(ar_cell_down == ar_cell[i])

        if np.size(ar_out) == 0:
            ar_out = ar_cell_label[ind]
        else:
            if np.size(ar_cell_label[ind]) > 0:
                ar_out = np.concatenate((ar_out, ar_cell_label[ind]))

    return ar_out

def all_up_cell(cell, ar_cell_down, ar_cell_label):
    """
    Calculate the network of uphill cells.

    Determines the labels of all cells that drain into the target
    cell. The cell connectivity is used to trace the network of cells.

    Parameters
    ----------
    cell : int
        The integer label identifying the target cell.
    ar_cell_down : array_like
        An array of cell labels, identifying the immediate downstream
        neighbour of each cell in the catchment.
    ar_cell_label : array_like
        The cell labels of the catchment. Ordered to match the information
        in `ar_cell_down`.

    Returns
    -------
    upcells : ndarray
        An ordered array of upstream cell labels. The order is in increasing
        distance from the target cell.

    """
    a = np.ones(1,int)*cell
    upcells = a

    while np.size(a)!= 0:
        a = direct_up_cell(a, ar_cell_down, ar_cell_label)
        upcells = np.concatenate((upcells, a))

    return upcells

##~~ Plot a field map ~~##
def field_map(ar_field, ar_coorx, ar_coory, X, picture_out, title, flip=0):
    """
    Plot an overview depicting the location of the subcatchment.

    """

    import matplotlib.numerix.ma as M

    max_val=max(ar_field)

    xmin=min(ar_coorx);xmax=max(ar_coorx)
    ymin=min(ar_coory);ymax=max(ar_coory)
    step=X
    nx=(xmax-xmin)/step+1
    ny=(ymax-ymin)/step+1

    ar_indx=np.array((ar_coorx-xmin)/step,int)
    ar_indy=np.array((ar_coory-ymin)/step,int)

    ar_map=sp.ones((ny,nx))*-99.9
    ar_map[ar_indy,ar_indx]=ar_field

    if flip==1:
        ar_map=np.flipud(ar_map)

    ar_map2 = M.masked_where(ar_map <0, ar_map)


    ut.check_file_exist(picture_out)

    pl.clf()
    pl.imshow(ar_map2, interpolation='Nearest',
              origin='lower', vmax=max_val,vmin=0)

    pl.title(title)
    pl.colorbar()
    pl.savefig(picture_out)

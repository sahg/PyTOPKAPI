""" model.py
Main programm of the TOPKAPI model.
"""

__author__ = "Theo Vischel"
__version__ = "$Revision: 1.0 $"
__date__ = "$Date: 08/08/2007 $"


#Python modules
import numpy as np

from ConfigParser import SafeConfigParser
config = SafeConfigParser()

#External modules from pytopkapi
#pretreatment: used for subroutines to read the column type files.
from pytopkapi import pretreatment as pm


def run(ini_file='zero_slope_management.ini'):

    ##================================##
    ##  Read the input file (*.ini)   ##
    ##================================##
    config.read(ini_file)
    print 'Read the file ',ini_file
    ##~~~~~~~~~~~ input files ~~~~~~~~~~~##
    #Param
    file_global_param=config.get('input_files','file_global_param')
    file_cell_param=config.get('input_files','file_cell_param')
    file_cell_param_out=config.get('output_files','file_cell_param_out')
    ##~~~~~~ numerical_values ~~~~~~##
    nb_param=config.getfloat('numerical_values','nb_param')

    ##============================##
    ## Pretreatment of input data ##
    ##============================##
    print 'Pretreatment of input data'
    #~~~~Read Global parameters file
    X,Dt,alpha_s,alpha_o,alpha_c,A_thres,W_min,W_max\
      =pm.read_global_parameters(file_global_param)
    #~~~~Read Cell parameters file
    ar_cell_label,ar_coorx,ar_coory,ar_lambda,ar_Xc,ar_dam,ar_tan_beta,ar_tan_beta_channel,ar_L,ar_Ks,\
    ar_theta_r,ar_theta_s,ar_n_o,ar_n_c,\
    ar_cell_down,ar_pVs_t0,ar_Vo_t0,ar_Qc_t0,ar_kc\
        =pm.read_cell_parameters(file_cell_param)

    #~~~~Number of cell in the catchment
    nb_cell=len(ar_cell_label)
    #~~~~Computation of cell order
    ar_label_sort=pm.sort_cell(ar_cell_label,ar_cell_down)
    #~~~~Computation of upcells
    li_cell_up=pm.direct_up_cell(ar_cell_label,ar_cell_down,ar_label_sort)
    #~~~~Computation of drained area
    ar_A_drained=pm.drained_area(ar_label_sort,li_cell_up,X)
    #~~~~Computation of model parameters from physical parameters
    ar_Vsm, ar_b_s, ar_b_o, ar_W, ar_b_c\
      =pm.compute_cell_param(X,ar_Xc,Dt,alpha_s,alpha_o,alpha_c,nb_cell,\
                              A_thres,W_max,W_min,\
                              ar_lambda,ar_tan_beta,ar_tan_beta_channel,ar_L,\
                              ar_Ks,ar_theta_r,ar_theta_s,ar_n_o,ar_n_c,\
                              ar_A_drained)

    print 'Treatment of cell slopes'
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

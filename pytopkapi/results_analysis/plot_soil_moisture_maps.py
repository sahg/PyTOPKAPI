""" soil_moisture_maps.py
Draw the maps of soil moisture from the outputs of the TOPKAPI simulations
"""
import numpy as np
import pylab as pl
import numpy.ma as M

import pytopkapi.pretreatment as pm
import pytopkapi.utils as ut

from ConfigParser import SafeConfigParser
config = SafeConfigParser()


def run(ini_file='plot_soil_moisture_maps.ini'):

    config.read(ini_file)
    print 'Read the file ',ini_file

    file_global_param=config.get('files','file_global_param')
    file_cell_param=config.get('files','file_cell_param')
    file_sim=config.get('files','file_sim')

    path_out=config.get('paths','path_out')

    fac_L=config.getfloat('calib_params','fac_L')
    fac_Ks=config.getfloat('calib_params','fac_Ks')
    fac_n_o=config.getfloat('calib_params','fac_n_o')
    fac_n_c=config.getfloat('calib_params','fac_n_c')

    t1=config.getfloat('flags','t1')
    t2=config.getfloat('flags','t2')
    variable=config.getfloat('flags','variable')

    ##~~~~~~PROGRAM~~~~~##
    # Create the folder if it doesn't exist
    ut.check_folder_exist(path_out)

    #~~~~Read Global parameters file
    X,Dt,alpha_s,alpha_o,alpha_c,A_thres,W_min,W_max\
      =pm.read_global_parameters(file_global_param)
    #~~~~Read Cell parameters file
    ar_cell_label, ar_coorx, \
    ar_coory, ar_lambda, \
    ar_Xc, ar_dam, \
    ar_tan_beta, ar_tan_beta_channel, \
    ar_L0, ar_Ks0, \
    ar_theta_r, ar_theta_s, \
    ar_n_o0, ar_n_c0, \
    ar_cell_down, ar_pVs_t0, \
    ar_Vo_t0, ar_Qc_t0, \
    ar_kc, psi_b, lamda = pm.read_cell_parameters(file_cell_param)
    #~~~~Number of cell in the catchment
    nb_cell=len(ar_cell_label)
    #~~~~Computation of cell order
    ar_label_sort=pm.sort_cell(ar_cell_label,ar_cell_down)
    #~~~~Computation of upcells
    li_cell_up=pm.direct_up_cell(ar_cell_label,ar_cell_down,ar_label_sort)
    #~~~~Computation of drained area
    ar_A_drained=pm.drained_area(ar_label_sort,li_cell_up,X)
    #~~~~Modifies the values of the parameters
    ar_L=ar_L0*fac_L
    ar_Ks=ar_Ks0*fac_Ks
    ar_n_o=ar_n_o0*fac_n_o
    ar_n_c=ar_n_c0*fac_n_c
    #~~~~Computation of model parameters from physical parameters
    ar_Vsm, ar_b_s, ar_b_o, ar_W, ar_b_c\
      =pm.compute_cell_param(X,ar_Xc,Dt,alpha_s,alpha_o,alpha_c,nb_cell,\
                              A_thres,W_max,W_min,\
                              ar_lambda,ar_tan_beta,ar_tan_beta_channel,ar_L,\
                              ar_Ks,ar_theta_r,ar_theta_s,ar_n_o,ar_n_c,\
                              ar_A_drained)

    #Read of data from the outputs of TOPKAPI in hdf5 format
    ndar_Vs=np.array(ut.read_one_array_hdf(file_sim,'/Soil/','V_s'))

    #Assign the variables
    if variable==1:
        im_out=path_out+'field_Vs_'
        tab=ndar_Vs
    elif variable==2:
        im_out=path_out+'field_Vo_'
        tab=ndar_Vo
    elif variable==3:
        im_out=path_out+'field_Vo_bin_'
        tab=ndar_Vo
        tab[tab>0]=1.
    elif variable==4:
        im_out=path_out+'field_SWI_'
        tab=ndar_Vs/ar_Vsm*100.

    # Plot the maps
    for t in np.arange(int(t1),int(t2+1)):
        print 'Map time-step ', t
        image_out=im_out+ut.string(t,len(str(t2)))+'.png'
        field_map_ndar(tab,t,ar_coorx,ar_coory,X,image_out,variable)

def field_map_ndar(ndar_field,t,ar_coorx,ar_coory,X,image_out,variable):

    ar_field=ndar_field[t,:]
    max_val=int(np.max(ndar_field))
    if variable==4:
        max_val=100.
    xmin=min(ar_coorx);xmax=max(ar_coorx)
    ymin=min(ar_coory);ymax=max(ar_coory)
    step=X
    nx=(xmax-xmin)/step+1
    ny=(ymax-ymin)/step+1

    ar_indx=np.array((ar_coorx-xmin)/step,int)
    ar_indy=np.array((ar_coory-ymin)/step,int)

    ar_map=np.ones((ny,nx))*-99.9
    ar_map[ar_indy,ar_indx]=ar_field

    ar_map2 = M.masked_where(ar_map <0, ar_map)
    ut.check_file_exist(image_out)

    pl.clf()
    pl.axes(axisbg='gray')
    pl.imshow(ar_map2, cmap=pl.cm.RdBu,
              interpolation='Nearest', origin='lower', vmax=max_val, vmin=0)
    pl.title('time step= '+ut.string(t,len(str(t))))
    pl.colorbar()
    pl.savefig(image_out)

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
############################################
############################################
sys.path.append('C:/Theo/topkapi_model/programme/TOPKAPI_Aug07/svn_working_copy/lib/')
sys.path.append('C:/Theo/topkapi_model/programme/TOPKAPI_Aug07/svn_working_copy/lib/results_analysis')

#Internal modules
import numpy as np
import pylab as pl
from numpy import ma
import scipy.io as io
import tables as h5
import os
import os.path

#External modules
#Utilities
import utils as ut
#pretreatment_module: used for subroutines to read the parameter files
import pretreatment as pm
import result_analysis as ra
####################################
### MODIFYING THE PARAMETER FILE ###
####################################

def from_param_to_subcatch_param(file_in, file_out,Xoutlet,Youtlet,nom_image,X=1000):
    """
    * Objective:
    Extract a sub catchment parameter file from a catchment parameter file
    by defining the position of the outlet.
    The final outlet cell is the closest channel cell from the given coordinates

    * Input
      - file_in: parameter file to be modified
      
    * Output
      - file_out: parameter file for the selected subcatchment
      - nom_image: picture of the selected subcatchment
    """

    nb_param=17.
    
    #Reading of parameter file
    print 'Reading parameter file'
    ar_cell_label,ar_coorx,ar_coory,ar_lambda,ar_dam,ar_tan_beta,ar_L,ar_Ks,\
    ar_theta_r,ar_theta_s,ar_n_o,ar_n_c,\
    ar_cell_down,ar_pVs_t0,ar_Vo_t0,ar_Qc_t0,ar_kc\
    =pm.read_cell_parameters(file_in)

    #Search for the cell close to the coordinates
    print 'Search for the outlet cell'
    cell_outlet=find_cell_coordinates(ar_cell_label,Xoutlet,Youtlet,ar_coorx,ar_coory,ar_lambda)
    print cell_outlet
    
    #Search for the catchment cells
    print 'Search for the catchment cells'
    subcatch_label=ra.all_up_cell(cell_outlet,ar_cell_down,ar_cell_label)

    #Select the subcatchmnent parameters
    print 'Select the subcatchmnent parameters'
    tab_param=np.zeros((len(subcatch_label),nb_param))
    new_label=np.arange(len(subcatch_label))
    tab_param[:,0]=new_label#ar_cell_label[subcatch_label]
    tab_param[:,1]=ar_coorx[subcatch_label]
    tab_param[:,2]=ar_coory[subcatch_label]
    tab_param[:,3]=ar_lambda[subcatch_label]
    tab_param[:,4]=ar_dam[subcatch_label]
    tab_param[:,5]=ar_tan_beta[subcatch_label]
    tab_param[:,6]=ar_L[subcatch_label]
    tab_param[:,7]=ar_Ks[subcatch_label]
    tab_param[:,8]=ar_theta_r[subcatch_label]
    tab_param[:,9]=ar_theta_s[subcatch_label]
    tab_param[:,10]=ar_n_o[subcatch_label]
    tab_param[:,11]=ar_n_c[subcatch_label]
    for i in range(len(subcatch_label)):
        if i==0.:
            tab_param[i,12]=-9999.
        else:
            tab_param[i,12]=new_label[np.where(ar_cell_label[subcatch_label]==ar_cell_down[subcatch_label][i])]
    tab_param[:,13]=ar_pVs_t0[subcatch_label]
    tab_param[:,14]=ar_Vo_t0[subcatch_label]
    tab_param[:,15]=ar_Qc_t0[subcatch_label]
    tab_param[:,16]=ar_kc[subcatch_label]

    #~~~~~~Write parameter file~~~~~~#
    #'help io.write_array' for more info
    f = file(file_out, 'w')
    io.write_array(f, tab_param)
    f.close()

    image_out=os.path.split(file_out)[0]+'/'+nom_image+'.png'
    ar_image=ar_cell_label*0.;ar_image[subcatch_label]=1.;ar_image[ar_lambda==1.]=10.;ar_image[cell_outlet]=5.
    ra.field_map(ar_image,ar_coorx,ar_coory,X,image_out,'Subcatchment')

def from_param_to_new_param(file_in,file_out,fac_L=1.,fac_Ks=1.,fac_n_o=1.,fac_n_c=1.,new_pVs_t0=50.,new_Vo_t0=0.,new_pVc_t0=0.,channel_lesotho=False):
    """
    * Objective
        Modifies the param file by multiplying some variables (L, Ks, n_o,n_c)
        or replacing by a given new value for the initial level of reservoir in percent
    * Input
      - file_in: parameter file to be modified
      - the other parameters are assigned by default and have to be changed according to user's choice.
    * Output
      - file_out: new parameter file
    * Comment
      The boolean channel_lesotho calls a subroutine that was used in the case of the Liebenbergsvlei catchment
      to link the tunnel output to the river by transforming the initially non-channel cells into channel cells
    """

    nb_param=17.
    
    #Reading of parameter file
    print 'Reading parameter file'
    ar_cell_label,ar_coorx,ar_coory,ar_lambda,ar_dam,ar_tan_beta,ar_L,ar_Ks,\
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
        print ar_pVs_t0[0]
    if new_Vo_t0!=ar_Vo_t0[0]:
        print 'Change pVs_t0'
        ar_Vo_t0=ar_Vo_t0*0.+new_Vo_t0
        print ar_Vo_t0[0]
    if new_Qc_t0!=ar_Qc_t0[0]:
        print 'Change pVc_t0'
        ar_Qc_t0=ar_Qc_t0*0.+new_Qc_t0

    if channel_lesotho:
        ar_lambda=channel_cell_tunnel(ar_cell_label,ar_coorx,ar_coory,ar_lambda,ar_cell_down)
        
    #~~~~~~Write parameter file~~~~~~#
    tab_param=np.zeros((len(ar_cell_label),nb_param))
    tab_param[:,0]=ar_cell_label
    tab_param[:,1]=ar_coorx
    tab_param[:,2]=ar_coory
    tab_param[:,3]=ar_lambda
    tab_param[:,4]=ar_dam
    tab_param[:,5]=ar_tan_beta
    tab_param[:,6]=ar_L
    tab_param[:,7]=ar_Ks
    tab_param[:,8]=ar_theta_r
    tab_param[:,9]=ar_theta_s
    tab_param[:,10]=ar_n_o
    tab_param[:,11]=ar_n_c
    tab_param[:,12]=ar_cell_down
    tab_param[:,13]=ar_pVs_t0
    tab_param[:,14]=ar_Vo_t0
    tab_param[:,15]=ar_Qc_t0
    tab_param[:,16]=ar_kc

    f = file(file_out, 'w')
    io.write_array(f, tab_param)
    f.close()

def from_param_to_new_param_catchVsi(file_in,file_in_global,file_out,file_h5,file_catchment_saturation,fac_L=1.,fac_Ks=1.,fac_n_o=1.,fac_n_c=1.,\
                                    new_pVs_t0=-99.,new_pVc_t0=-99.,channel_lesotho=False):
    """

    * Objective
        Modifies the param file by multiplying some variables (L, Ks, n_o,n_c)
        or replacing by a given new value for the initial level of reservoir in percent
        The initial soil saturation is affected from the evolutive saturation file
    * Input
      - file_in: parameter file to be modified
      - the other parameters are assigned by default and have to be changed according to user's choice.
    * Output
      - file_out: new parameter file
    * Comment
      The boolean channel_lesotho calls a subroutine that was used in the case of the Liebenbergsvlei catchment
      to link the tunnel output to the river by transforming the initially non-channel cells into channel cells
    """
    nb_param=17.

    
    #--Read and compute the parameters to have the values of parameter and ar_Vsm
    #~~~~Read Global parameters file
    print 'Reading parameter file'
    X,Dt,alpha_s,alpha_o,alpha_c,A_thres,W_min,W_max\
      =pm.read_global_parameters(file_in_global)
    #~~~~Read Cell parameters file
    ar_cell_label,ar_coorx,ar_coory,ar_lambda,ar_dam,ar_tan_beta,ar_L,ar_Ks,\
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
    #~~~~Computation of model parameters from physical parameters
    ar_Vsm, ar_b_s, ar_b_o, ar_W, ar_b_c\
      =pm.compute_cell_param(X,Dt,alpha_s,alpha_o,alpha_c,nb_cell,\
                              A_thres,W_max,W_min,\
                              ar_lambda,ar_tan_beta,ar_L,\
                              ar_Ks,ar_theta_r,ar_theta_s,ar_n_o,ar_n_c,\
                              ar_A_drained)

    #Read the soil volume file
    ndar_Vs=np.array(ut.read_one_array_hdf(file_h5,'/Soil/','V_s'))

    #Read the file of catchment saturation rates
    f=file(file_catchment_saturation,'r')
    tab_rate=io.read_array(f)
    f.close()

    #Look for the rate closest to the expected mean value
    if new_pVs_t0>tab_rate[0]:
        ind=0
        print 'rate=100%'
    elif new_pVs_t0<tab_rate[-1]:
        ind=-1
        print 'rate=0%'
    else:
        loop=True
        i=-1
        while loop:
            i=i+1
            if new_pVs_t0>=tab_rate[i] and new_pVs_t0>tab_rate[i+1]:
                ind=i+1
                loop=False
                print 'rate=',tab_rate[i]
    print ind
    ar_Vs=ndar_Vs[ind,:]
    ar_Vsi=ar_Vs/ar_Vsm*100.
    print ar_Vsi
    #~~~~~~Change in parameters~~~~~~#
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
    if new_pVc_t0>=0.:
        print 'Change pVc_t0'
        ar_pVc_t0=ar_pVc_t0*0.+new_pVc_t0
    if new_pVs_t0>=0.:
        print 'Change pVs_t0'
        ar_pVs_t0=ar_Vsi
        print ar_pVs_t0[0]

    if channel_lesotho:
        ar_lambda,ar_n_c=channel_cell_tunnel(ar_cell_label,ar_coorx,ar_coory,ar_lambda,ar_cell_down,ar_n_c)
    
    #~~~~~~Write parameter file~~~~~~#
    tab_param=np.zeros((len(ar_cell_label),nb_param))
    tab_param[:,0]=ar_cell_label
    tab_param[:,1]=ar_coorx
    tab_param[:,2]=ar_coory
    tab_param[:,3]=ar_lambda
    tab_param[:,4]=ar_dam
    tab_param[:,5]=ar_tan_beta
    tab_param[:,6]=ar_L
    tab_param[:,7]=ar_Ks
    tab_param[:,8]=ar_theta_r
    tab_param[:,9]=ar_theta_s
    tab_param[:,10]=ar_n_o
    tab_param[:,11]=ar_n_c
    tab_param[:,12]=ar_cell_down
    tab_param[:,13]=ar_pVs_t0
    tab_param[:,14]=ar_Vo_t0
    tab_param[:,15]=ar_Qc_t0
    tab_param[:,16]=ar_kc

    f = file(file_out, 'w')
    io.write_array(f, tab_param)
    f.close()

def from_param_to_new_param_initial_simu(file_in,file_in_global,file_out,file_h5,indice_extrac_file,fac_L=1.,fac_Ks=1.,fac_n_o=1.,fac_n_c=1.,channel_lesotho=False):
    """
    * Objective
        Modifies the param file by multiplying some variables (L, Ks, n_o,n_c)
        and/or by replacing the initial level of reservoirs by values extracted from a TOPKAPI simulation file.
    * Input
      - file_in: parameter file to be modified
      - the other parameters are assigned by default and have to be changed according to user's choice.
    * Output
      - file_out: new parameter file
    * Comment
      The boolean channel_lesotho calls a subroutine that was used in the case of the Liebenbergsvlei catchment
      to link the tunnel output to the river by transforming the initially non-channel cells into channel cells
    """
    nb_param=17.

    #--Read and compute the parameters to have the values of parameter and ar_Vsm
    #~~~~Read Global parameters file
    print 'Reading parameter file'
    X,Dt,alpha_s,alpha_o,alpha_c,A_thres,W_min,W_max\
      =pm.read_global_parameters(file_in_global)
    #~~~~Read Cell parameters file
    ar_cell_label,ar_coorx,ar_coory,ar_lambda,ar_dam,ar_tan_beta,ar_L,ar_Ks,\
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
    #~~~~Computation of model parameters from physical parameters
    ar_Vsm, ar_b_s, ar_b_o, ar_W, ar_b_c\
      =pm.compute_cell_param(X,Dt,alpha_s,alpha_o,alpha_c,nb_cell,\
                              A_thres,W_max,W_min,\
                              ar_lambda,ar_tan_beta,ar_L,\
                              ar_Ks,ar_theta_r,ar_theta_s,ar_n_o,ar_n_c,\
                              ar_A_drained)

    #Read the soil volume file
    ndar_Vs=np.array(ut.read_one_array_hdf(file_h5,'/Soil/','V_s'))
    #Read the overland volume file
    ndar_Vo=np.array(ut.read_one_array_hdf(file_h5,'/Overland/','V_o'))
    #Read the channel dischargefile
    ndar_Qc=np.array(ut.read_one_array_hdf(file_h5,'/Channel/','Qc_out'))


    ar_Vs=ndar_Vs[indice_extrac_file,:]
    ar_pVs_t0=ar_Vs/ar_Vsm*100.
    ar_Vo_t0=ndar_Vo[indice_extrac_file,:]
    ar_Qc_t0=ndar_Qc[indice_extrac_file,:]
    #~~~~~~Change in parameters~~~~~~#
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
        
    if channel_lesotho:
        ar_lambda,ar_n_c=channel_cell_tunnel(ar_cell_label,ar_coorx,ar_coory,ar_lambda,ar_cell_down,ar_n_c)
    
    #~~~~~~Write parameter file~~~~~~#
    tab_param=np.zeros((len(ar_cell_label),nb_param))
    tab_param[:,0]=ar_cell_label
    tab_param[:,1]=ar_coorx
    tab_param[:,2]=ar_coory
    tab_param[:,3]=ar_lambda
    tab_param[:,4]=ar_dam
    tab_param[:,5]=ar_tan_beta
    tab_param[:,6]=ar_L
    tab_param[:,7]=ar_Ks
    tab_param[:,8]=ar_theta_r
    tab_param[:,9]=ar_theta_s
    tab_param[:,10]=ar_n_o
    tab_param[:,11]=ar_n_c
    tab_param[:,12]=ar_cell_down
    tab_param[:,13]=ar_pVs_t0
    tab_param[:,14]=ar_Vo_t0
    tab_param[:,15]=ar_Qc_t0
    tab_param[:,16]=ar_kc

    f = file(file_out, 'w')
    io.write_array(f, tab_param)
    f.close()


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

def distance(x1,y1,x2,y2):
    """
    Compute the distance between two points
    """
    dist=((x1-x2)**2+(y1-y2)**2)**0.5
    return dist

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

def channel_cell_tunnel(ar_cell_label,ar_coorx,ar_coory,ar_lambda,ar_cell_down,ar_n_c,Xoutlet=-255032.83,Youtlet=-3149857.34):
    """
    Change the ar_lambda to connnect the tunnel to the channel. Zero values are changed in 1 values.
    The values of the manning coefficient for the changed cells is taken as the manning value of the closest channel cell.
    """
    cell=find_cell_coordinates(ar_cell_label,Xoutlet,Youtlet,ar_coorx,ar_coory,ar_lambda,channel=False)
    hillslope=True
    li_ind=[]
    while hillslope:
        ind=np.where(ar_cell_label==cell)
        if ar_lambda[ind]==1.:
            hillslope=False
            last_ind=ind
        else:
            li_ind.append(ind)
            print cell,ind
            ar_lambda[ind]=1.
            cell=ar_cell_down[ind]
    for i in li_ind:
        ar_n_c[i]=ar_n_c[last_ind]
    
    return ar_lambda,ar_n_c
        

if __name__ == '__main__':


##    #### TO EXTRACT A SUBCATCHMENT ####
##    path_main='C:/Theo/liebenbergsvlei/topkapi_model/parameters/parameter_files_Aug07/'
##    file_in=path_main+'cell_param_Aug07_Vsi80%.dat'
##    file_out=path_main+'cell_param_sub_Aug07_Vsi80%.dat'
##    Xoutlet=-259210.38
##    Youtlet=-3065828.87
##    nom_image='subC8H020'
##    from_param_to_subcatch_param(file_in, file_out,Xoutlet,Youtlet,nom_image)
    

##    #### TO CHANGE THE INITIAL SOIL MOISTURE (CONSTANT) ####
##    path_main='C:/Theo/liebenbergsvlei/topkapi_model/parameters/'
##    file_in=path_main+'cell_param_20Feb07_Vsi80%.dat'
##    file_out=path_main+'cell_param_Vsi100%.dat'
##    from_param_to_new_param(file_in,file_out,new_pVs_t0=100.)


##    #### TO CHANGE THE INITIAL SOIL MOISTURE (VARIABLE) ####
##    path_main='C:/Theo/liebenbergsvlei/topkapi_model/parameters/parameter_files_Aug07/'
##    file_in=path_main+'cell_param_Aug07_Vsi80%.dat'
##    file_in_global=path_main+'global_param_6h.dat'
##    file_h5='C:/Theo/liebenbergsvlei/topkapi_model/field_initial_soil_moisture/Event1_6h_L1.1_Ks101.0_no1.0_nc1.0_100%_zerorain.h5'
##    file_catchment_saturation='C:/Theo/liebenbergsvlei/topkapi_model/field_initial_soil_moisture/Catchment_soil_moisture_evolution.dat'
##
##    tab=[10,20,30,40,50,60,70,80,90,100]
##    for i in tab:
##        print ''
##        print i
##        file_out=path_main+'cell_param_Aug07_catchVsi'+str(i)+'%.dat'
##        from_param_to_new_param_catchVsi(file_in,file_in_global,file_out,file_h5,file_catchment_saturation,new_pVs_t0=i,channel_lesotho=True)


##    #### TO MODIFY INITIAL VALUES FROM A SIMULATION FILE ####
##    path_main='C:/Theo/liebenbergsvlei/topkapi_model/parameters/parameter_files_Aug07/'
##    file_in=path_main+'cell_param_sub_Aug07_Vsi80%.dat'
##    file_in_global='C:/Theo/liebenbergsvlei/topkapi_model/parameters/global_param_6h.dat'
##    file_out=path_main+'cell_param_sub_Aug07_intialized.dat'
##    file_h5='C:/Theo/liebenbergsvlei/topkapi_model/simulation_after_calibration/simulations/Event3_6h_subC8H020_L1.1_Ks101.0_no1.0_nc7.0_40%.h5'
##    indice_extrac_file=180
##    from_param_to_new_param_initial_simu(file_in,file_in_global,file_out,file_h5,indice_extrac_file,fac_L=1.,fac_Ks=1.,fac_n_o=1.,fac_n_c=1.,channel_lesotho=True)

    

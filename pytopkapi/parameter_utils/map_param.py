from TOPKAPI import pretreatment as pm
from TOPKAPI import utils as ut
import numpy as np

def field_map(ar_field,ar_coorx,ar_coory,X,image_out,title,flip=0,min_val=0.,max_val=0.):

    import pylab as pl
    import matplotlib.numerix.ma as M

    #max_val=max(ar_field)
    
    xmin=min(ar_coorx);xmax=max(ar_coorx)
    ymin=min(ar_coory);ymax=max(ar_coory)
    step=X
    nx=(xmax-xmin)/step+1
    ny=(ymax-ymin)/step+1
    
    ar_indx=np.array((ar_coorx-xmin)/step,int)
    ar_indy=np.array((ar_coory-ymin)/step,int)
    
    ar_map=np.ones((ny,nx))*-99.9
    ar_map[ar_indy,ar_indx]=ar_field
    
    if flip==1:
        ar_map=np.flipud(ar_map)
        
    ar_map2 = M.masked_where(ar_map <0, ar_map)

        
    ut.check_file_exist(image_out)

    pl.clf()
    pl.imshow(ar_map2,interpolation='Nearest',origin='lower',vmax=max_val,vmin=min_val)
    pl.title(title)
    pl.colorbar()
    pl.savefig(image_out)

def field_map2(ar_field,ar_coorx,ar_coory,X,image_out,title,flip=0,min_val=0.,max_val=0.):

    import pylab as pl
    import matplotlib.numerix.ma as M

    #max_val=max(ar_field)
    
    xmin=min(ar_coorx);xmax=max(ar_coorx)
    ymin=min(ar_coory);ymax=max(ar_coory)
    step=X
    nx=(xmax-xmin)/step+1
    ny=(ymax-ymin)/step+1
    
    ar_indx=np.array((ar_coorx-xmin)/step,int)
    ar_indy=np.array((ar_coory-ymin)/step,int)
    
    ar_map=np.ones((ny,nx))*-99.9
    ar_map[ar_indy,ar_indx]=ar_field
    
    if flip==1:
        ar_map=np.flipud(ar_map)
        
    ar_map2=ar_map

        
    ut.check_file_exist(image_out)

    pl.clf()
    pl.imshow(ar_map2,interpolation='Nearest',origin='lower',vmax=max_val,vmin=min_val)
    pl.title(title)
    pl.colorbar()
    pl.savefig(image_out)


#if __name__ == '__main__':
##~~~~~~~~~~~ INPUTS FILES ~~~~~~~~~~~##
file_cell_param='E:/post_doc/liebenbergsvlei/topkapi_model_Oct07/parameters/cell_parameter/8D/slope_GIS_channel/cell_param_Vsi100.0_slope.dat'
file_global_param='E:/post_doc/liebenbergsvlei/topkapi_model_Oct07/parameters/global_parameter/global_param_Wmin5.0_Wmax40.0.dat'
path_out='E:/post_doc/liebenbergsvlei/topkapi_model_Oct07/parameters/cell_parameter/8D/slope_GIS_channel/'

#~~~~Read Global parameters file
X,Dt,alpha_s,alpha_o,alpha_c,A_thres,W_min,W_max\
  =pm.read_global_parameters(file_global_param)
#~~~~Read Cell parameters file
ar_cell_label,ar_coorx,ar_coory,ar_lambda,ar_Xc,ar_dam,ar_tan_beta,ar_tan_beta_channel,ar_L,ar_Ks,\
ar_theta_r,ar_theta_s,ar_n_o,ar_n_c,\
ar_cell_down,ar_pVs_t0,ar_Vo_t0,ar_Qc_t0,ar_kc\
    =pm.read_cell_parameters(file_cell_param)

#Lambda
image_out=path_out+'field_lambda.png'
field_map(ar_lambda,ar_coorx,ar_coory,X,image_out,'Channel cells',max_val=max(ar_lambda))
#ar_Xc
image_out=path_out+'field_Xc.png'
field_map(ar_Xc,ar_coorx,ar_coory,X,image_out,'Channel length',max_val=max(ar_Xc))
#Ks
image_out=path_out+'field_KS.png'
field_map(ar_Ks,ar_coorx,ar_coory,X,image_out,'Ks',max_val=max(ar_Ks))
#Slopes
image_out=path_out+'field_tan_beta.png'
field_map(ar_tan_beta,ar_coorx,ar_coory,X,image_out,r'$Slope \ tan(\beta)$',max_val=max(ar_tan_beta))
#Slopes zero
image_out=path_out+'field_tan_beta_zero.png'
tab=ar_tan_beta*1
tab[tab==0.]=-10
tab[tab>0]=10
field_map(tab,ar_coorx,ar_coory,X,image_out,r'$Slope \ tan(\beta) \ zero$',min_val=min(tab),max_val=max(tab))
#Slopes
image_out=path_out+'field_tan_beta_channel.png'
field_map(ar_tan_beta_channel,ar_coorx,ar_coory,X,image_out,r'$Slope \ tan(\beta)$',max_val=max(ar_tan_beta))
#Slopes zero
image_out=path_out+'field_tan_beta_channelzero.png'
tab=ar_tan_beta_channel*1
tab[tab==0.]=-10
tab[tab>0]=10
field_map(tab,ar_coorx,ar_coory,X,image_out,r'$Slope \ tan(\beta) \ zero$',min_val=min(tab),max_val=max(tab))
#diff slope
image_out=path_out+'field_diff_slope_GIS-channel.png'
tab=ar_tan_beta-ar_tan_beta_channel
field_map2(tab,ar_coorx,ar_coory,X,image_out,r'$Diff tan(\beta) \ GIS-Channel$',min_val=min(tab),max_val=max(tab))
#diff slope
image_out=path_out+'field_diff_slope_GIS-channel2.png'
tab=ar_tan_beta-ar_tan_beta_channel
tab2=tab*1.
tab2[tab<-10]=-10.
tab2[(tab>-10) & (tab<0)]=1.
tab2[tab>=0]=2.
##tab2[tab<0]=1.
##tab2[tab>0]=2.
field_map(tab2,ar_coorx,ar_coory,X,image_out,r'$Diff tan(\beta) \ GIS-Channel$',min_val=min(tab2),max_val=max(tab2))

#no
image_out=path_out+'field_no.png'
field_map(ar_n_o,ar_coorx,ar_coory,X,image_out,'n overland',max_val=max(ar_n_o))
#nc
image_out=path_out+'field_nc.png'
field_map(ar_n_c,ar_coorx,ar_coory,X,image_out,'n channel',max_val=max(ar_n_c))
#L
image_out=path_out+'field_L.png'
field_map(ar_L,ar_coorx,ar_coory,X,image_out,'Soil depth',max_val=max(ar_L))
#Theta_s
image_out=path_out+'field_theta_s.png'
field_map(ar_theta_s,ar_coorx,ar_coory,X,image_out,'Saturated soil moisture',max_val=max(ar_theta_s))
#Theta_r
image_out=path_out+'field_theta_r.png'
field_map(ar_theta_r,ar_coorx,ar_coory,X,image_out,'Residual soil moisture',max_val=max(ar_theta_r))

#Vsi
image_out=path_out+'field_pVs_t0.png'
field_map(ar_pVs_t0,ar_coorx,ar_coory,X,image_out,'Initial soil moisture %',max_val=max(ar_pVs_t0))

#label
image_out=path_out+'field_label.png'
field_map(ar_cell_label,ar_coorx,ar_coory,X,image_out,'Cell label',max_val=max(ar_cell_label))

    

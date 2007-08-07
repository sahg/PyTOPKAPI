
""" result_analysis.py

Analysis of outputs of TOPKAPI model
"""

__author__ = "Theo Vischel"
__version__ = "$Revision: 1.0 $"
__date__ = "$Date: 13/10/2006 $"

import numpy as np
import scipy as sp
import pylab as pl
import utils as ut
import scipy.io as io
import tables as h5
import pretreatment as pm
import datetime as dt
from matplotlib.dates import drange,date2num,num2date,YearLocator, MonthLocator,DayLocator, HourLocator,DateFormatter

#####################################
##~~ Compute all the upper cells ~~##
#####################################
#````````````````````````````````````````````````````
def direct_up_cell(ar_cell,ar_cell_down,ar_cell_label):
    ar_out=sp.array([],int)
    for i in range(len(ar_cell)):
        ind=np.where(ar_cell_down==ar_cell[i])
        if np.size(ar_out)==0:
            ar_out=ar_cell_label[ind]
        else:
            if np.size(ar_cell_label[ind])>0:
                ar_out=np.concatenate((ar_out,ar_cell_label[ind]))
    return ar_out

def all_up_cell(cell,ar_cell_down,ar_cell_label):
    """ all_up_cell
        Return an array of the label of all the cells drained by a given cell (cell=label of cell)
        Example b=all_up_cell(13,ar_cell_down,ar_cell_label)
                b = array([13,  8, 14,  4,  9,  2,  5])
    """
    a=sp.ones(1,int)*cell
    b=a
    while np.size(a)!=0:
        a=direct_up_cell(a,ar_cell_down,ar_cell_label)
        b=np.concatenate((b,a))
    return b

##################################
##~~ Read outputs of TOPKAPI ~~~##
##################################
#````````````````````````````````````````````````````
def read_topkapi_output(file_in):

    h5file_in=h5.openFile(file_in,mode='r')
    
    group='/Soil/'
    node = h5file_in.getNode(group+'Qin')
    ndar_a_s=node.read()
    node = h5file_in.getNode(group+'Volume')
    ndar_Vs=node.read()
    node = h5file_in.getNode(group+'Qout')
    ndar_Qs_out=node.read()
    
    group='/Overland/'
    node = h5file_in.getNode(group+'Qin')
    ndar_a_o=node.read()
    node = h5file_in.getNode(group+'Volume')
    ndar_Vo=node.read()
    node = h5file_in.getNode(group+'Qout')
    ndar_Qo_out=node.read()
    
    group='/Channel/'
    node = h5file_in.getNode(group+'Qin')
    ndar_a_c=node.read()
    node = h5file_in.getNode(group+'Volume')
    ndar_Vc=node.read()
    node = h5file_in.getNode(group+'Qout')
    ndar_Qc_out=node.read()
    
    group='/Other_variables/'
    node = h5file_in.getNode(group+'Q_to_next_cell')
    ndar_Q_to_next_cell=node.read()
    node = h5file_in.getNode(group+'Qc_cell_up')
    ndar_Qc_cell_up=node.read()

    group='/Evap/'
    node = h5file_in.getNode(group+'ETa')
    ndar_ETa=node.read()

    h5file_in.close()
    
    return ndar_a_s,ndar_Vs,ndar_Qs_out,\
           ndar_a_o,ndar_Vo,ndar_Qo_out,\
           ndar_a_c,ndar_Vc,ndar_Qc_out,\
           ndar_Q_to_next_cell,\
           ndar_Qc_cell_up,ndar_ETa

##################################
##~~~ Continuity validation ~~~~##
##################################

##---     Main programm    ----##
def continuity(file_in,path_out,cell,ndar_rain,ar_lambda,X,Dt,\
               ar_label_sort,ar_cell_label,ar_cell_down):
    import tables as h5
    
    #Read of data from the outputs of TOPKAPI in hdf5 format
    ndar_a_s,ndar_Vs,ndar_Qs_out,\
    ndar_a_o,ndar_Vo,ndar_Qo_out,\
    ndar_a_c,ndar_Vc,ndar_Qc_out,\
    ndar_Q_to_next_cell,\
    ndar_Qc_cell_up,ndar_ETa=read_topkapi_output(file_in)
    
    
    #Continuity analysis
    h5file=h5.openFile(path_out+'continuity.h5',mode='w',title='Numerical_continuity')

    ## 1.At each time step
    group_et=h5file.createGroup('/','Each_time_step','Each time step')
    ## 2.Over all the simulation
    group_os=h5file.createGroup('/','Overall_the_simulation','Overall_the_simulation')

    #At the cell scale
    group_ts=h5file.createGroup(group_et,'Cells','Cells')
    group_as=h5file.createGroup(group_os,'Cells','Cells')

    #For soil only
    ndar_In=ndar_a_s
    ndar_Out=ndar_a_o+ndar_Qs_out
    ndar_Store=ndar_Vs
    #1.
    ndar_Vol_err_cell_s,ndar_Rel_err_cell_s\
                =cell_continuity_each_time_step(ndar_In,ndar_Out,ndar_Store,Dt)
    h5file.createArray(group_ts,'Vol_err_s',ndar_Vol_err_cell_s,'Vol_err_cell_s[time,cell]')
    h5file.createArray(group_ts,'Rel_err_s',ndar_Rel_err_cell_s,'Rel_err_cell_s[time,cell]')
    #2.
    ar_Vol_err_cell_s,ar_Rel_err_cell_s\
                =cell_continuity_simulation(ndar_In,ndar_Out,ndar_Store,Dt)
    h5file.createArray(group_as,'Vol_err_s',ar_Vol_err_cell_s,'Vol_err_cell_s[time,cell]')
    h5file.createArray(group_as,'Rel_err_s',ar_Rel_err_cell_s,'Rel_err_cell_s[time,cell]')


    #For overland only
    ndar_In=ndar_a_o
    ndar_Out=ndar_Qo_out
    ndar_Store=ndar_Vo
    #1.
    ndar_Vol_err_cell_o,ndar_Rel_err_cell_o\
                =cell_continuity_each_time_step(ndar_In,ndar_Out,ndar_Store,Dt)
    h5file.createArray(group_ts,'Vol_err_o',ndar_Vol_err_cell_o,'Vol_err_cell_o[time,cell]')
    h5file.createArray(group_ts,'Rel_err_o',ndar_Rel_err_cell_o,'Rel_err_cell_o[time,cell]')
    #2.
    ar_Vol_err_cell_o,ar_Rel_err_cell_o\
                =cell_continuity_simulation(ndar_In,ndar_Out,ndar_Store,Dt)
    h5file.createArray(group_as,'Vol_err_o',ar_Vol_err_cell_o,'Vol_err_cell_o[time,cell]')
    h5file.createArray(group_as,'Rel_err_o',ar_Rel_err_cell_o,'Rel_err_cell_o[time,cell]')

    #For Channel only
    ndar_In=ndar_a_c
    ndar_Out=ndar_Qc_out
    ndar_Store=ndar_Vc
    #1.
    ndar_Vol_err_cell_c,ndar_Rel_err_cell_c\
                =cell_continuity_each_time_step(ndar_In,ndar_Out,ndar_Store,Dt)
    h5file.createArray(group_ts,'Vol_err_c',ndar_Vol_err_cell_c,'Vol_err_cell_c[time,cell]')
    h5file.createArray(group_ts,'Rel_err_c',ndar_Rel_err_cell_c,'Rel_err_cell_c[time,cell]')
    #2.
    ar_Vol_err_cell_c,ar_Rel_err_cell_c\
                =cell_continuity_simulation(ndar_In,ndar_Out,ndar_Store,Dt)
    h5file.createArray(group_as,'Vol_err_c',ar_Vol_err_cell_c,'Vol_err_cell_c[time,cell]')
    h5file.createArray(group_as,'Rel_err_c',ar_Rel_err_cell_c,'Rel_err_cell_c[time,cell]')

    #For Soil and Overland
    ndar_In=ndar_a_s
    ndar_Out=ndar_Qs_out+ndar_Qo_out+ndar_ETa*1e3/Dt
    ndar_Store=ndar_Vs+ndar_Vo
    #1.
    ndar_Vol_err_cell_s_o,ndar_Rel_err_cell_s_o\
                =cell_continuity_each_time_step(ndar_In,ndar_Out,ndar_Store,Dt)
    h5file.createArray(group_ts,'Vol_err_s_o',ndar_Vol_err_cell_s_o,'Vol_err_cell_s_o[time,cell]')
    h5file.createArray(group_ts,'Rel_err_s_o',ndar_Rel_err_cell_s_o,'Rel_err_cell_s_o[time,cell]')
    #2.
    ar_Vol_err_cell_s_o,ar_Rel_err_cell_s_o\
                =cell_continuity_simulation(ndar_In,ndar_Out,ndar_Store,Dt)
    h5file.createArray(group_as,'Vol_err_s_o',ar_Vol_err_cell_s_o,'Vol_err_cell_s_o[time,cell]')
    h5file.createArray(group_as,'Rel_err_s_o',ar_Rel_err_cell_s_o,'Rel_err_cell_s_o[time,cell]')

    #For the global cell (soil+overland+channel)
    ndar_In=ndar_a_s+ndar_Qc_cell_up
    ndar_Out=ndar_Qc_out+ndar_Q_to_next_cell+ndar_ETa*1e3/Dt
    ndar_Store=ndar_Vs+ndar_Vo+ndar_Vc
    #1.
    ndar_Vol_err_cell,ndar_Rel_err_cell\
                =cell_continuity_each_time_step(ndar_In,ndar_Out,ndar_Store,Dt)
    h5file.createArray(group_ts,'Vol_err_cell',ndar_Vol_err_cell,'Vol_err_cell[time,cell]')
    h5file.createArray(group_ts,'Rel_err_cell',ndar_Rel_err_cell,'Rel_err_cell[time,cell]')
    #2.
    ar_Vol_err_cell,ar_Rel_err_cell\
                =cell_continuity_simulation(ndar_In,ndar_Out,ndar_Store,Dt)
    h5file.createArray(group_as,'Vol_err_cell',ar_Vol_err_cell,'Vol_err_cell[time,cell]')
    h5file.createArray(group_as,'Rel_err_cell',ar_Rel_err_cell,'Rel_err_cell[time,cell]')

   
    print ' '
    print 'Max continuity cell error in volume (m3)'
    print 'Cell     : ',max(abs(ar_Vol_err_cell))
    
    print 'Max continuity cell relative error(%)'
    print 'Cell     : ',max(abs(ar_Rel_err_cell))
    
    ## At catchment scale

    #At each time step
    group=h5file.createGroup(group_et,'Catchment','Catchment array')
    ar_Vol_err_catchment,ar_Rel_err_catchment=catchment_continuity_each_time_step(ndar_rain, ndar_ETa, \
                                        ndar_Vs, ndar_Vo, ndar_Vc,\
                                        ndar_Q_to_next_cell, ndar_Qc_out,\
                                        Dt, X, ar_label_sort)
    h5file.createArray(group,'Rel_err_catchment',ar_Rel_err_catchment,'ar_Rel_err_catchment')
    h5file.createArray(group,'Vol_err_catchment',ar_Vol_err_catchment,'ar_Vol_err_catchment')

    
    #over all the simulation period
    group=h5file.createGroup('/Overall_the_simulation','Catchment','Catchment_array')

    li_Vol_err_catchment=[]
    li_Rel_err_catchment=[]
    nb_time_step=len(ndar_rain[:,0])
    for t in np.arange(nb_time_step):
        Vol_err_catchment,Rel_err_catchment=\
                    catchment_continuity(ndar_rain, ndar_ETa, \
                                         ndar_Vs, ndar_Vo, ndar_Vc,\
                                         ndar_Q_to_next_cell, ndar_Qc_out,\
                                         Dt, X, t,\
                                         ar_label_sort)
        li_Vol_err_catchment.append(Vol_err_catchment)
        li_Rel_err_catchment.append(Rel_err_catchment)


    h5file.createArray(group,'Rel_err_catchment',Rel_err_catchment,'Rel_err_catchment')
    h5file.createArray(group,'Vol_err_catchment',Vol_err_catchment,'Vol_err_catchment')

    print 'Volume Continuity Error at catchment scale (%)'
    print ' %1.8f' % Vol_err_catchment
    print 'Relative Continuity Error at catchment scale (%)'
    print ' %1.8f' % Rel_err_catchment

    h5file.close()

    pl.clf()
    pl.plot(li_Rel_err_catchment,'k-')
    pl.xlabel('Time step')
    pl.ylabel('Relative Error catchment (%)')
    #pl.ylim(-1,1)
    pl.savefig(path_out+'Rel_err_catchment_continuity.png')

    pl.clf()
    pl.plot(li_Vol_err_catchment,'k-')
    pl.xlabel('Time step')
    pl.ylabel('Volume Error')
    pl.savefig(path_out+'Vol_err_catchment_continuity.png')

    pl.clf()
    pl.plot(ar_Vol_err_catchment,'k-')
    pl.xlabel('Time step')
    pl.ylabel('Volume Error')
    pl.savefig(path_out+'Vol_err_catchment_continuity_ets.png')


##---     Subroutines    ----##

#````````````````````````````````````````````````````
def cell_continuity_simulation(ndar_In,ndar_Out,ndar_Store,Dt):
    Input=sum(ndar_In)*Dt
    Output=sum(ndar_Out)*Dt
    Storage=ndar_Store[-1,:]-ndar_Store[0,:]

    ar_Vol_err=Output+Storage-Input
    ar_Rel_err=ar_Vol_err/Input
    
    return ar_Vol_err,ar_Rel_err

#````````````````````````````````````````````````````
def catchment_continuity(ndar_rain, ndar_ET, \
                         ndar_Vs, ndar_Vo, ndar_Vc,\
                         ndar_Q_to_next_cell, ndar_Qc_out,\
                         Dt, X, time_step,\
                         ar_label_sort):
    #ndar_rain and ndar_ET must given in mm/time_step
    #ndar_Q in m3/s
    #ndar_V in m3
    Input= sum(sum(ndar_rain[:time_step+1,0:]))*1e-3*X**2-sum(sum(ndar_ET[:time_step+1,0:]))*1e-3*X**2
    cell_outlet=ar_label_sort[-1]
    Output= sum(ndar_Q_to_next_cell[:time_step+1,cell_outlet] + \
                   ndar_Qc_out[:time_step+1,cell_outlet])*Dt

    Storage=  sum(ndar_Vs[time_step+1,:]-ndar_Vs[0,:])+\
              sum(ndar_Vo[time_step+1,:]-ndar_Vo[0,:])+\
              sum(ndar_Vc[time_step+1,:]-ndar_Vc[0,:])
    Vol_err = Output-Input+Storage   
    Rel_err = Vol_err/Input*100.
    
    return Vol_err,Rel_err
#````````````````````````````````````````````````````
def catchment_continuity_each_time_step(ndar_rain, ndar_ET, \
                                        ndar_Vs, ndar_Vo, ndar_Vc,\
                                        ndar_Q_to_next_cell, ndar_Qc_out,\
                                        Dt, X,ar_label_sort):
    #ndar_rain and ndar_ET must given in mm/time_step
    #ndar_Q in m3/s
    #ndar_V in m3
    Input= np.sum(ndar_rain,axis=1)*1e-3*X**2-np.sum(ndar_ET,axis=1)*1e-3*X**2

    cell_outlet=ar_label_sort[-1]
    Output= (ndar_Q_to_next_cell[:,cell_outlet] + \
                   ndar_Qc_out[:,cell_outlet])*Dt

    Storage=  np.sum(ndar_Vs[1:,:]-ndar_Vs[:-1,:],axis=1)+\
              np.sum(ndar_Vo[1:,:]-ndar_Vo[:-1,:],axis=1)+\
              np.sum(ndar_Vc[1:,:]-ndar_Vc[:-1,:],axis=1)

    ar_Vol_err = Output-Input+Storage   
    ar_Rel_err = ar_Vol_err/Input*100.
    
    return ar_Vol_err,ar_Rel_err

def cell_continuity_each_time_step(ndar_In,ndar_Out,ndar_Store,Dt):
    Input=ndar_In*Dt
    Output=ndar_Out*Dt
    Storage=ndar_Store[1:,:]-ndar_Store[:-1,:]
     
    ndar_Vol_err=Output+Storage-Input
    ndar_Rel_err=ndar_Vol_err/Input

    return ndar_Vol_err,ndar_Rel_err

#````````````````````````````````````````````````````
def subcatchment_continuity(ndar_rain, ndar_ET,\
                            ndar_Vs, ndar_Vo, ndar_Vc,\
                            ndar_Q_to_next_cell, ndar_Qc_out,\
                            Dt, X, time_step,\
                            cell,ar_cell_down,ar_cell_label):
    #ndar_rain and ndar_ET must given in mm/time_step
    #ndar_Q in m3/s
    #ndar_V in m3
    #selection of the catchment cells
    ind_up=all_up_cell(cell,ar_cell_down,ar_cell_label)
    
    if len(ind_up)>1:
        Input= sum(sum(ndar_rain[:time_step+1,ind_up]))*1e-3*X**2-sum(sum(ndar_ET[:time_step+1,ind_up]))*1e-3*X**2
        Output= sum(ndar_Q_to_next_cell[:time_step+1,cell]+\
                    ndar_Qc_out[:time_step+1,cell])*Dt

        Storage=  sum(ndar_Vs[time_step+1,ind_up]-ndar_Vs[0,ind_up])+\
                  sum(ndar_Vo[time_step+1,ind_up]-ndar_Vo[0,ind_up])+\
                  sum(ndar_Vc[time_step+1,ind_up]-ndar_Vc[0,ind_up])
        Vol_err = Output-Input+Storage   
        Rel_err = Vol_err/Input*100.
    else:
        Vol_err = 0. 
        Rel_err = 0.
    return Vol_err,Rel_err

##################################
##~~~   Arrays comparison   ~~~~##
##################################

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

###################################################
##~~ Suvsurface contribution to channel flows ~~~##
###################################################
def Q_subsurface_contribution(cell,t1,t2,\
                              ndar_Qc_out,ndar_a_c,\
                              ndar_Q_to_channel_sub,\
                              ind_up,li_cell_up,nb_cell):

    #Compute the ratio of channel flow that is due to subsurface flows
    li_ratio=[]
    for t in range(t1,t2):
        li_ratio.append(\
            ratio_subsurface_contribution(ndar_Qc_out[t,],\
                                          ndar_a_c[t,],\
                                          ndar_Q_to_channel_sub[t,],\
                                          ind_up,li_cell_up,nb_cell))

    ar_Qc_out=ndar_Qc_out[t1:t2,cell]
    ar_ratio=np.array(li_ratio)
    ar_Qc_sub=ar_ratio*ar_Qc_out

    return ar_Qc_sub

def ratio_subsurface_contribution(ar_Qc_out,ar_a_c,ar_Q_to_channel_sub,ind_up,\
                                  li_cell_up,nb_cell):
    '''
    Ratio computed at a given time step
    Dimension of array= [nb_cell] nb_cell= total number of upper cells
    '''
    ind_up_sort=ind_up[::-1]
    ar_ratio=np.ones(nb_cell)*-99.9
    for i in ind_up_sort:
        ind_up_local=li_cell_up[i]
        if ar_a_c[i] !=0:
            ar_ratio[i]=(sum(ar_ratio[ind_up_local]*ar_Qc_out[ind_up_local])+ar_Q_to_channel_sub[i])\
                         /ar_a_c[i]
        else:
            ar_ratio[i]=0.
        ratio_fin=ar_ratio[i]
    return ratio_fin



##################################
##~~ Rainfall-Runoff graphic ~~~##
##################################
def mean_spatial_rainfall(ndar_rain,ind_up,t1,t2):
    """
    Compute the average rainfall over a catchment delimited by the a list of cell labels
    """
    li_rain=[]
    for t in range(t1,t2):
        a=sum(ndar_rain[t,ind_up])
        li_rain.append(a)
    ar_rain=np.array(li_rain)/len(ind_up)
    return ar_rain

def rainfall_runoff_graphic(ndar_rain, ndar_Qc_out,ndar_Qc_sub,ar_Qobs,\
                            cell,t1,t2,date1,date2,image_out,\
                            ar_cell_down,ar_cell_label\
                            ,Xcell,Dt,mm=False,Qsub=False,Qobs=True):
    
    
    #Look for the cells drained by a given cell (cell=label)
    print '**Compute ind_up**'
    if cell==np.min(ar_cell_label):
        ind_up=ar_cell_label
        outlet=True
    else:
        ind_up=all_up_cell(cell,ar_cell_down,ar_cell_label)
        outlet=False

    #Compute the average rainfall over the such delimited catchment
    if outlet:
        ar_rain=np.average(ndar_rain[t1:t2,:],axis=1)
    else:
        ar_rain=np.average(ndar_rain[t1:t2,ind_up],axis=1)
        #ar_rain=mean_spatial_rainfall(ndar_rain,ind_up,t1,t2)

    ar_Qc_out=ndar_Qc_out[t1:t2,cell]
    if Qsub:
        ar_Qc_sub=ndar_Qc_sub[t1:t2,cell]
    
    
    #Graphic
    #X=np.arange(t1,t2)
 
    #conversion from m3/s to mm
    if mm:
        Y1=ar_Qc_out*(Dt/(1e-3*len(ind_up)*Xcell**2))
        Y3=ar_rain
        if Qsub:
            Y2=ar_Qc_sub*Dt/(1e-3*len(ind_up)*Xcell**2)
        if Qobs:
            Y4=ar_Qobs*Dt/(1e-3*len(ind_up)*Xcell**2)

        Qlabel=r'$Q \ (mm)$'
    else:
        Y1=ar_Qc_out
        if Qsub:
            Y2=ar_Qc_sub
        Y3=ar_rain
        if Qobs:
            Y4=ar_Qobs
        Qlabel=r'$Q \ (m^3/s)$'

    print '**plot the graph**'
    pl.clf()

    #Put the label as dates
    ax=pl.subplot(111)
    delta = dt.timedelta(hours=Dt)
    X= drange(date1, date2, delta)
    years    = YearLocator()
    months   = MonthLocator()
    days     = DayLocator()
    hours    = HourLocator()
    d=num2date(X[-1])-num2date(X[0])
    time_length=d.days
    lim_day0=3;lim_day1=61;lim_day2=3*365
    if time_length>lim_day1 and time_length<lim_day2:
        ax.xaxis.set_major_locator(months)
        Fmt = DateFormatter("%m/%y")
        ax.xaxis.set_major_formatter(Fmt)
        labels = ax.get_xticklabels()
        #pl.setp(labels,'rotation', 90)
    elif time_length>=lim_day2:
        ax.xaxis.set_major_locator(years)
        ax.xaxis.set_minor_locator(months)
        Fmt = DateFormatter("%m/%y")
        ax.xaxis.set_major_formatter(Fmt)
        labels = ax.get_xticklabels()
        pl.setp(labels,'rotation', 90)
    elif time_length>lim_day0 and time_length<lim_day1:
        ax.xaxis.set_major_locator(months)
        ax.xaxis.set_minor_locator(days)
        Fmt = DateFormatter("%d/%m/%y")
        ax.xaxis.set_major_formatter(Fmt)
    else:
        ax.xaxis.set_major_locator(hours)
        #ax.xaxis.set_minor_locator(hours)
        Fmt = DateFormatter("%H")
        ax.xaxis.set_major_formatter(Fmt)
 

    pl.ylabel(r'$Rainfall \  (mm)$',fontsize=18)
    #pl.bar(X,Y3,width=1,color='b',edgecolor='w')
    pl.plot(X,Y3,'b-')
    X_fill=np.concatenate((X[0:1],X,X[-1:]))
    Y3_fill=np.concatenate((sp.zeros(1),Y3,sp.zeros(1)))
    p=pl.fill(X_fill,Y3_fill,facecolor='b',edgecolor='b',alpha=0.2)
    pl.ylim(max(Y3)*2,min(Y3))
    pl.legend(['Rainfall'],'center left')
    leg = pl.gca().get_legend()
    leg.draw_frame(False)
    ltext  = leg.get_texts()
    pl.setp(ltext, fontsize=10)


    ax2 = pl.twinx()
    ax2.yaxis.tick_right()
    pl.ylabel(Qlabel,fontsize=18)
    if Qobs:
        pl.plot(X,Y4,'g-',linewidth=1.)
    pl.plot(X,Y1,'r-',linewidth=1.5)
    if Qsub:
        X_fill=np.concatenate((X[0:1],X,X[-1:]))
        Y2_fill=np.concatenate((sp.zeros(1),Y2,sp.zeros(1)))
        p=pl.fill(X_fill,Y2_fill,facecolor='r',edgecolor='r',alpha=0.2)
    


    pl.legend(['Q obs','Q sim'],'center right')
    leg = pl.gca().get_legend()
    leg.draw_frame(False)
    ltext  = leg.get_texts()
    pl.setp(ltext, fontsize=10)
      

    pl.savefig(image_out)


def plot_hydrograph(ndar_Qc_out,ndar_Qc_sub,cell,t1,t2,path_out):
    ut.check_folder_exist(path_out)

    Y1=ndar_Qc_out[t1:t2,cell]
    Y2=ndar_Qc_sub[t1:t2,cell]
    
    #Graphic
    X=np.arange(t1,t2)
    
    pl.clf()
    pl.title('Cell '+str(cell))
    pl.ylabel(r'$Q \ (m^3/s)$',fontsize=18)
    pl.xlabel('Time step',fontsize=18)
    pl.plot(X,Y1,'r-',linewidth=2)
    X_fill=np.concatenate((X[0:1],X,X[-1:]))
    Y2_fill=np.concatenate((sp.zeros(1),Y2,sp.zeros(1)))
    p=pl.fill(X_fill,Y2_fill,facecolor='r',edgecolor='r',alpha=0.2)

    pl.savefig(path_out+'Hydrograph_cell'+str(cell)+'.png')

##################################
##~~ Soil moisture graphic ~~~##
##################################
def mean_catchment_moisture(ndar_Vs,ar_Vsm,t1,t2,Dt,date1,date2,image_out,percentage=True):
    """
    Plot the mean catchment soil moisture from t1 to t2
    """
    if percentage:
        ndar_Vs_sat=ndar_Vs/ar_Vsm*100.
        Y=np.average(ndar_Vs_sat[t1:t2,:],axis=1)
        ytitle=r'$Mean catchment soil saturation \  (%)$'
    else:
        Y=np.average(ndar_Vs[t1:t2,:],axis=1)
        ytitle=r'$Mean catchment soil moisture \  (m^3)$'
    
    pl.clf()
    #Set up the date axes
    ax=pl.subplot(111)
    delta = dt.timedelta(hours=Dt)
    X= drange(date1, date2, delta)
    years    = YearLocator()
    months   = MonthLocator()
    days     = DayLocator()
    hours    = HourLocator()
    d=num2date(X[-1])-num2date(X[0])
    time_length=d.days
    lim_day0=3;lim_day1=61;lim_day2=3*365
    if time_length>lim_day1 and time_length<lim_day2:
        ax.xaxis.set_major_locator(months)
        Fmt = DateFormatter("%m/%y")
        ax.xaxis.set_major_formatter(Fmt)
        labels = ax.get_xticklabels()
        pl.setp(labels,'rotation', 90)
    elif time_length>=lim_day2:
        ax.xaxis.set_major_locator(years)
        ax.xaxis.set_minor_locator(months)
        Fmt = DateFormatter("%m/%y")
        ax.xaxis.set_major_formatter(Fmt)
        labels = ax.get_xticklabels()
        pl.setp(labels,'rotation', 90)
    elif time_length>lim_day0 and time_length<lim_day1:
        ax.xaxis.set_major_locator(months)
        ax.xaxis.set_minor_locator(days)
        Fmt = DateFormatter("%d/%m/%y")
        ax.xaxis.set_major_formatter(Fmt)
    else:
        ax.xaxis.set_major_locator(hours)
        #ax.xaxis.set_minor_locator(hours)
        Fmt = DateFormatter("%H")
        ax.xaxis.set_major_formatter(Fmt)

    #plot the graph
    pl.ylabel(ytitle,fontsize=18)
    pl.plot(X,Y,'k-')
    pl.ylim((0.,100.))
    pl.savefig(image_out)

    return Y

def cell_moisture_plot(ndar_Vs,ar_Vsm,tab_cell_label,t1,t2,Dt,date1,date2,image_out,percentage=True):
    """
    Plot the soil moisture from t1 to t2 for selected cells
    tab_cell_label is a table containing the label of a selection of cells to be considered
    """
    if percentage:
        ndar_Vs_sat=ndar_Vs/ar_Vsm*100.
        Y0=ndar_Vs_sat[t1:t2,:]
        nrow=np.shape(Y0)[0]
        ncol=len(tab_cell_label)
        Y=np.zeros((nrow,ncol))
        i=-1
        for cell in tab_cell_label:
            i=i+1
            Y[:,i]=Y0[:,cell]
        ytitle=r'$Mean catchment soil saturation \  (%)$'
    else:
        Y0=ndar_Vs[t1:t2,:]
        nrow=np.shape(Y0)[0]
        ncol=len(tab_cell_label)
        Y=np.zeros((nrow,ncol))
        i=-1
        for cell in tab_cell_label:
            i=i+1
            Y[:,i]=Y0[:,cell]
        ytitle=r'$Mean catchment soil moisture \  (m^3)$'

    tab_leg=np.str(tab_cell_label)
    tab_leg=tab_leg[1:-1].split()
    pl.clf()
    #Set up the date axes
    ax=pl.subplot(111)
    delta = dt.timedelta(hours=Dt)
    X= drange(date1, date2, delta)
    years    = YearLocator()
    months   = MonthLocator()
    days     = DayLocator()
    hours    = HourLocator()
    d=num2date(X[-1])-num2date(X[0])
    time_length=d.days
    lim_day0=3;lim_day1=61;lim_day2=3*365
    if time_length>lim_day1 and time_length<lim_day2:
        ax.xaxis.set_major_locator(months)
        Fmt = DateFormatter("%m/%y")
        ax.xaxis.set_major_formatter(Fmt)
        labels = ax.get_xticklabels()
        pl.setp(labels,'rotation', 90)
    elif time_length>=lim_day2:
        ax.xaxis.set_major_locator(years)
        ax.xaxis.set_minor_locator(months)
        Fmt = DateFormatter("%m/%y")
        ax.xaxis.set_major_formatter(Fmt)
        labels = ax.get_xticklabels()
        pl.setp(labels,'rotation', 90)
    elif time_length>lim_day0 and time_length<lim_day1:
        ax.xaxis.set_major_locator(months)
        ax.xaxis.set_minor_locator(days)
        Fmt = DateFormatter("%d/%m/%y")
        ax.xaxis.set_major_formatter(Fmt)
    else:
        ax.xaxis.set_major_locator(hours)
        #ax.xaxis.set_minor_locator(hours)
        Fmt = DateFormatter("%H")
        ax.xaxis.set_major_formatter(Fmt)

    #plot the graph
    for i in range(len(tab_cell_label)):
        pl.plot(X,Y[:,i])
    pl.ylabel(ytitle,fontsize=18)

    pl.legend(tab_leg,'center right')
    leg = pl.gca().get_legend()
    leg.draw_frame(False)
    ltext  = leg.get_texts()
    pl.setp(ltext, fontsize=10)
      
    pl.savefig(image_out)


def cell_SWI(ndar_Vs,ar_Vsm,tab_cell_label,t1,t2,Dt,date1,date2,file_out,fileh5=False):
    """
    Compute the SWI for each cell in tab_cell_label
    t0         22.3 21.2 18.5   ...  ...
    ..
    ..
    t2
    """
    ndar_Vs_sat=ndar_Vs/ar_Vsm*100.
    Y0=ndar_Vs_sat[t1:t2,:]
    nrow=np.shape(Y0)[0]
    ncol=len(tab_cell_label)
    Y=np.zeros((nrow,ncol))
    i=-1
    for cell in tab_cell_label:
        print cell
        i=i+1
        Y[:,i]=Y0[:,cell]

    #Write the outfile of soil moisture
    if h5:
        h5file=h5.openFile(file_out,mode='a',title='Cell results')
        group=h5file.createGroup('/','SWI','SWI arrays')
        h5file.createArray(group,'SWI',Y,'SWI array')
        h5file.close()
    else:
        f=file(file_out,'w')
        io.write_array(f,Y)
        f.close()

    

##################################
##~~~~~~~ Field maps ~~~~~~~~##
##################################
def field_map_ndar(ndar_field,t,ar_coorx,ar_coory,X,image_out):
    import matplotlib.numerix.ma as M

    ar_field=ndar_field[t,:]
    max_val=np.max(ndar_field)
    
    xmin=min(ar_coorx);xmax=max(ar_coorx)
    ymin=min(ar_coory);ymax=max(ar_coory)
    step=X
    nx=(xmax-xmin)/step+1
    ny=(ymax-ymin)/step+1
    
    ar_indx=np.array((ar_coorx-xmin)/step,int)
    ar_indy=np.array((ar_coory-ymin)/step,int)
    
    ar_map=sp.ones((ny,nx))*-99.9
    ar_map[ar_indy,ar_indx]=ar_field
    
    ar_map2 = M.masked_where(ar_map <0, ar_map)
    ut.check_file_exist(image_out)
    
    pl.clf()
    pl.imshow(ar_map2,interpolation='Nearest',origin='lower',vmax=max_val,vmin=0)
    pl.title('time step= '+ut.string(t,5))
    pl.colorbar()
    pl.savefig(image_out)
    
def field_map(ar_field,ar_coorx,ar_coory,X,image_out,title,flip=0):

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

        
    ut.check_file_exist(image_out)
    
    pl.clf()
    pl.imshow(ar_map2,interpolation='Nearest',origin='lower',vmax=max_val,vmin=0)
    pl.title(title)
    pl.colorbar()
    pl.savefig(image_out)

def field_map_ndar_stream_option(ndar_field,ndar_Q,ar_lambda,stream_option,t,ar_coorx,ar_coory,X,image_out):
    import matplotlib.numerix.ma as M

    ar_field=ndar_field[t,:]
    max_val=np.max(ndar_field)
    
    xmin=min(ar_coorx);xmax=max(ar_coorx)
    ymin=min(ar_coory);ymax=max(ar_coory)
    step=X
    nx=(xmax-xmin)/step+1
    ny=(ymax-ymin)/step+1
    
    ar_indx=np.array((ar_coorx-xmin)/step,int)
    ar_indy=np.array((ar_coory-ymin)/step,int)
    
    ar_map=sp.ones((ny,nx))*-99.9
    ar_map[ar_indy,ar_indx]=ar_field
    
    ar_map2 = M.masked_where(ar_map <0, ar_map)
    ut.check_file_exist(image_out)
    
    pl.clf()
    pl.imshow(ar_map2,interpolation='Nearest',origin='lower',vmax=max_val,vmin=0)
    pl.title('time step= '+ut.string(t,3))
    pl.colorbar()
    pl.savefig(image_out)


##################################
##~~       Aggregation       ~~~##
##################################
def agreg_Qout_hourly_to_delta(file_1h_Q,file_out_Q,cell,nb_time_step):
    """
    Nb_time_step is the number of hourly time steps to be agregated
    """
    #OUTPUT FILES
    file_h5='C:/Theo/liebenbergsvlei/topkapi_model/results/Event1/res_20Feb07_Vsi10%_alpha_s1.1.h5'

    #~~~~Read the output file
    ndar_Qc_out=ut.read_one_array_hdf(file_h5,'/Channel/','Qc_out')

    ar_Q=ndar_Qc_out[:,cell]
    len_Q=int(len(ar_Q)/nb_time_step)
    print len_Q
    tab=np.zeros((len_Q))
    for i in range(len_Q):
        tab[i]=np.average(ar_Q[i*nb_time_step:(i+1)*nb_time_step])

    f=file(file_out_Q,'w')
    io.write_array(f,tab)
    f.close()

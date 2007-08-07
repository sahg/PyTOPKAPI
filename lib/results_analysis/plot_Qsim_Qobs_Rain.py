import sys
sys.path.append('C:/Theo/topkapi_model/programme/TOPKAPI_Aug07/svn_working_copy/lib/')

import numpy as np
import scipy as sp
import pylab as pl
import datetime as dt
from matplotlib.dates import drange,date2num,num2date,YearLocator, MonthLocator,DayLocator, HourLocator,DateFormatter
import tables as h5

import utils as ut
reload(ut)


def graph(input):
    #############################################################################################
    path_fileQ='C:/Theo/topkapi_model/programme/TOPKAPI_Aug07/Example/run_the_model/results/'
    path_fileQobs='C:/Theo/topkapi_model/programme/TOPKAPI_Aug07/Example/analyse_the_results/'
    path_fileP='C:/Theo/topkapi_model/programme/TOPKAPI_Aug07/Example/run_the_model/forcing_variables/'
    path_out='C:/Theo/topkapi_model/programme/TOPKAPI_Aug07/Example/analyse_the_results/graphics/'
    #create path_out if it does'nt exist
    ut.check_folder_exist(path_out)

    tab_fileQ=['sample_event_simulation_results.h5']
    file_Qobs='C8H020_6h_Obs_Example.dat'
    fileP='rainfields.h5'
    event='sample_event'
    title=event

    Qobs=True
    Pobs=True
    nash=True
            
    tab_leg=['Model']

    image=path_out+tab_fileQ[0].replace('.h5','_Qsim_Qobs_P_color.png')

    #tab_col=['k','g','b','y','r','g','b','y','0.5']
    tab_col=['k','g','b','y','r','g','b','y','r']
    tab_style=['-','-','-','-','--','--','--','--','-']
    #tab_width=['1','2']
    tab_width=['1','1']
    #color_P='0.5'
    color_P='b'
    transparency_P=0.5#(0 for invisible)


    #Read the obs
    #Qobs
    ar_date,ar_Qobs=read_delta_flow(path_fileQobs+file_Qobs)
    date1=dt.datetime(ar_date[0,0],ar_date[0,1],ar_date[0,2],ar_date[0,3],ar_date[0,4])
    date1prim=dt.datetime(ar_date[1,0],ar_date[1,1],ar_date[1,2],ar_date[1,3],ar_date[1,4])
    date2=dt.datetime(ar_date[-1,0],ar_date[-1,1],ar_date[-1,2],ar_date[-1,3],ar_date[-1,4])
    delta=date1prim-date1
    X= drange(date1, date2+delta, delta)
    ndar_Q=np.zeros((len(X),len(tab_fileQ)))
    ar_nash=np.zeros(len(tab_fileQ))
    #Rain
    if Pobs:
        file_rain=path_fileP+fileP
        h5file_in=h5.openFile(file_rain,mode='r')
        group='/'+event+'/'
        node = h5file_in.getNode(group+'rainfall')
        ndar_rain=node.read()
        h5file_in.close()
        #Compute the mean catchment rainfall
        ar_rain=np.average(ndar_rain,axis=1)

    #Read the simulated data Q
    n=-1
    for file in tab_fileQ:
        n=n+1
        file_h5=path_fileQ+file
        ndar_Qc_out=ut.read_one_array_hdf(file_h5,'/Channel/','Qc_out')
        ndar_Q[:,n]=ndar_Qc_out[:,0]


    ##Graph
    pl.clf()
    #Set up the date axes
    lim_day0=3;lim_day1=61;lim_day2=3*365
    ax=pl.subplot(111)

    d=num2date(X[-1])-num2date(X[0])
    time_length=d.days
    years    = YearLocator()
    months   = MonthLocator()
    days     = DayLocator()
    ax.xaxis.set_major_locator(months)
    Fmt = DateFormatter("%m/%y")
    ax.xaxis.set_major_formatter(Fmt)
    labels = ax.get_xticklabels()
    #pl.setp(labels,'rotation', 90)
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
        #pl.setp(labels,'rotation', 90)
    elif time_length>lim_day0 and time_length<lim_day1:
        ax.xaxis.set_major_locator(months)
        ax.xaxis.set_minor_locator(days)
        Fmt = DateFormatter("%d/%m/%y")
        ax.xaxis.set_major_formatter(Fmt)
        #pl.setp(labels,'rotation', 90)
    else:
        ax.xaxis.set_major_locator(hours)
        ax.xaxis.set_minor_locator(hours)
        Fmt = DateFormatter("%H")
        ax.xaxis.set_major_formatter(Fmt)

    #plot the graph
    if Qobs:
        pl.plot(X,ar_Qobs,color=tab_col[-1],linestyle=tab_style[-1],linewidth=tab_width[-1])
        tab_leg.append(('Observation'))
        tab_leg=tab_leg[::-1]
        
    for n in range(len(tab_fileQ)):
        pl.plot(X,ndar_Q[:,n],color=tab_col[n],linestyle=tab_style[n],linewidth=tab_width[n])

    for n in range(len(tab_fileQ)):
        if nash:
            ar_nash[n]=ut.Nash(ndar_Q[:,n],ar_Qobs)
            pl.plot(X[0:1],ndar_Q[0:1,n],'w:')
            tab_leg.append(('Eff= '+str(ar_nash[n])[0:5]))
            
    pl.xlim(X[0],X[-1])
    ytitle=r'$Q \  (m^3/s)$'
    pl.ylabel(ytitle,fontsize=18)
    pl.title(title)
    pl.legend(tab_leg,'center')
    leg = pl.gca().get_legend()
    leg.draw_frame(False)
    ltext  = leg.get_texts()
    pl.setp(ltext, fontsize=14)


    ax2 = pl.twinx()
    ax2.yaxis.tick_right()
    #Plot the rain first
    pl.ylabel(r'$Rainfall \  (mm)$',fontsize=18,color=color_P)
    #pl.plot(X,ar_rain,'b-')
    X_fill=np.concatenate((X[0:1],X,X[-1:]))
    rain_fill=np.concatenate((sp.zeros(1),ar_rain,sp.zeros(1)))
    p=pl.fill(X_fill,rain_fill,facecolor=color_P,edgecolor=color_P,alpha=transparency_P)
    pl.ylim(max(ar_rain)*2,min(ar_rain))
    pl.xlim(X[0],X[-1])

    pl.savefig(image)


#SUBROUTINE
#``````````````````````````````````````````   
def read_delta_flow(file_name):
    import scipy.io as io
    f=file(file_name,'r')
    tab=io.read_array(f)
    f.close()
    ar_date=tab[:,0:-1]
    ar_Q=tab[:,-1]
    return ar_date,ar_Q



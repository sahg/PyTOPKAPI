#General module importation
import numpy as np
import scipy as sp
import pylab as pl
import datetime as dt
from matplotlib.dates import drange,date2num,num2date,YearLocator, MonthLocator,DayLocator, HourLocator,DateFormatter
import tables as h5
from ConfigParser import SafeConfigParser
config = SafeConfigParser()

#Personnal module importation
from TOPKAPI import utils as ut
#reload(ut)


def run(ini_file='plot_Qsim_Qobs_Rain.ini'):
    
    config.read(ini_file)
    print 'Read the file ',ini_file

    file_Qsim=config.get('files','file_Qsim')
    file_Qobs=config.get('files','file_Qobs')
    file_rain=config.get('files','file_rain')
    image_out=config.get('files','image_out')
    
    group_name=config.get('groups','group_name')

    Qobs=config.getboolean('flags','Qobs')
    Pobs=config.getboolean('flags','Pobs')
    nash=config.getboolean('flags','nash')
            
    tab_leg=['Model']
    #tab_col=['k','0.5']
    tab_col=['k','r']
    tab_style=['-','-']
    #tab_width=['1','2']
    tab_width=['1','1']
    #color_P='0.5'
    color_P='b'
    transparency_P=0.5#(0 for invisible)

    #create path_out if it does'nt exist
    ut.check_file_exist(image_out)

    #Read the obs
    #Qobs
    ar_date,ar_Qobs=read_delta_flow(file_Qobs)
    date1=dt.datetime(ar_date[0,0],ar_date[0,1],ar_date[0,2],ar_date[0,3],ar_date[0,4])
    date1prim=dt.datetime(ar_date[1,0],ar_date[1,1],ar_date[1,2],ar_date[1,3],ar_date[1,4])
    date2=dt.datetime(ar_date[-1,0],ar_date[-1,1],ar_date[-1,2],ar_date[-1,3],ar_date[-1,4])
    delta=date1prim-date1
    X= drange(date1, date2+delta, delta)
    #Rain
    if Pobs:
        h5file_in=h5.openFile(file_rain,mode='r')
        group='/'+group_name+'/'
        node = h5file_in.getNode(group+'rainfall')
        ndar_rain=node.read()
        h5file_in.close()
        #Compute the mean catchment rainfall
        ar_rain=np.average(ndar_rain,axis=1)

    #Read the simulated data Q
    file_h5=file_Qsim
    ndar_Qc_out=ut.read_one_array_hdf(file_h5,'/Channel/','Qc_out')
    ar_Qsim=ndar_Qc_out[:,0]


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
        

    pl.plot(X,ar_Qsim,color=tab_col[0],linestyle=tab_style[0],linewidth=tab_width[0])

    if nash:
        nash_value=ut.Nash(ar_Qsim,ar_Qobs)
        pl.plot(X[0:1],ar_Qsim[0:1],'w:')
        tab_leg.append(('Eff= '+str(nash_value)[0:5]))
            
    pl.xlim(X[0],X[-1])
    ytitle=r'$Q \  (m^3/s)$'
    pl.ylabel(ytitle,fontsize=18)
    pl.title(group_name)
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

    pl.savefig(image_out)


#SUBROUTINE
#``````````````````````````````````````````   
def read_delta_flow(file_name):
    import scipy.io as io

    tab=io.read_array(file_name)
    ar_date=np.array(tab[:,0:-1], np.int32) # date components are integer numbers
    ar_Q=tab[:,-1]
    return ar_date,ar_Q



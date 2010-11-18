#General module importation
import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
from matplotlib.dates import drange,num2date,YearLocator, MonthLocator,DayLocator,DateFormatter
import tables as h5
from ConfigParser import SafeConfigParser
config = SafeConfigParser()

#Personnal module importation
import pytopkapi.utils as ut

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
    ar_Qsim=ndar_Qc_out[1:,0]


    ##Graph
    plt.clf()
    #Set up the date axes
    lim_day0=3;lim_day1=61;lim_day2=3*365
    ax=plt.subplot(111)

    d=num2date(X[-1])-num2date(X[0])
    time_length=d.days
    years    = YearLocator()
    months   = MonthLocator()
    days     = DayLocator()
    ax.xaxis.set_major_locator(months)
    Fmt = DateFormatter("%m/%y")
    ax.xaxis.set_major_formatter(Fmt)
    labels = ax.get_xticklabels()
    #plt.setp(labels,'rotation', 90)
    if time_length>lim_day1 and time_length<lim_day2:
        ax.xaxis.set_major_locator(months)
        Fmt = DateFormatter("%m/%y")
        ax.xaxis.set_major_formatter(Fmt)
        labels = ax.get_xticklabels()
        #plt.setp(labels,'rotation', 90)
    elif time_length>=lim_day2:
        ax.xaxis.set_major_locator(years)
        ax.xaxis.set_minor_locator(months)
        Fmt = DateFormatter("%m/%y")
        ax.xaxis.set_major_formatter(Fmt)
        labels = ax.get_xticklabels()
        #plt.setp(labels,'rotation', 90)
    elif time_length>lim_day0 and time_length<lim_day1:
        ax.xaxis.set_major_locator(months)
        ax.xaxis.set_minor_locator(days)
        Fmt = DateFormatter("%d/%m/%y")
        ax.xaxis.set_major_formatter(Fmt)
        #plt.setp(labels,'rotation', 90)
    else:
        ax.xaxis.set_major_locator(hours)
        ax.xaxis.set_minor_locator(hours)
        Fmt = DateFormatter("%H")
        ax.xaxis.set_major_formatter(Fmt)

    #plot the graph
    if Qobs:
        plt.plot(X,ar_Qobs,color=tab_col[-1],linestyle=tab_style[-1],linewidth=tab_width[-1])
        tab_leg.append(('Observation'))
        tab_leg=tab_leg[::-1]
        

    plt.plot(X,ar_Qsim,color=tab_col[0],linestyle=tab_style[0],linewidth=tab_width[0])

    if nash:
        nash_value=ut.Nash(ar_Qsim,ar_Qobs)
        plt.plot(X[0:1],ar_Qsim[0:1],'w:')
        tab_leg.append(('Eff= '+str(nash_value)[0:5]))
            
    plt.xlim(X[0],X[-1])
    ytitle=r'$Q \  (m^3/s)$'
    plt.ylabel(ytitle,fontsize=18)
    plt.title(group_name)
    plt.legend(tab_leg,'center')
    leg = plt.gca().get_legend()
    leg.draw_frame(False)
    ltext  = leg.get_texts()
    plt.setp(ltext, fontsize=14)


    ax2 = plt.twinx()
    ax2.yaxis.tick_right()
    #Plot the rain first
    plt.ylabel(r'$Rainfall \  (mm)$',fontsize=18,color=color_P)
    #plt.plot(X,ar_rain,'b-')
    X_fill=np.concatenate((X[0:1],X,X[-1:]))
    rain_fill=np.concatenate((np.zeros(1),ar_rain,np.zeros(1)))
    p=plt.fill(X_fill,rain_fill,facecolor=color_P,edgecolor=color_P,alpha=transparency_P)
    plt.ylim(max(ar_rain)*2,min(ar_rain))
    plt.xlim(X[0],X[-1])

    plt.savefig(image_out)


#SUBROUTINE
#``````````````````````````````````````````   
def read_delta_flow(file_name):
    import scipy.io as io

    tab=np.loadtxt(file_name)
    ar_date=np.array(tab[:,0:-1], np.int32) # date components are integer numbers
    ar_Q=tab[:,-1]
    return ar_date,ar_Q



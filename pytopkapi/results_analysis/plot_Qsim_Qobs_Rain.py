import datetime as dt
from ConfigParser import SafeConfigParser

import numpy as np
import tables as h5
import matplotlib.pyplot as plt
from matplotlib.dates import date2num

import pytopkapi.utils as ut

def run(ini_file='plot_Qsim_Qobs_Rain.ini'):
    config = SafeConfigParser()
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

    tab_col=['k','r']
    tab_style=['-','-']
    tab_width=['1','1']
    color_P='b'
    transparency_P=0.5#(0 for invisible)

    #create path_out if it does'nt exist
    ut.check_file_exist(image_out)

    #Read the obs
    #Qobs
    ar_date, ar_Qobs = read_observed_flow(file_Qobs)

    delta = date2num(ar_date[1]) - date2num(ar_date[0])

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
    ndar_Qc_out=ut.read_one_array_hdf(file_h5,'Channel','Qc_out')
    ar_Qsim=ndar_Qc_out[1:,0]

    ##Graph
    fig, ax = plt.subplots()

    lines = []
    tab_leg = []
    if Qobs:
        lines += ax.plot(ar_date, ar_Qobs,
                         color=tab_col[-1],
                         linestyle=tab_style[-1], linewidth=tab_width[-1])
        tab_leg.append(('Observation'))
        tab_leg = tab_leg[::-1]

    lines += ax.plot(ar_date, ar_Qsim,
                     color=tab_col[0],
                     linestyle=tab_style[0], linewidth=tab_width[0])
    tab_leg.append('Model')

    if nash:
        nash_value = ut.Nash(ar_Qsim,ar_Qobs)
        lines += ax.plot(ar_date[0:1], ar_Qsim[0:1], 'w:')
        tab_leg.append(('Eff = '+str(nash_value)[0:5]))

    ax.set_xlim(ar_date[0], ar_date[-1])
    ytitle=r'$Q \  (m^3/s)$'
    ax.set_ylabel(ytitle, fontsize=18)
    ax.set_title(group_name)

    ax2 = ax.twinx()

    ax2.set_ylabel(r'$Rainfall \ (mm)$', fontsize=18, color=color_P)
    ax2.bar(ar_date, ar_rain, width=delta,
            facecolor='blue', edgecolor='blue', alpha=transparency_P)
    ax2.set_ylim(max(ar_rain)*2, min(ar_rain))

    ax2.legend(lines, tab_leg, loc='upper right', fancybox=True)
    leg = ax2.get_legend()
    leg.get_frame().set_alpha(0.75)

    # rotate and align the tick labels so they look better,
    # unfortunately autofmt_xdate doesn't work with twinx due to a bug
    # in matplotlib <= 1.0.0 so we do it manually
    ## fig.autofmt_xdate()

    bottom=0.2
    rotation=30
    ha='right'

    for ax in fig.get_axes():
        if hasattr(ax, 'is_last_row') and ax.is_last_row():
            for label in ax.get_xticklabels():
                label.set_ha(ha)
                label.set_rotation(rotation)
        else:
            for label in ax.get_xticklabels():
                label.set_visible(False)
            ax.set_xlabel('')

    fig.subplots_adjust(bottom=bottom)

    fig.savefig(image_out)
    plt.show()

def read_observed_flow(file_name):
    """Read the observed flow from a data file.

    """
    date = np.loadtxt(file_name, dtype=np.int, usecols=(0, 1, 2, 3, 4))
    dates = [dt.datetime(yr, mon, dy, hr, mn) for yr, mon, dy, hr, mn in date]

    Q = np.loadtxt(file_name, usecols=(5,))

    return dates, Q

"""
*** OBJECTIVE

1. creating the TOPKAPI parameter file from GIS binary files (binary grid 
   format).
   -main routine "creat_param_file"

*** COMMENT
In the present program, to write the parameter file, generally the grid are 
tranformed in an array by rearranging the cell order from the West to East 
and North to South, as in the following example:

The grid format is like:

GIS bingrid file 
-9999  -9999  -9999  -9999  -9999  -9999  -9999
-9999      0      1  -9999  -9999  -9999  -9999
-9999      2      3     4   -9999  -9999  -9999
    5      6      7     8       9  -9999  -9999 
-9999     10     11    12      13     14  -9999
-9999  -9999     15    16      17  -9999  -9999
-9999  -9999  -9999  -9999     18  -9999  -9999
-9999  -9999  -9999 -9999   -9999  -9999  -9999

The corresponding array extracted and ordered from West to East, North 
to South is:
0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18

"""
import sys
from shutil import copyfile

#Python modules
import numpy as np
import pylab as pl
from numpy import ma
import scipy.io as io
import tables as h5
from ConfigParser import SafeConfigParser
config = SafeConfigParser()

#External modules from TOPKAPI
from TOPKAPI import utils as ut
from TOPKAPI import pretreatment as pm
from TOPKAPI import arcfltgrid

def bingrid_to_label(file_grid, file_label='label.dat', write_file=False):
    """
    * Objective
      Replace values of the grid by the label of the cells.
      The label are assigned from 0 to nb_cell, from West to East, North to South.
    * Input
      file_bin_grid: A binary grid file (whatever it is).
    * Output
      tab: a 1D array with nb_cell components.
    """
    tab, headers = arcfltgrid.read(file_grid)
    
    nrows=np.shape(tab)[0]
    ncols=np.shape(tab)[1]
    tab=np.reshape(tab,ncols*nrows)
    ind=np.where(tab>-99)
    tab_label=np.arange(len(tab[ind]))
    for i in tab_label:
        tab[ind[0][i]]=i
    tab=np.reshape(tab,(nrows,ncols))
    tab=tab.astype('int32')

    if write_file:
        f = file(file_label, 'w')
        io.write_array(f, tab)
        f.close()

    return tab

def from_asciigrid_to_label(file_ascii_grid,file_label='label.dat',write_file=False):
    """
    * Objective
      Replace values of the grid by the label of the cells.
      The label are assigned from 0 to nb_cell, from West to East, North to South.
    * Input
      file_bin_grid: A binary grid file (whatever it is).
    * Output
      tab: a 1D array with nb_cell components.
    """
    tab=np.loadtxt(file_ascii_grid)
    nrows=np.shape(tab)[0]
    ncols=np.shape(tab)[1]
    tab=np.reshape(tab,ncols*nrows)
    ind=np.where(tab>-99)
    tab_label=np.arange(len(tab[ind]))
    for i in tab_label:
        tab[ind[0][i]]=i
    tab=np.reshape(tab,(nrows,ncols))
    tab=tab.astype('int32')

    if write_file:
        f = file(file_label, 'w')
        io.write_array(f, tab)
        f.close()

    return tab

def from_flowacc_to_stream(file_flowacc_grid,file_stream_grid,threshold_cell):
    tab_flowacc=np.loadtxt(file_flowacc_grid)
    nrows=np.shape(tab_flowacc)[0]
    ncols=np.shape(tab_flowacc)[1]  

    ar_flowacc=np.reshape(tab_flowacc,ncols*nrows)
    ar_stream=np.array(ar_flowacc)
    
    ar_stream[ar_flowacc>=threshold_cell]=1
    ar_stream[np.where((ar_flowacc<threshold_cell) & (ar_flowacc>-1))]=0

    total_cell=len(ar_flowacc[ar_flowacc>-1.])
    stream_cell=len(ar_stream[ar_stream==1])

    print 'total_cell=',total_cell,'stream_cell=',stream_cell,'Drainage density=',float(stream_cell)/float(total_cell)
    
    tab_stream=np.reshape(ar_stream,(nrows,ncols))
    
    f = file(file_stream_grid, 'w')
    io.write_array(f, tab_stream)
    f.close()


def compute_slope_8D(file_flowdir, file_DEM, 
                     file_slope_degree, file_slope, Xcell=1000.0):
    """Compute the channel slopes from 8D flow direction and a DEM.

    Calculate the slope from the centre of each cell in the catchment DEM
    to it's downstream neighbour. The calculated slopes and the tangents of 
    the slopes are written to seperate text files. These slopes may then be 
    used for the channel slopes required by the TOPKAPI model.
        
    Parameters
    ----------
    file_flowdir : string
        Precipitation input to the cell during the current time-step (mm).
    file_DEM : string
        The length of the current time-step in seconds.
    file_slope_degree : string
        The lateral dimension of the grid-cell (in m).
    file_slope : string
        The total contribution from each cell to it's downstream neighbour as a 
        result of subsurface and overland fluxes calculated during the previous 
        timestep (m^3/s).
    Xcell : `float`
        List of integer indices into `ar_Q_to_next_cell`. The indices point to 
        the cells upstream of the current cell.
        
    """
    tab_DEM=np.loadtxt(file_DEM)
    tab_dir=np.loadtxt(file_flowdir)
    tab_label=from_asciigrid_to_label(file_DEM)
    tab_slope_degree=np.array(tab_label)
    tab_slope=np.array(tab_label)

    nrows=np.shape(tab_DEM)[0]
    ncols=np.shape(tab_DEM)[1]

    for i in range(nrows):
        for j in range(ncols):

            label=tab_label[i,j]
            direction=tab_dir[i,j]
            if label>0:
                if direction==64:
                    dist=Xcell
                    x=i-1
                    y=j
                if direction==1:
                    dist=Xcell
                    x=i
                    y=j+1
                if direction==4:
                    dist=Xcell
                    x=i+1
                    y=j
                if direction==16:
                    dist=Xcell
                    x=i
                    y=j-1
                if direction==32:
                    dist=Xcell*(2**0.5)
                    x=i-1
                    y=j-1
                if direction==128:
                    dist=Xcell*(2**0.5)
                    x=i-1
                    y=j+1
                if direction==2:
                    dist=Xcell*(2**0.5)
                    x=i+1
                    y=j+1
                if direction==8:
                    dist=Xcell*(2**0.5)
                    x=i+1
                    y=j-1
                if tab_label[x,y]>=0:
                    if tab_DEM[x,y]<tab_DEM[i,j]:
                        tab_slope_degree[i,j]=np.arctan((tab_DEM[i,j]-tab_DEM[x,y])/dist)*180./np.pi
                        tab_slope[i,j]=(tab_DEM[i,j]-tab_DEM[x,y])/dist
                    if tab_DEM[x,y]==tab_DEM[i,j]:
                        tab_slope_degree[i,j]=0.
                        tab_slope[i,j]=0.
                    if tab_DEM[x,y]>tab_DEM[i,j]:
                        print 'Problem negative slope cell',tab_label[i,j],direction,tab_label[x,y],tab_DEM[i,j],tab_DEM[x,y]
                        tab_slope_degree[i,j]=np.arctan((tab_DEM[i,j]-tab_DEM[x,y])/dist)*180./np.pi
                        tab_slope[i,j]=(tab_DEM[i,j]-tab_DEM[x,y])/dist
                else:
                    print 'Problem cell external to the catchment...'
                    print tab_label[i,j],direction,tab_label[x,y],tab_DEM[i,j],tab_DEM[x,y]

    f = file(file_slope_degree, 'w')
    io.write_array(f, tab_slope_degree)
    f.close()
    
    f = file(file_slope, 'w')
    io.write_array(f, tab_slope)
    f.close()

def create_channel_slope_file(file_flowdir, file_DEM, file_slope_degree):
    """Compute the channel slopes from 8D flow direction and a DEM.

    Calculate the slope from the centre of each cell in the catchment DEM
    to it's downstream neighbour. The calculated slopes are written to 
    an ArcGIS Float grid file. This file is required by the TOPKAPI model
    for the generation of a parameter file.
        
    Parameters
    ----------
    file_flowdir : string
        Name of an ArcGIS binary Float file containing the flow 
        direction raster.
    file_DEM : string
        Name of an ArcGIS binary Float file containing the DEM raster.
    file_slope_degree : string
        Name of an ArcGIS binary Float file for the output raster.
        
    """
    dem, headers = arcfltgrid.read(file_DEM)
    nrows=np.shape(dem)[0]
    ncols=np.shape(dem)[1]
    Xcell = headers[4]

    flowdir, headers = arcfltgrid.read(file_flowdir)
    
    tab_label = bingrid_to_label(file_DEM)
    tab_slope_degree = np.array(tab_label, np.float32)

    for i in range(nrows):
        for j in range(ncols):

            label = tab_label[i,j]
            direction = flowdir[i,j]
            
            if label > 0:
                if direction==64:
                    dist=Xcell
                    x=i-1
                    y=j
                if direction==1:
                    dist=Xcell
                    x=i
                    y=j+1
                if direction==4:
                    dist=Xcell
                    x=i+1
                    y=j
                if direction==16:
                    dist=Xcell
                    x=i
                    y=j-1
                if direction==32:
                    dist=Xcell*(2**0.5)
                    x=i-1
                    y=j-1
                if direction==128:
                    dist=Xcell*(2**0.5)
                    x=i-1
                    y=j+1
                if direction==2:
                    dist=Xcell*(2**0.5)
                    x=i+1
                    y=j+1
                if direction==8:
                    dist=Xcell*(2**0.5)
                    x=i+1
                    y=j-1
                
                if tab_label[x,y] >= 0:
                    if dem[x,y] < dem[i,j]:
                        tab_slope_degree[i,j] = np.arctan((dem[i,j]-dem[x,y])
                                                          /dist) * 180./np.pi
                    if dem[x,y] == dem[i,j]:
                        tab_slope_degree[i,j] = 0.0
                    if dem[x,y] > dem[i,j]:
                        print 'Negative slope cell', tab_label[i,j], \
                            direction, tab_label[x,y], dem[i,j], dem[x,y]
                        tab_slope_degree[i,j] = np.arctan((dem[i,j]-dem[x,y])
                                                          /dist) * 180./np.pi
                else:
                    tab_slope_degree[i,j] = -9999
                    print 'Downslope cell external to the catchment...'
                    print tab_label[i,j], direction, \
                          tab_label[x,y], dem[i,j], dem[x,y]

    # write slope file
    tab_slope_degree.tofile(file_slope_degree)

    if file_DEM[-4:] == '.flt':
        copyfile(file_DEM[:-4] + '.hdr', file_slope_degree[:-4] + '.hdr')
    else:
        copyfile(file_DEM + '.hdr', file_slope_degree[:-4] + '.hdr')
    

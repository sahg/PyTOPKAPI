"""
*** OBJECTIVE

1. creating the TOPKAPI parameter file from GIS binary files (binary grid format).
   -main routine "creat_param_file"

*** COMMENT
In the present program, to write the parameter file, generally the grid are tranformed in an array
by rearranging the cell order from the West to East and North to South, as in the following example:

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

The corresponding array extracted and ordered from West to East, North to South is:
0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18

"""
import sys

#Python modules
import numpy as np
import pylab as pl
from numpy import ma
import scipy.io as io
import tables as h5
from ConfigParser import SafeConfigParser
config = SafeConfigParser()

#External modules from TOPKAPI
#Utilities
from TOPKAPI import utils as ut
#pretreatment: used for subroutines to read the column type files.
from TOPKAPI import pretreatment as pm


####### READING THE BINARY FILES ###########

def read_bin_data(fname):
    """Read data from a binary file into an array.
    
    Read the data from a binary file created by wgrib into
    a Numpy array. This assumes that wgrib has output the data
    in float format and seems to work for the UM outputs that
    are being sent to us by SAWS. The first and last values
    read from the file are trimmed since they are superflous.

    """
    
    f = open(fname, "rb")
    raw = f.read()
    f.close()
    data = np.fromstring(raw, 'f')
    if sys.byteorder == 'big':
            data = data.byteswap()

    return data

def read_arc_bin(bingrid_name):
    """Read the data field from a grib file into an array.

    This function uses wgrib to write and then read the GRIB1
    data to an intermediate binary format, which is later
    deleted. It is currently very specific to the GRIB1 data
    files produced by the UM.

    """

    li_headers=read_headers_arc_bin(bingrid_name)

    rows = li_headers[1] # fixed for UM grid
    cols = li_headers[0] # fixed for UM grid
    
    bin_name = bingrid_name + '.flt'

    a = read_bin_data(bin_name)

    a = a.reshape(rows, cols)

    return a

def read_headers_arc_bin(bingrid_name):
    """Read the ascii headers of the binary grid file:
    ncols         62
    nrows         121
    xllcorner     -288595.47161281
    yllcorner     -3158065.5722693
    cellsize      1000
    NODATA_value  -9999
    byteorder     LSBFIRST
    """

    hdr_name = bingrid_name + '.hdr'
    f=open(hdr_name,'r')
    tab_read=f.readlines()
    f.close()

    li_headers=[]
    i=-1
    for line in tab_read:
        i=i+1
        donnees=line.split()
        if i<6:
            li_headers.append(float(donnees[1]))
        else:
            li_headers.append(donnees[1])
            
    return li_headers

def arc_bin_plot(bin_name, fig_name, title='GRIB Plot'):
    """Create a plot of the data in a GRIB1 file."""
    
    a = read_arc_bin(bin_name)
    
    a_mask = ma.masked_where(a < 0, a)
    pl.imshow(a_mask, interpolation='nearest')
    pl.colorbar()
    pl.title(title)
    pl.savefig(fig_name)
    pl.close()


##def from_DEM_to_cell_down_and_slope(file_DEM_grid,file_label):
##    """
##    """
##    tab_DEM=io.read_array(file_DEM_grid)
##    nrows=np.shape(tab_DEM)[0]
##    ncols=np.shape(tab_DEM)[1]
##    tab_label=io.read_array(file_label)
##
##    tab_dir=np.array(tab_DEM)
##    tab_slope=np.array(tab_DEM)
##
##    
##    #Image of tab_DEM
##    tab_DEM=np.array(tab_DEM)
##    
##    ar=np.reshape(tab_DEM,ncols*nrows)
##    ar_DEM=ar[ar>-99.]
##
##    ar_cell_down=np.zeros(len(ar_DEM))
##    ar_tan_beta=np.zeros(len(ar_DEM))
##
##    ncell=-1
##    nb_sink=0
##    for i in range(nrows):
##        for j in range(ncols):
##            label=tab_label[i,j]
##            li_select_cells=[]
##            li_select_DEM=[]
##            li_select_dir=[]
##            if label>0:
##                ncell=ncell+1
##                #North
##                if i>0:
##                    x=i-1
##                    y=j
##                    if tab_label[x,y]>-1:
##                        li_select_cells.append(tab_label[x,y])
##                        li_select_DEM.append(tab_DEM[x,y])
##                        li_select_dir.append(64)
##                #East
##                if j<max(range(ncols)):
##                    x=i
##                    y=j+1
##                    if tab_label[x,y]>-1:
##                        li_select_cells.append(tab_label[x,y])
##                        li_select_DEM.append(tab_DEM[x,y])
##                        li_select_dir.append(1)
##                #South
##                if i<max(range(nrows)):
##                    x=i+1
##                    y=j
##                    if tab_label[x,y]>-1:
##                        li_select_cells.append(tab_label[x,y])
##                        li_select_DEM.append(tab_DEM[x,y])
##                        li_select_dir.append(4)
##                #West
##                if j>0:
##                    x=i
##                    y=j-1
##                    if tab_label[x,y]>-1:
##                        li_select_cells.append(tab_label[x,y])
##                        li_select_DEM.append(tab_DEM[x,y])
##                        li_select_dir.append(16)
##                    
##                
##                ar_select_cells=np.array(li_select_cells)
##                ar_select_DEM=np.array(li_select_DEM)
##                ar_select_dir=np.array(li_select_dir)
##
##                min_DEM=min(ar_select_DEM)
##                ind_min_DEM=np.where(ar_select_DEM==min_DEM)[0][0]
##                
##
##                if min_DEM>=tab_DEM[i,j]:
##                    nb_sink=nb_sink+1
##                    if min_DEM>tab_DEM[i,j]:
##                        print 'SINK FOR CELL',tab_label[i,j]
##                    else:
##                        print 'SINK EQUAL FOR CELL',tab_label[i,j]
##                    if len(ar_select_cells)==1. and tab_label[i,j]>0.:
##                        print 'Problem with cell number',tab_label[i,j],i,j
##                        stop
##                    tab_DEM[i,j]=np.average(ar_select_DEM)
##                    
##                tab_dir[i,j]=ar_select_dir[ind_min_DEM]
##                tab_slope[i,j]=np.arctan((tab_DEM[i,j]-min_DEM)/1000.)*180./np.pi
##                
##                ar_cell_down[ncell]=ar_select_cells[ind_min_DEM]
##                ar_tan_beta[ncell]=tab_slope[i,j]
## 
##                tab_dir[0,0]=1
##                tab_slope[0,0]=tab_slope[1,0]
##                
##    return tab_DEM,tab_slope,tab_dir,ar_cell_down,ar_tan_beta,nb_sink


def from_flowacc_to_stream(file_flowacc_grid,file_stream_grid,threshold_cell):
    tab_flowacc=io.read_array(file_flowacc_grid)
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


def compute_slope_8D(file_flowdir,file_DEM,file_slope_degree,file_slope,Xcell=1000.):
    tab_DEM=io.read_array(file_DEM)
    tab_dir=io.read_array(file_flowdir)
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
    tab=io.read_array(file_ascii_grid)
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

if __name__ == '__main__':
##    nb_sink=1
##    it=0
##    while nb_sink>0:
##        file_DEM_grid='C:/Theo/liebenbergsvlei/topkapi_model/parameters/GIS_ascii/DEM_fill4D_'+str(it)+'.TXT'
##        file_label='C:/Theo/liebenbergsvlei/topkapi_model/parameters/GIS_ascii/label.dat'
##        print file_DEM_grid
##        it=it+1
##        file_DEM_fill='C:/Theo/liebenbergsvlei/topkapi_model/parameters/GIS_ascii/DEM_fill4D_'+str(it)+'.TXT'
##        file_slope='C:/Theo/liebenbergsvlei/topkapi_model/parameters/GIS_ascii/slope_degrees_'+str(it)+'.TXT'
##        file_dir='C:/Theo/liebenbergsvlei/topkapi_model/parameters/GIS_ascii/flowdir_4D_'+str(it)+'.TXT'
##        tab_DEM,tab_slope,tab_dir,ar_cell_down,ar_tan_beta,nb_sink=from_DEM_to_cell_down_and_slope(file_DEM_grid,file_label)
##        #'help io.write_array' for more info
##        f = file(file_DEM_fill, 'w')
##        io.write_array(f, tab_DEM)
##        f.close()
##        f = file(file_slope, 'w')
##        io.write_array(f, tab_slope)
##        f.close()
##        f = file(file_dir, 'w')
##        io.write_array(f, tab_dir)
##        f.close()
    
##    file_flowacc_grid='C:/Theo/liebenbergsvlei/topkapi_model/parameters/GIS_ascii/flowacc8d_cor.txt'
##    file_stream_grid='C:/Theo/liebenbergsvlei/topkapi_model/parameters/GIS_ascii/stream8d_cor.txt'
##    threshold_cell=25.
##    from_flowacc_to_stream(file_flowacc_grid,file_stream_grid,threshold_cell)

##    file_flowdir='C:/Theo/liebenbergsvlei/topkapi_model_Aug07/parameters/GIS_ascii_files/flowdir_8d_cor.txt'
##    file_DEM='C:/Theo/liebenbergsvlei/topkapi_model_Aug07/parameters/GIS_ascii_files/DEM_fill8D_cor.txt'
##    file_label='C:/Theo/liebenbergsvlei/topkapi_model_Aug07/parameters/GIS_ascii_files/label.dat'
##    file_slope_degree='C:/Theo/liebenbergsvlei/topkapi_model_Aug07/parameters/GIS_ascii_files/slope_fill8D_degree.txt'
##    file_slope='C:/Theo/liebenbergsvlei/topkapi_model_Aug07/parameters/GIS_ascii_files/slope_fill8D_tan.txt'
##    compute_slope_8D(file_flowdir,file_DEM,file_label,file_slope_degree,file_slope,Xcell=1000.)

    file_flowdir='E:/post_doc/liebenbergsvlei/topkapi_model_Oct07/parameters/GIS_ascii_files/flowdir_4d_moh.txt'
    file_DEM='E:/post_doc/liebenbergsvlei/topkapi_model_Oct07/parameters/GIS_ascii_files/DEM_original.txt'
    file_label='E:/post_doc/liebenbergsvlei/topkapi_model_Oct07/parameters/GIS_ascii_files/label.dat'
    file_slope_degree='E:/post_doc/liebenbergsvlei/topkapi_model_Oct07/parameters/GIS_ascii_files/slope_channel_4D_degree_real.txt'
    file_slope='E:/post_doc/liebenbergsvlei/topkapi_model_Oct07/parameters/GIS_ascii_files/slope_channel_4D_tan_real.txt'
    compute_slope_8D(file_flowdir,file_DEM,file_label,file_slope_degree,file_slope,Xcell=1000.)

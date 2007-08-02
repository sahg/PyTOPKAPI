"""
*** OBJECTIVE

1. creating the TOPKAPI parameter file from GIS binary files (binary grid format).
   -part "FROM GIS TO PARAMETER FILE"
   -main routine "creat_param_file"
2. modifying the created parameter files in order to
   a. extract the parameter for a subcatchment
   b. changing the value of one of the parameters (for instance for calibration or sensitivity analyses)
   - part "MODIFYING THE PARAMETER FILES"
   - main routines

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

#Internal modules
import numpy as np
import pylab as pl
from numpy import ma
import scipy.io as io
import tables as h5

#External modules
#Utilities
import utils as ut
#pretreatment_module: used for subroutines to read the parameter files
import pretreatment_module as pm


##################################
### FROM GIS TO PARAMETER FILE ###
##################################

def creat_param_file(file_bin_streamnet,file_bin_beta,file_bin_flowdir,\
                     file_bin_GLCC,file_bin_SIRI,\
                     file_bin_WRC90,file_bin_strahler,\
                     file_table_GLCC_manning,file_table_SIRI_soil,\
                     file_table_WRC90_soil,file_table_strahler_manning,\
                     Vs_t0=80.,Vc_t0=0.,kc=1\
                     file_out):
    """
    * objective:
      Create a parameter file from the catchment data (GIS maps), in association with the TABLES povided in the litterature.
      
    * Input:
      +Binary grid files (from GIS)
      - file_bin_streamnet: grid of the streamnet (code: 1 for channel cell, 0 for hillslope cell)
      - file_bin_beta: grid of cell slopes in degress.
      - file_bin_flowdir: grid of flow directions
      - file_bin_GLCC: grid file containing the GLCC land use/cover codes
      - file_bin_SIRI: grid of soil properties (code given by SIRI)
      - file_bin_WRC90: binary grid file containing the WRC90 soil property codes (Here only three codes are considered (cf. comments): 3 for Loamy Sand, 2 for Sandy Loam, 1 for Clay)
      - file_bin_strahler: grid of strahler order of the channel cells.
      +Tables from litterature for transformaing grids into TOPKAPI physical parameters
      - file_table_GLCC_manning: an ASCII file containing a table of correspondance between the GLCC land use/cover codes
                                and the values of manning's coef proposed in different references
                                (Chow et al., 1998; Maidment,1993 - give a range of values and the MUSIC report)
      - file_table_SIRI_soil: an ASCII file containing a table of correspondance between the SIRI codes
                              and the values of a selection of the soil parameters proposed by SIRI
                              (Land-type, depth A, depth B, WP A, WP B, FC A, FC B, Por A, Por B)
      - file_table_WRC90_soil: an ASCII file containing a table of correspondance between the WRC90 codes
                               and the values of Ks (permeability) and theta_r (residual soil moisture) from Maidment(1993)
      - file_table_strahler_manning: an ASCII file containing a table of correspondance between the strahler order
                                     and the values of manning strickler proposed by Liu and Todini (2002)
      +Constant parameters
      - Vs_t0: Constant value for the initial soil saturation of each cell.
      - Vc_t0: Constant value for the initial channel water 
      - kc: Crop coefficient
      
    * Output:
      - file_out: parameter file (ASCII column format) containing:
       label X Y lambda L Ks Theta_r Theta_s n_o n_c cell_down pVs_t0 pVc_t0 kc


      
    * Comment:
     1. !!!!VERY IMPORTANT!!!!
     The routine refers to several subroutines listed below that must be carefully read before running the programm. This programm is helpfull for
     creating the parameter file but it is not automated.
     Some valuable informations are required inside these subroutines especially for data that are assigned using Tables.
     The Tables are indeed simple ASCII files that were created for the special case of the Liebenbergsvlei catchment.
     The application to another catchment might require the modification of the Tables and thus the subroutines.
     One must refer to the headers of each subroutine for the detailed information.

     2. Note that here, the initial soil moisture value, as well as the initial channel saturation are constant.
        Assigning a spatially variable initial soil moisture can be done through the routine "from_param_to_new_param_catchVsi" in this file.

    """

    #~~~~Number of parameter (used to fix the final array dimension)~~~#
    nb_param=15.

    #~~~~~Paremeters directly read ~~~~~~#
    ar_lambda=from_grid_to_param(file_bin_streamnet)
    #beta is given in degres -->tan(beta) is computed
    ar_tan_beta=np.tan(np.pi/180.*from_grid_to_param(file_bin_beta))

    #~~~~~Parameters extracted from GIS and estimated from TABLES~~~~~#
    ### !!!USER MUST CHECK THESE 4 SUBROUTINES BEFORE RUNNING THE CODE!!! ###
    ar_n_o=from_GLCC_to_manning(file_bin_GLCC,file_table_GLCC_manning)
    ar_L,ar_Theta_s=from_SIRI_to_soil_properties(file_bin_SIRI,file_table_SIRI_soil)
    ar_Theta_r,ar_Ks=from_WRC90_to_soil_properties(file_bin_WRC90,file_table_WRC90_soil)
    ar_n_c=from_Strahler_to_channel_manning(file_bin_strahler,file_table_strahler_manning,np.array(ar_lambda))

    #~~~~~Parameters computed~~~~~#
    ar_label=np.arange(len(ar_lambda))
    ar_coorX,ar_coorY=from_bingrid_to_coordinate(file_bin_streamnet)
    ar_cell_down=from_flowdir_to_celldown_4D(file_bin_flowdir)
    
    #~~~~~ Initial values and crop factor (CONSTANT) ~~~~~~#
    #Creation of a vector zero
    ar_tab=ar_label*0.
    #Constant value for the initial soil moisture
    ar_pVs_t0=ar_tab+Vs_t0
    #Constant value for the initial channel water volume
    ar_pVc_t0=ar_tab+Vc_t0
    #Constant value for the crop coefficient
    ar_kc=ar_tab+kc

    #~~~~~~Write parameter file~~~~~~#
    tab_param=np.zeros((len(ar_tab),nb_param))
    tab_param[:,0]=ar_label
    tab_param[:,1]=ar_coorX
    tab_param[:,2]=ar_coorY
    tab_param[:,3]=ar_lambda
    tab_param[:,4]=ar_tan_beta
    tab_param[:,5]=ar_L
    tab_param[:,6]=ar_Ks
    tab_param[:,7]=ar_Theta_r
    tab_param[:,8]=ar_Theta_s
    tab_param[:,9]=ar_n_o
    tab_param[:,10]=ar_n_c
    tab_param[:,11]=ar_cell_down
    tab_param[:,12]=ar_pVs_t0
    tab_param[:,13]=ar_pVc_t0
    tab_param[:,14]=ar_kc

    #'help io.write_array' for more info
    f = file(file_out, 'w')
    io.write_array(f, tab_param)
    f.close()

##############################################    
###  SUBROUTINE USED IN "creat_param_file" ###
##############################################

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


############# OTHER SUBROUTINES ###################
def from_grid_to_param(file_bin_grid):
    """
    * Objective
    Directly extract the values at the cell locations from the grid
    * Input:
    GIS binary grid
    * Output
    A 1D array with nb_cell values. Cells are ordered from West to East, North to South.
    """
    tab=read_arc_bin(file_bin_grid)
    nrows=np.shape(tab)[0]
    ncols=np.shape(tab)[1]
    tab=np.reshape(tab,ncols*nrows)
    ind=np.where(tab>-99)
    ar_param=tab[ind]
    
    return ar_param

def from_GLCC_to_manning(file_bin_GLCC,file_table_GLCC_manning):
    """
    * Objective:
      Estimation of the overland manning's coefficient for each catchment cell from the GLCC land cover/use map
    * Input
      - file_bin_GLCC is the binary grid file containing the GLCC land use/cover codes
      - file_table_GLCC_manning is an ASCII file containing a table of correspondance between the GLCC land use/cover codes
      and the values of manning's coef proposed in different references (Chow et al., 1998; Maidment,1993 - give a range of values and the MUSIC report)
    * Ouput
      This routine returns a 1D array (ar_n_o) containing the values of the Manning's coefficient for each cell. Cells are ordered from West to East, North to South.
    * Comment:
      This routine has to be used carefully since the code is dependant to:
       1. the format of the Table (file_table_GLCC_manning)
       2. the choice of the user to use a given reference 
    """
    #Read the binary grid file of GLCC land use type
    tab=read_arc_bin(file_bin_GLCC)
    nrows=np.shape(tab)[0]
    ncols=np.shape(tab)[1]
    tab=np.reshape(tab,ncols*nrows)
    ind=np.where(tab>-99)
    ar_GLCC=tab[ind]

    #Read the Table file within a header line
    tab=pm.read_column_input(file_table_GLCC_manning,5)
    ar_code=tab[:,0]
    ar_n_chow=tab[:,1]
    ar_n_maidment_inf=tab[:,2]
    ar_n_maidment_sup=tab[:,3]
    ar_n_music=tab[:,4]

    ar_n_o=np.array(ar_GLCC)
    for i in ar_code:
        ind=np.where(ar_GLCC==i)
        #!!!! TO BE CHANGED ACCORDING TO USER CHOICE !!!#
        if i>1:
            ar_n_o[ind]=ar_n_chow[np.where(ar_code==i)][0]
        if i==1:
            ar_n_o[ind]=ar_n_music[np.where(ar_code==i)][0]
        
    return ar_n_o

def from_SIRI_to_soil_properties(file_bin_SIRI,file_table_SIRI_soil):
    """
    * Objective:
      Extraction of the parameters L (soil depth) and theta_s (porosity or humidity at saturation) for each catchment cell from the SIRI map
    * Input
      - file_bin_SIRI is the binary grid file containing the SIRI soil property codes
      - file_table_SIRI_soil is an ASCII file containing a table of correspondance between the SIRI codes
      and the values of a selection of the soil parameters proposed by SIRI (Land-type, depth A, depth B, WP A, WP B, FC A, FC B, Por A, Por B)
    * Ouput
      This routine returns two 1D array (ar_L, ar_theta_s) containing respectively the values of L and theta_s for each cell. Cells are ordered from West to East, North to South.
    * Comment:
      This routine has to be used carefully since the code is dependant to:
       1. the format of the Table (file_table_SIRI_soil)
       2. the choice of the user to use only one soil layer A, or combining the two soil layer A+B characteristics provided by SIRI.
       Here L=LA+LB and Theta_S=average(theta_sA,theta_sB)
    """
    #Read the binary grid file of GLCC land use type
    tab=read_arc_bin(file_bin_SIRI)
    nrows=np.shape(tab)[0]
    ncols=np.shape(tab)[1]
    tab=np.reshape(tab,ncols*nrows)
    ind=np.where(tab>-99)
    ar_SIRI=tab[ind]
    
    #Read the Table file within a header line
    tab=pm.read_column_input(file_table_SIRI_soil,9)
    ar_code=tab[:,0]
    ar_depthA=tab[:,1]
    ar_depthB=tab[:,2]
    ar_porA=tab[:,7]
    ar_porB=tab[:,8]

    ar_L=np.array(ar_SIRI)
    ar_theta_s=np.array(ar_SIRI)
    for i in ar_code:
        ind=np.where(ar_L==i)
        #!!!! TO BE CHANGED ACCORDING TO USER CHOICE !!!#
        ar_L[ind]=ar_depthA[np.where(ar_code==i)][0]+ar_depthB[np.where(ar_code==i)][0]
        ar_theta_s[ind]=0.5*(ar_porA[np.where(ar_code==i)][0]+ar_porB[np.where(ar_code==i)][0])
    
    return ar_L, ar_theta_s

def from_WRC90_to_soil_properties(file_bin_WRC90,file_table_WRC90_soil):
    """
    * Objective:
      Extraction of the parameters L (soil depth) and theta_s (porosity or humidity at saturation) for each catchment cell from the SIRI map
    * Input
      - file_bin_WRC90 is the binary grid file containing the WRC90 soil property codes (Here only three 3 for Loamy Sand, 2 for Sandy Loam, 1 for Clay)
      - file_table_WRC90_soil is an ASCII file containing a table of correspondance between the WRC90 codes and the values of Ks (permeability) and theta_r (residual soil moisture)
    * Ouput
      This routine returns two 1D array (ar_theta_r, ar_theta_s) containing respectively the values of Ks and theta_r for each cell. Cells are ordered from West to East, North to South.
    """
    #Read the binary grid file of GLCC land use type
    tab=read_arc_bin(file_bin_WRC90)
    nrows=np.shape(tab)[0]
    ncols=np.shape(tab)[1]
    tab=np.reshape(tab,ncols*nrows)
    ind=np.where(tab>-99)
    ar_WRC90=tab[ind]

    #Read the Table file within a header line
    tab=pm.read_column_input(file_table_WRC90_soil,9)
    ar_code=tab[:,0]
    ar_theta_r_moy=tab[:,3]
    ar_theta_r_ect=tab[:,4]
    ar_conduct=tab[:,8]

    ar_theta_r=np.array(ar_WRC90)
    ar_Ks=np.array(ar_WRC90)
    for i in ar_code:
        ind=np.where(ar_WRC90==i)
        #!!!! TO BE CHANGED FOR PARAMETER ADJUSTMENT !!!#
        ar_theta_r[ind]=ar_theta_r_moy[np.where(ar_code==i)][0]
        ar_Ks[ind]=ar_conduct[np.where(ar_code==i)][0]
         
    return ar_theta_r,ar_Ks

def from_Strahler_to_channel_manning(file_bin_strahler,file_table_strahler_manning,ar_lambda):
    """
    * objective:
      Extraction of the channel mannings for each channel cell within the catchment from the strahler order map
    * Input
      - file_bin_strahler is the binary grid file containing the strahler order of the channel cells
      - file_table_strahler_manning is an ASCII file containing a table of correspondance between the strahler order
      and the values of manning strickler proposed by Liu and Todini (2002)
      - ar_lambda is an array (dimension equal to the number of cell) with value 1 for channel cells, 0 otherwise. Cells being ordered from West to East, North to South.
    * Ouput
      This routine returns a 1D array (ar_n_c), ar_theta_s) containing the values of the manning coefficient for each channel cell. Cells are ordered from West to East, North to South.
    """

    #Read the binary grid file of GLCC land use type
    tab=read_arc_bin(file_bin_strahler)
    nrows=np.shape(tab)[0]
    ncols=np.shape(tab)[1]
    tab=np.reshape(tab,ncols*nrows)
    ind=np.where(tab>-99)
    ar_strahler=tab[ind]
    ar_lambda[ar_lambda==1]=ar_strahler
    
    #Read the Table file within a header line
    tab=pm.read_column_input(file_table_strahler_manning,2)
    ar_code=tab[:,0]
    ar_manning=tab[:,1]

    ar_n_c=np.array(ar_lambda)
    for i in ar_code:
        ind=np.where(ar_n_c==i)
        #!!!! TO BE CHANGED FOR PARAMETER ADJUSTMENT !!!#
        ar_n_c[ind]=ar_manning[np.where(ar_code==i)][0]

    return ar_n_c


def from_flowdir_to_celldown_4D(file_bin_flowdir):
    """
    * objective:
      Associate to each cell the label of the celldown (outcell) according to the 4D flow directions:
    * Input
      - file_bin_flowdir: 4D flow direction binary file (from GIS) with codes: 1 right, 16 left, 64 up, 4 down
    * Output
      - ar_cell_down: an array (dimension equal to the number of cell) containing the label of the outcells.
    * Commment
     Some errors can occur in the direction files that should be detected by the routine. A manual correction is required if it does happen.
    """

    tab_flowdir=read_arc_bin(file_bin_flowdir)
    tab_label=from_bingrid_to_label(file_bin_flowdir)
    nrows=np.shape(tab_label)[0]
    ncols=np.shape(tab_label)[1]
    ar_cell_down=np.zeros(np.shape(np.where(tab_label>=0))[1])-99
    n=0
    num=0
    outlet=0
    for i in range(nrows):
        for j in range(ncols):
            flow=tab_flowdir[i,j]
            OK=0
            if tab_label[i,j] >= 0:
                num=num+1
                if flow==1:
                    x=i
                    y=j+1
                elif flow==16:
                    x=i
                    y=j-1
                elif flow==64:
                    x=i-1
                    y=j
                elif flow==4:
                    x=i+1
                    y=j
                else:
                    print 'ERROR'
                if tab_label[x,y]>=0:
                    ar_cell_down[tab_label[i,j]]=tab_label[x,y]
                else:
                    if outlet==0:
                        ar_cell_down[tab_label[i,j]]=tab_label[x,y]
                        outlet=1
                    else:
                        n=n+1
                        #Print the data to detect where the errors are (manual correction to be done...)
                        print n,num, flow,i,j,tab_label[i,j],x,y,tab_label[x,y]
    return ar_cell_down


def from_bingrid_to_label(file_bin_grid,write_file=False):
    """
    * Objective
      Replace values of the grid by the label of the cells.
      The label are assigned from 0 to nb_cell, from West to East, North to South.
    * Input
      file_bin_grid: A binary grid file (whatever it is).
    * Output
      tab: a 1D array with nb_cell components.
    """
    tab=read_arc_bin(file_bin_grid)
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
        f = file('c:/Theo/liebenbergsvlei/topkapi_model/parameters/label.dat', 'w')
        io.write_array(f, tab)
        f.close()

    return tab

def from_bingrid_to_coordinate(file_bin_grid):
    """
    * Objective
    Compute the coordinates of the cells relatively to
    the low-left corner from any binary grid file
    * Input
      file_bin_grid: A binary grid file (whatever it is).
    * Output
      ar_coorX,ar_coorY: two 1D array with nb_cell components.
    """
    li_headers=read_headers_arc_bin(file_bin_grid)
    ncols = li_headers[0] # fixed for UM grid
    nrows = li_headers[1] # fixed for UM grid
    xmin = li_headers[2]
    ymin = li_headers[3]
    cellsize= li_headers[4]
    
    ar_line_coorX=xmin+np.arange(ncols)*cellsize
    ar_line_coorY=ymin+np.arange(nrows)*cellsize
    ar_line_coorY=ar_line_coorY[::-1]
    
    mat_coorX=np.zeros((nrows,ncols))
    mat_coorY=np.zeros((nrows,ncols))
    for i in range(nrows):
        mat_coorX[i,:]=ar_line_coorX
    for j in range(ncols):
        mat_coorY[:,j]=ar_line_coorY

    tab=read_arc_bin(file_bin_grid)
    ind=np.where(tab>-99.)

    ar_coorX=mat_coorX[ind]
    ar_coorY=mat_coorY[ind]
    
    return ar_coorX,ar_coorY

def matrix_plot(matrix, fig_name, title='GRID Plot'):
    """Create a plot of the data in a GRIB1 file."""
    
    a=matrix
    
    a_mask = ma.masked_where(a < 0, a)
    pl.imshow(a_mask, interpolation='nearest')
    pl.colorbar()
    pl.title(title)
    pl.savefig(fig_name)
    pl.close()


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
    import pretreatment_module as pm
    import result_analysis as ra
    import os
    import os.path
    nb_param=15.
    
    #Reading of parameter file
    print 'Reading parameter file'
    ar_cell_label,ar_coorx,ar_coory,ar_lambda,ar_tan_beta,ar_L,ar_Ks,\
    ar_theta_r,ar_theta_s,ar_n_o,ar_n_c,\
    ar_cell_down,ar_pVs_t0,ar_pVc_t0,ar_kc\
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
    tab_param[:,4]=ar_tan_beta[subcatch_label]
    tab_param[:,5]=ar_L[subcatch_label]
    tab_param[:,6]=ar_Ks[subcatch_label]
    tab_param[:,7]=ar_theta_r[subcatch_label]
    tab_param[:,8]=ar_theta_s[subcatch_label]
    tab_param[:,9]=ar_n_o[subcatch_label]
    tab_param[:,10]=ar_n_c[subcatch_label]
    for i in range(len(subcatch_label)):
        if i==0.:
            tab_param[i,11]=-9999.
        else:
            tab_param[i,11]=new_label[np.where(ar_cell_label[subcatch_label]==ar_cell_down[subcatch_label][i])]
    tab_param[:,12]=ar_pVs_t0[subcatch_label]
    tab_param[:,13]=ar_pVc_t0[subcatch_label]
    tab_param[:,14]=ar_kc[subcatch_label]

    #~~~~~~Write parameter file~~~~~~#
    #'help io.write_array' for more info
    f = file(file_out, 'w')
    io.write_array(f, tab_param)
    f.close()

    image_out=os.path.split(file_out)[0]+'/'+nom_image+'.png'
    ar_image=ar_cell_label*0.;ar_image[subcatch_label]=1.;ar_image[ar_lambda==1.]=10.;ar_image[cell_outlet]=5.
    ra.field_map(ar_image,ar_coorx,ar_coory,X,image_out,'Subcatchment')

def from_param_to_new_param(file_in,file_out,fac_L=1.,fac_Ks=1.,fac_n_o=1.,fac_n_c=1.,new_pVs_t0=50.,new_pVc_t0=0.,channel_lesotho=False):
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
    import pretreatment_module as pm
    import result_analysis as ra
    import os
    import os.path
    
    nb_param=15.
    
    #Reading of parameter file
    print 'Reading parameter file'
    ar_cell_label,ar_coorx,ar_coory,ar_lambda,ar_tan_beta,ar_L,ar_Ks,\
    ar_theta_r,ar_theta_s,ar_n_o,ar_n_c,\
    ar_cell_down,ar_pVs_t0,ar_pVc_t0,ar_kc\
    =pm.read_cell_parameters(file_in)

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
    if new_pVs_t0!=ar_pVs_t0[0]:
        print 'Change pVs_t0'
        ar_pVs_t0=ar_pVs_t0*0.+new_pVs_t0
        print ar_pVs_t0[0]
    if new_pVc_t0!=ar_pVc_t0[0]:
        print 'Change pVc_t0'
        ar_pVc_t0=ar_pVc_t0*0.+new_pVc_t0

    if channel_lesotho:
        ar_lambda=channel_cell_tunnel(ar_cell_label,ar_coorx,ar_coory,ar_lambda,ar_cell_down)
        
    #~~~~~~Write parameter file~~~~~~#
    tab_param=np.zeros((len(ar_cell_label),nb_param))
    tab_param[:,0]=ar_cell_label
    tab_param[:,1]=ar_coorx
    tab_param[:,2]=ar_coory
    tab_param[:,3]=ar_lambda
    tab_param[:,4]=ar_tan_beta
    tab_param[:,5]=ar_L
    tab_param[:,6]=ar_Ks
    tab_param[:,7]=ar_theta_r
    tab_param[:,8]=ar_theta_s
    tab_param[:,9]=ar_n_o
    tab_param[:,10]=ar_n_c
    tab_param[:,11]=ar_cell_down
    tab_param[:,12]=ar_pVs_t0
    tab_param[:,13]=ar_pVc_t0
    tab_param[:,14]=ar_kc

    f = file(file_out, 'w')
    io.write_array(f, tab_param)
    f.close()

def from_param_to_new_param_catchVsi(file_in,file_in_global,file_out,file_h5,file_catchment_saturation,fac_L=1.,fac_Ks=1.,fac_n_o=1.,fac_n_c=1.,new_pVs_t0=-99.,new_pVc_t0=-99.,channel_lesotho=False):
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
    import pretreatment_module as pm
    import result_analysis as ra
    import os
    import os.path
    nb_param=15.

    
    #--Read and compute the parameters to have the values of parameter and ar_Vsm
    #~~~~Read Global parameters file
    print 'Reading parameter file'
    X,Dt,alpha_s,alpha_o,alpha_c,A_thres,W_min,W_max\
      =pm.read_global_parameters(file_in_global)
    #~~~~Read Cell parameters file
    ar_cell_label,ar_coorx,ar_coory,ar_lambda,ar_tan_beta,ar_L,ar_Ks,\
    ar_theta_r,ar_theta_s,ar_n_o,ar_n_c,\
    ar_cell_down,ar_pVs_t0,ar_pVc_t0,ar_kc\
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
    ar_Vsm, ar_b_s, ar_Vs_t0, ar_b_o, ar_W, ar_b_c, ar_Vc_t0\
      =pm.compute_cell_param(X,Dt,alpha_s,alpha_o,alpha_c,nb_cell,\
                              A_thres,W_max,W_min,\
                              ar_lambda,ar_tan_beta,ar_L,\
                              ar_Ks,ar_theta_r,ar_theta_s,ar_n_o,ar_n_c,\
                              ar_A_drained,ar_pVs_t0,ar_pVc_t0)

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
    tab_param[:,4]=ar_tan_beta
    tab_param[:,5]=ar_L
    tab_param[:,6]=ar_Ks
    tab_param[:,7]=ar_theta_r
    tab_param[:,8]=ar_theta_s
    tab_param[:,9]=ar_n_o
    tab_param[:,10]=ar_n_c
    tab_param[:,11]=ar_cell_down
    tab_param[:,12]=ar_pVs_t0
    tab_param[:,13]=ar_pVc_t0
    tab_param[:,14]=ar_kc

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

##    path_main='C:/Theo/liebenbergsvlei/topkapi_model/parameters/'
##    path_GIS=path_main+'GIS_bin_files/'
##    path_table=path_main+'tables_lieben/'
##    creat_param_file(path_GIS+'streamnet',path_GIS+'slope',path_GIS+'flowdir',\
##                     path_GIS+'land_use_GLCC', path_GIS+'soil_type_SIRI',\
##                     path_GIS+'soil_text_WRC90',path_GIS+'strahler_order',\
##                     path_table+'manning_from_GLCC_lieben.dat',path_table+'soil_properties_from_SIRI_lieben.dat',\
##                     path_table+'soil_physical_properties_from_WRC90_lieben.dat',path_table+'Strahler_manning.dat',\
##                     path_main+'cell_param_20Feb07_Vsi80%.dat')

##    path_main='C:/Theo/liebenbergsvlei/topkapi_model/parameters/'
##    file_in=path_main+'cell_param_20Feb07_Vsi10%.dat'
##    file_out=path_main+'subC8H020/cell_param_sub_Vsi10%.dat'
##    Xoutlet=-259210.38
##    Youtlet=-3065828.87
##    nom_image='subC8H020'
##    from_param_to_subcatch_param(file_in, file_out,Xoutlet,Youtlet,nom_image)
    
##    TO CHANGE THE INITIAL SOIL MOISTURE (CONSTANT)
##    path_main='C:/Theo/liebenbergsvlei/topkapi_model/parameters/'
##    file_in=path_main+'cell_param_20Feb07_Vsi80%.dat'

##    file_out=path_main+'cell_param_Vsi100%.dat'
##    from_param_to_new_param(file_in,file_out,new_pVs_t0=100.)
##    
##    file_out=path_main+'cell_param_C8H020_Vsi10%.dat'
##    from_param_to_new_param(file_in,file_out,new_pVs_t0=10.)
##    file_out=path_main+'cell_param_C8H020_Vsi20%.dat'
##    from_param_to_new_param(file_in,file_out,new_pVs_t0=20.)
##    file_out=path_main+'cell_param_C8H020_Vsi50%.dat'
##    from_param_to_new_param(file_in,file_out,new_pVs_t0=50.)
##    file_out=path_main+'cell_param_C8H020_Vsi80%.dat'
##    from_param_to_new_param(file_in,file_out,new_pVs_t0=80.)
##
##    file_out=path_main+'cell_param_C8H020_Ks53_Vsi10%.dat'
##    from_param_to_new_param(file_in,file_out,fac_Ks=53.0,new_pVs_t0=10.)
##    file_out=path_main+'cell_param_C8H020_Ks53_Vsi20%.dat'
##    from_param_to_new_param(file_in,file_out,fac_Ks=53.0,new_pVs_t0=20.)
##    file_out=path_main+'cell_param_C8H020_Ks53_Vsi50%.dat'
##    from_param_to_new_param(file_in,file_out,fac_Ks=53.0,new_pVs_t0=50.)
##    file_out=path_main+'cell_param_C8H020_Ks53_Vsi80%.dat'
##    from_param_to_new_param(file_in,file_out,fac_Ks=53.0,new_pVs_t0=80.)
##    file_out=path_main+'cell_param_C8H020_Ks53_Vsi100%.dat'
##    from_param_to_new_param(file_in,file_out,fac_Ks=53.0,new_pVs_t0=100.)

    path_main='C:/Theo/liebenbergsvlei/topkapi_model/parameters/'
    file_in=path_main+'cell_param_Vsi100%.dat'
    file_in_global=path_main+'global_param_6h.dat'
    file_h5='C:/Theo/liebenbergsvlei/topkapi_model/field_initial_soil_moisture/Event1_6h_L1.1_Ks101.0_no1.0_nc1.0_100%_zerorain.h5'
    file_catchment_saturation='C:/Theo/liebenbergsvlei/topkapi_model/field_initial_soil_moisture/Catchment_soil_moisture_evolution.dat'

    path_main='C:/Theo/liebenbergsvlei/topkapi_model/parameters/'
    tab=[10,20,30,40,50,60,70,80,90,100]
    for i in tab:
        print ''
        print i
        file_out=path_main+'cell_catchVsi'+str(i)+'%.dat'
        from_param_to_new_param_catchVsi(file_in,file_in_global,file_out,file_h5,file_catchment_saturation,new_pVs_t0=i,channel_lesotho=True)

##    path_main='C:/Theo/liebenbergsvlei/topkapi_model/parameters/'
##    file_in=path_main+'cell_C8H020_catchVsi30%.dat'
##    #Reading of parameter file
##    print 'Reading parameter file'
##    ar_cell_label,ar_coorx,ar_coory,ar_lambda,ar_tan_beta,ar_L,ar_Ks,\
##    ar_theta_r,ar_theta_s,ar_n_o,ar_n_c,\
##    ar_cell_down,ar_pVs_t0,ar_pVc_t0,ar_kc\
##    =pm.read_cell_parameters(file_in)
##
##    Xoutlet=-255032.83;Youtlet=-3149857.34
##    cell=find_cell_coordinates(ar_cell_label,Xoutlet,Youtlet,ar_coorx,ar_coory,ar_lambda,channel=False)

    

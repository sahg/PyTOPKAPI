""" model_module.py

Functions required by the TOPKAPI model for the management
of input and output of cells:

NB:The subroutines solving the differential equation are not in this module
   (see ODE_module for more information)

"""

__author__ = "Theo Vischel"
__version__ = "$Revision: 1.0 $"
__date__ = "$Date: 09/10/2006 $"

import scipy as sp
import numpy as np

##        ROUTINES FOR SOIL STORE
#```````````````````````````````````````````
def initial_volume_soil(ar_pVs_t0,ar_Vsm):
    """ initialize_volume_soil
    Compute the intial content of water (volume of water) from the
    initial saturated values given in the parameter file.
    """
    ar_Vs_t0=ar_pVs_t0/100.*ar_Vsm
    return ar_Vs_t0

def input_soil(P,Dt,X,ar_Q_to_next_cell,ar_cell_up):
    """ input_soil
        Compute the total input of soil cell:
        - Rainfall input
        + Upstream from subsurface
        + Upstream from overland
        
        Return: a_s
    """
    #Transformation of P in mm to P_flux in m^3/s
    P_flux=(P*(10**-3)/Dt)*X**2
    #Case 1: cell without up)
    ind=ar_cell_up[ar_cell_up>-90.]
    ar_sel=ar_Q_to_next_cell[ind]
    if ar_sel[ar_sel<0].size!=0:
        print ''
        print 'STOP-ERROR: By computing -->Upstream from cell n. ',ind[ar_sel<0],' missing'
        print ''
        a_s='Error on upstream'
        return a_s
    else:
        a_s=P_flux+ar_sel.sum()
        return a_s


#```````````````````````````````````````````
def output_soil(Vs_t0, Vs_t1_prim, Vsm, a_s, b_s, alpha_s, Dt):
    Q_max=b_s*Vsm**alpha_s
    if Vs_t1_prim> Vsm:
        Qs_out=Q_max
        Vs_out=Vsm
    else:
        Qs_out=Qout_computing(Vs_t0,Vs_t1_prim, a_s, Dt)
        Vs_out=Vs_t1_prim
    if Qs_out<0:
        print 'a=',a_s,'Vs_t1_prim=',Vs_t1_prim,'Vs_t0=',Vs_t0,'Vsm=',Vsm
        print 'Qs=',Qs_out
    return Qs_out,Vs_out


#```````````````````````````````````````````
def output_soil_parak(Vs_t0, Vs_t1_prim, Vsm, b_s, alpha_s):
    if Vs_t1_prim> Vsm:
        Vs_out=Vsm
    else:
        Vs_out=Vs_t1_prim
    Qs_out=Qout_computing2(Vs_t0,Vs_out, b_s, alpha_s)
    
    return Qs_out

#```````````````````````````````````````````
def input_overland(Vs_t0,Vs_prim,Vsm,a_s,b_s,alpha_s,Dt):
    Q_prim=Qout_computing2(Vs_t0,Vs_prim,a_s,Dt)
    Q_max=b_s*Vsm**alpha_s
    a_o=max(0,Q_prim-Q_max)

    return a_o, Q_max

#```````````````````````````````````````````
def Qout_computing(V_t0,V_t1_prim,a,Dt):
    """ Qout_computing
        Compute the output flows Qout from the computed water volumes:

        Return: Qout
    """
    Qout=a-(V_t1_prim-V_t0)/Dt
    return Qout

#```````````````````````````````````````````
def Qout_computing2(V_t0,V_t1,b,alpha):
    """ Qout_computing
        Compute the output flows Qout from the computed water volumes:

        Return: Qout
    """
    Qout=b/2.*(V_t1**alpha+V_t0**alpha)
    
    return Qout

#```````````````````````````````````````````
def flow_partitioning(Lambda,Qs_out,Qo_out,W,X,Xc):
    """ flow partitioning
        Partition the outflows from soil and overland into:
        - next cell
        - channel if channel cell

        Return: Q_to_next_cell, Q_to_channel, Q_to_cannel_sub
    """
    if Lambda!=0:
        Q_to_next_cell=(1-Lambda*W*X/(X**2))*(Qs_out+Qo_out)    
        Q_to_channel=(Lambda*W*X/(X**2))*(Qs_out+Qo_out)
        Q_to_channel_sub=(Lambda*W*X/(X**2))*(Qs_out)
    else:
        Q_to_next_cell=(Qs_out+Qo_out)    
        Q_to_channel=0.
        Q_to_channel_sub=0.

    return Q_to_next_cell,Q_to_channel,Q_to_channel_sub

#```````````````````````````````````````````
def input_channel(ar_Qc_out,Q_to_channel,ar_cell_up):
    """ input_channel
        Compute the total input of channel inchannel cell:
        + Upstream from channel cells
        + Flows from soil partitionig
        
        Return: a_c
    """
    ind=ar_cell_up[ar_cell_up>-90.]
    if len(ind)>0:
        ar_Qc_cell_up=ar_Qc_out[ind]
        a_c=Q_to_channel+ar_Qc_cell_up.sum()
    else:
        ar_Qc_cell_up=np.array([0.])
        a_c=Q_to_channel
    
    return a_c, ar_Qc_cell_up.sum()

#```````````````````````````````````````````
def manning_depth(Q,W,X,n):
    """
    Compute with the manning equation (high width hypothesis) the water depth h corresponding to
     the flow Q for the channel characteristics: width W, manning n.
    The volume is then returned V=hWX
    """
    h=(n*Q/(W**(3./2.)))**(2./5.)

    return h

def initial_volume_channel(ar_Q,ar_W,X,ar_n_c):
    """
    The initial volume ar_Vc_t0 is computed for all the channel cells, from a known initial flow.
    """
    nb_cell=len(ar_W)
    ar_Vc_t0=np.zeros(nb_cell)
    for i in  range(nb_cell):
        if ar_W[i]>0.:
            h=manning_depth(ar_Q[i],ar_W[i],X,ar_n_c[i])
            ar_Vc_t0[i]=h*ar_W[i]*X

    return ar_Vc_t0

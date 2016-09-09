"""
Functions required by the TOPKAPI model for the management of input
and output of cells.

.. note::

  The subroutines solving the differential equation are not in this
  module (see :mod:`~TOPKAPI.ode` for more information)

"""

import numpy as np

##        ROUTINES FOR SOIL STORE
#```````````````````````````````````````````
def initial_volume_soil(ar_pVs_t0, ar_Vsm):
    """ initialize_volume_soil
    Compute the intial content of water (volume of water) from the
    initial saturated values given in the parameter file.
    """
    ar_Vs_t0=ar_pVs_t0/100.*ar_Vsm
    return ar_Vs_t0

def input_soil(P, Dt, X, ar_Q_to_next_cell, ar_cell_up):
    """Compute the total input to a soil store.

    Calculate the total rate of input to a single soil store. This
    comprises the sum of rainfall input, subsurface contribution from
    upstream cells and the overland contribution from upstream cells.

    Parameters
    ----------
    P : scalar
        Precipitation input to the cell during the current time-step
        (:math:`mm`)
    Dt : scalar
        The length of the current time-step in seconds
    X : scalar
        The lateral dimension of the grid-cell (:math:`m`)
    ar_Q_to_next_cell : (N,) Numpy array
        The total contribution from each cell to it's downstream
        neighbour as a result of subsurface and overland fluxes
        calculated during the previous timestep (:math:`m^3/s`)
    ar_cell_up : list of `int`
        List of integer indices into `ar_Q_to_next_cell`. The indices
        point to the cells upstream of the current cell

    Returns
    -------
    a_s : scalar
        The input flux to the soil store during the current time-step
        in :math:`m^3/s`

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

def output_soil(Vs_t0, Vs_t1_prim, Vsm, a_s, b_s, alpha_s, Dt):
    """Compute the outflow from and volume in a soil store.

    Calculate the outflow rate and final volume in a soil store after
    transferring volume greater than the saturated soil moisture
    content to the overland store of the same model cell.

    Parameters
    ----------
    Vs_t0 : scalar
        Volume of water in the soil store at the beginning of the
        current time- step (:math:`m^3`)
    Vs_t1_prim : scalar
        Volume of water in the soil store at the end of the current
        time-step as a result of combined input to and drainage from
        the non-linear reservoir (:math:`m^3`)
    Vsm : scalar
        Volume of soil moisture for the store under saturated
        conditions (:math:`m^3`)
    a_s : scalar
        The input flux to the soil store during the current time-step
        (:math:`m^3/s`)
    b_s : scalar
        The constant term of the non-linear differential equation
        :math:`dV_s/dt = a_s - b_sV_s^{\\alpha_s}`
    alpha_s : scalar
        The dimensionless pore-size distribution parameter for the
        soil store
    Dt : scalar
        The length of the current time-step (:math:`s`)

    Returns
    -------
    Qs_out : scalar
        Rate of flow out of the soil store during the time-step
        (:math:`m^3/s`)
    Vs_out : scalar
        Volume remaining in the soil store at the end of the time-step
        (:math:`m^3`)

    """
    if Vs_t1_prim > Vsm:
        Q_max = b_s*Vsm**alpha_s

        Qs_out=Q_max
        Vs_out=Vsm
    else:
        Qs_out = Qout_computing(Vs_t0, Vs_t1_prim, a_s, Dt)
        Vs_out = Vs_t1_prim

    if Qs_out < 0:
        print 'a=', a_s, 'Vs_t1_prim=', Vs_t1_prim, 'Vs_t0=', Vs_t0, 'Vsm=', Vsm
        print 'Qs=', Qs_out

    return Qs_out, Vs_out

def output_soil_parak(Vs_t0, Vs_t1_prim, Vsm, b_s, alpha_s):
    if Vs_t1_prim> Vsm:
        Vs_out=Vsm
    else:
        Vs_out=Vs_t1_prim
    Qs_out=Qout_computing2(Vs_t0,Vs_out, b_s, alpha_s)

    return Qs_out

def input_overland(Vs_t0,Vs_prim,Vsm,a_s,b_s,alpha_s,Dt):
    Q_prim=Qout_computing2(Vs_t0,Vs_prim,a_s,Dt)
    Q_max=b_s*Vsm**alpha_s
    a_o=max(0,Q_prim-Q_max)

    return a_o, Q_max

def Qout_computing(V_t0, V_t1_prim, a, Dt):
    """Compute the outflow `Qout` from a generic water store.

    Calculate the mean outflow rate during the current time-step from
    a generic water store. The outflow is calculated as the difference
    between the inflow and rate of change in storage during the
    time-step.

    Parameters
    ----------
    V_t0 : scalar
        Volume of water at the start of the current time-step
        (:math:`m^3`)
    V_t1_prim : scalar
        Volume of water at the end of the current time-step
        (:math:`m^3`)
    a : scalar
        The inflow rate to the water store during the time-step
        (:math:`m^3/s`)
    Dt : scalar
        The length of the current time-step (:math:`s`)

    Returns
    -------
    Qout : scalar
        The outflow rate from the water store during the time-step
        (:math:`m^3/s`)

    """
    b = (V_t1_prim - V_t0)/Dt

    if np.isclose(a, b):
        Qout = 0
    else:
        Qout = a - b

    return Qout

def Qout_computing2(V_t0,V_t1,b,alpha):
    """ Compute the output flows Qout from the computed water volumes:

    Returns
    -------
    Qout : scalar
        The outflow rate from the water store during the time-step
        (:math:`m^3/s`)

    """
    Qout=b/2.*(V_t1**alpha+V_t0**alpha)

    return Qout

def flow_partitioning(Lambda, Qs_out, Qo_out, W, X, Xc):
    """Partition the outflows from soil and overland stores.

    Calculates the correct partitioning of outflows from the soil and
    overland stores into contributions going to the downstream cell
    and/or the channel store of the current cell (if this exists).

    Parameters
    ----------
    Lambda : int
        Switch indicating whether the current cell contains a
        channel. A value of `1` indicates a channel cell, `0`
        indicates no channel
    Qs_out : scalar
        Outflow from the soil store during the current time-step
        (:math:`m^3/s`)
    Qo_out : scalar
        Outflow from the overland store during the current time-step
        (:math:`m^3/s`)
    W : scalar
        Width of the channel (:math:`m`)
    X : scalar
        The lateral dimension of the grid-cell (:math:`m`)
    Xc : scalar
        The length of the channel cell, this can be different from `X`
        if the channel runs along the cell diagonal (:math:`m`)

    Returns
    -------
    Q_to_next_cell : scalar
        Combined outflow from soil and overland to the downstream cell
        (:math:`m^3/s`)
    Q_to_channel : scalar
        Combined outflow from soil and overland to the channel store
        (:math:`m^3/s`)
    Q_to_channel_sub : scalar
        Outflow from soil store to the channel store (:math:`m^3/s`)

    """
    if Lambda != 0: #ToDo: check for invalid values of Lambda.
        Q_to_next_cell = (1-Lambda*W*X/(X**2))*(Qs_out+Qo_out)
        Q_to_channel = (Lambda*W*X/(X**2))*(Qs_out+Qo_out)
        Q_to_channel_sub = (Lambda*W*X/(X**2))*(Qs_out)
    else:
        Q_to_next_cell = (Qs_out+Qo_out)
        Q_to_channel = 0.
        Q_to_channel_sub = 0.

    return Q_to_next_cell, Q_to_channel, Q_to_channel_sub

def input_channel(ar_Qc_out, Q_to_channel, ar_cell_up):
    """Compute the total inflow to the channel of a channel cell.

    Calculate the total inflow to the channel as the sum of inflows
    from upstream cells and inflows from the soil and overland stores
    in the current cell.

    Parameters
    ----------
    ar_Qc_out : (N,) Numpy array
        Array of channel outflows from each cell in the catchment for
        the current time-step (:math:`m^3/s`)
    Q_to_channel : scalar
        Combined outflow from soil and overland to the channel store
        (:math:`m^3/s`)
    ar_cell_up : list of `int`
        List of integer indices into `ar_Qc_out`. The indices point to
        the cells upstream of the current cell.

    Returns
    -------
    a_c : scalar
        The total inflow to the channel store during the current
        time-step (:math:`m^3/s`)
    Qc_cell_up : scalar
        The contribution to the total inflow to the channel store from
        upstream cells during the current time-step (:math:`m^3/s`)

    """
    ind=ar_cell_up[ar_cell_up>-90.]
    if len(ind)>0:
        ar_Qc_cell_up=ar_Qc_out[ind]
        a_c=Q_to_channel+ar_Qc_cell_up.sum()
    else:
        ar_Qc_cell_up=np.array([0.])
        a_c=Q_to_channel

    Qc_cell_up = ar_Qc_cell_up.sum()

    return a_c, Qc_cell_up

def manning_depth(Q,W,X,n):
    """Compute Manning depth for flow `Q`.

    Compute with the manning equation (high width hypothesis) the
    water depth h corresponding to the flow Q for the channel
    characteristics: width W, manning n. The volume is then returned
    V=hWX

    """
    h=(n*Q/(W**(3./2.)))**(2./5.)

    return h

def initial_volume_channel(ar_Q,ar_W,X,ar_n_c):
    """Compute initial channel volume.

    The initial volume ar_Vc_t0 is computed for all the channel cells,
    from a known initial flow.

    """
    nb_cell=len(ar_W)
    ar_Vc_t0=np.zeros(nb_cell)
    for i in  range(nb_cell):
        if ar_W[i]>0.:
            h=manning_depth(ar_Q[i],ar_W[i],X,ar_n_c[i])
            ar_Vc_t0[i]=h*ar_W[i]*X

    return ar_Vc_t0

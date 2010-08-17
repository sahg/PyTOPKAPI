""" evap_module.py

Routines used to compute the evapotranspiration losses in TOPKAPI

"""

def evapot_soil_Liu_and_Todini_ETc(Vs0,Vsm,kc,ETr,X):
    """
    The evapotranspiration taken up from the soil:
      - at the reference crop rate ETc (potential evapotranspiration rate)
      - without overland runoff infiltration
    """
    #ETr is in mm
    #From Reference crop evapotranspiration
    #to crop evapotranspiration
    ETc=kc*ETr
    #to actual evapotranspiration
    ETa=ETc
    #From mm to m
    ETa=ETa*1e-3
    #Compute the new soil volume 
    Vs1=max(0.,Vs0-ETa*(X**2))
    #From m to mm
    ETa=ETa*1e3
    return ETa, Vs1

def evapot_soil_Liu_and_Todini(Vs0,Vsm,kc,ETr,X):
    """
    The evapotranspiration taken up from the soil:
      - at the rate ks*ETc ks being the soil saturation rate
      - without overland runoff infiltration
    """
    #ETr is in mm
    #From Reference crop evapotranspiration
    #to crop evapotranspiration
    ETc=kc*ETr
    #to actual evapotranspiration
    ks=Vs0/Vsm
    ETa=ks*ETc
    #From mm to m
    ETa=ETa*1e-3
    #Compute the new soil volume 
    Vs1=max(0.,Vs0-ETa*(X**2))
    #From m to mm
    ETa=ETa*1e3
    return ETa, Vs1

def evapot_soil_overland(Vo0, Vs0, Vsm, kc, ETr, X):
    """Compute the evapotranspiration from a model cell.
    
    The evapotranspiration loss is first taken from the overland store, if the
    storage in the overland store cannot satify the ET demand then the water is
    extracted from the soil store. In cases where the ET demand is greater than
    the available water, both the soil and overland store are totally drained.
      
    Parameters
    ----------
    Vo0 : scalar
        Volume in the overland store before ET removal (m^3).
    Vs0 : scalar
        Volume in the soil store before ET removal (m^3).
    Vsm : scalar
        Volume of soil moisture for the store under saturated conditions (m^3).
    kc : scalar
        Dimensionless crop co-efficient for the model cell.
    ETr : scalar
        Reference crop ET for the cell during the current time-step (mm).
    X : scalar
        The lateral dimension of the grid-cell (in m).
    
    Returns
    -------
    ETa : scalar
        The actual ET removed from the soil store, calculated as kc*ks*ETr with
        ks the relative saturation of the soil store before ET is removed (mm).
    Vs1 : scalar
        Volume in the soil store after ET removal (m^3).
    Vo1 : scalar
        Volume in the overland store after ET removal (m^3).
    
    """
    #ETr is in mm
    #From Reference crop evapotranspiration
    #to crop evapotranspiration
    ETc = kc*ETr
    #to actual evapotranspiration
    ks = Vs0/Vsm
    ETa = ks*ETc
    #From mm to m
    ETa = ETa*1e-3
    if Vo0 > 0:
        if Vo0-ETa*(X**2) >= 0:
            Vo1 = Vo0-ETa*(X**2)
            Vs1 = Vs0
        else:
            Vo1 = 0.
            Vs1 = max(0., Vs0-(ETa*(X**2)-Vo0))
    else:
        Vo1 = 0.
        Vs1 = max(0., Vs0-ETa*(X**2))
    #From m to mm
    ETa = ETa*1e3
    
    return ETa, Vs1, Vo1

def evapor_channel(Vc0, ETo, W, X):
    """Compute the evaporation from a channel store.
    
    Calculate the evaporation loss from a channel store. The function attempts
    to remove the demand from the open water potential evaporation from the 
    channel store, if there isn't sufficient water in the store then the total 
    volume in the channel is removed.
      
    Parameters
    ----------
    Vc0 : scalar
        Volume in the channel store before evaporation removal (m^3).
    ETo : scalar
        The potential evaporation from an open water surface (mm).
    W : scalar
        Width of the channel (m).
    X : scalar
        The length of the channel cell, this can be different from the cell
        dimension if the channel runs along the cell diagonal (m).
    
    Returns
    -------
    ET_channel : scalar
        The actual depth of water removed from the channel store (mm).
    Vc1 : scalar
        Volume in the channel store after evaporation removal (m^3).

    """
    ETo = ETo*1e-3
    if Vc0-ETo*W*X > 0:
        Vc1 = Vc0-ETo*W*X
        Vevap_c = ETo*W*X
    else:
        Vc1 = 0
        Vevap_c = Vc0
    ET_channel = Vevap_c*1e3/(W*X)
    
    return ET_channel, Vc1

def intercept_rain_ET(P,ETr,kc):
    """
    The evapotranspiration is taken from the precipitation:
      - at the reference crop rate ETc
    """
    ETc=kc*ETr
    if P-ETc>=0:
        Pn=P-ETc
        ETcn=0.
        ETa=ETc
    else:
        Pn=0.
        ETcn=(ETc-P)
        ETa=P
    ETrn=ETcn/kc
    return Pn,ETrn,ETa

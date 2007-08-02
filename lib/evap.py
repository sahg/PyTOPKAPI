""" evap_module.py

Routines used to compute the evapotranspiration losses in TOPKAPI

"""

__author__ = "Theo Vischel"
__version__ = "$Revision: 1.0 $"
__date__ = "$Date: 01/15/2007 $"

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

def evapot_soil_overland(Vo0,Vs0,Vsm,kc,ETr,X):
    """
    The evapotranspiration taken up from the soil:
      - at the rate ks*ETc ks being the soil saturation rate
      - with overland runoff infiltration
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
    if Vo0>0:
        if Vo0-ETa*(X**2)>=0:
            Vo1=Vo0-ETa*(X**2)
            Vs1=Vs0
        else:
            Vo1=0.
            Vs1=max(0.,Vs0-(ETa*(X**2)-Vo0))
    else:
        Vo1=0.
        Vs1=max(0.,Vs0-ETa*(X**2))
    #From m to mm
    ETa=ETa*1e3
    return ETa, Vs1,Vo1

def evapor_channel(Vc0,ETo,W,X):
    """
    The evaporation is taken from the channel:
      - at the rate ETo being the potential evaporation
    """
    ETo=ETo*1e-3
    if Vc0-ETo*W*X>0:
        Vc1=Vc0-ETo*W*X
        Vevap_c=ETo*W*X
    else:
        Vc1=0
        Vevap_c=Vc0
    ET_channel=Vevap_c*1e3/(W*X)
    return ET_channel,Vc1

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

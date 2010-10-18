"""Infiltration module.

"""
import numpy as np
from scipy.optimize import fsolve

def green_ampt_cum_infiltration(F, psi, dtheta, K, t):
    """The Green-Ampt cumulative infiltration equation.

    """
    tmp = psi*dtheta
    
    # np.log(x) computes ln(x)
    return F - tmp*np.log(1 + F/tmp) - K*t

if __name__ == '__main__':
    psi = 16.7
    dtheta = 0.34
    K = 0.65
    t = 1
    
    F = K*t # initial guess
    
    print fsolve(green_ampt_cum_infiltration,
                 F, args=(psi, dtheta, K, t), full_output=True)

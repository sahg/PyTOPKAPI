"""Infiltration module.

"""
import numpy as np
from scipy.optimize import fsolve

def _green_ampt_cum_eq(F, psi, dtheta, K, t):
    """The Green-Ampt cumulative infiltration equation

    """
    tmp = psi*dtheta

    # np.log(x) computes ln(x)
    return F - tmp*np.log(1 + F/tmp) - K*t

def green_ampt_cum_infiltration(psi, dtheta, K, t):
    """Compute the Green-Ampt cumulative infiltration

    Compute the potential cumulative infiltration up to time `t`,
    using Green-Ampt.

    Parameters
    ----------
    psi : array_like
        Soil suction head at wetting front.
    dtheta : array_like
        Ratio of initial effective saturation to effective porosity.
    K : array_like
        Saturated hydraulic conductivity.
    t : array_like
        Time since beginning of event

    Returns
    -------
    soln : array_like
        Cumulative infiltration up to time `t`.

    Raises
    ------
    ValueError - If no solution can be found.

    """

    F = K*t # initial guess

    soln, infodict, ierr, mesg = fsolve(_green_ampt_cum_eq, F,
                                        args=(psi, dtheta, K, t),
                                        full_output=True)

    if ierr == 1:
        return soln
    else:
        raise ValueError(mesg)

def test_basic_green_ampt():
    """Test the Green-Ampt cumulative infiltration solution"""

    psi = 16.7
    dtheta = 0.34
    K = 0.65
    t = 1

    result = green_ampt_cum_infiltration(psi, dtheta, K, t)

    assert np.allclose(result, [3.16641923])

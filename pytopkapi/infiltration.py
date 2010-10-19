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

def green_ampt_cum_infiltration(psi, eff_theta, eff_sat, K, t):
    """Compute the Green-Ampt cumulative infiltration

    Compute the potential cumulative infiltration up to time `t`,
    using Green-Ampt.

    Parameters
    ----------
    psi : scalar
        Soil suction head at wetting front.
    eff_theta : scalar
        Effective porosity.
    eff_sat : scalar
        Effective saturation.
    K : scalar
        Saturated hydraulic conductivity.
    t : scalar
        Time since beginning of event

    Returns
    -------
    soln : scalar
        Cumulative infiltration up to time `t`.

    Raises
    ------
    ValueError - If no solution can be found.

    """

    dtheta = (1 - eff_sat)*eff_theta

    F = K*t # initial guess

    soln, infodict, ierr, mesg = fsolve(_green_ampt_cum_eq, F,
                                        args=(psi, dtheta, K, t),
                                        full_output=True)

    if ierr == 1:
        return soln
    else:
        raise ValueError(mesg)

def green_ampt_infiltration_rate(psi, eff_theta, eff_sat, K, t):
    """Compute the Green-Ampt infiltration rate

    Compute the potential infiltration rate at time `t`, using
    Green-Ampt.

    Parameters
    ----------
    psi : scalar
        Soil suction head at wetting front.
    eff_theta : scalar
        Effective porosity.
    eff_sat : scalar
        Effective saturation.
    K : scalar
        Saturated hydraulic conductivity.
    t : scalar
        Time since beginning of event

    Returns
    -------
    ft : scalar
        Infiltration rate at time `t`.

    Raises
    ------
    ValueError - If no solution can be found.

    """

    F = green_ampt_cum_infiltration(psi, eff_theta, eff_sat, K, t)

    dtheta = (1 - eff_sat)*eff_theta

    return K*((psi*dtheta)/F + 1)

def test_basic_green_ampt():
    """Test the Green-Ampt infiltration solution

    This test implements Example 4.3.1 pg 116 from Chow et al. 1988,
    'Applied Hydrology'

    """

    psi = 16.7
    eff_theta = 0.486
    eff_sat = 0.3
    K = 0.65
    t = 1

    cum_infil = green_ampt_cum_infiltration(psi, eff_theta, eff_sat, K, t)

    assert np.allclose(cum_infil, [3.16721373])

    infil_rate = green_ampt_infiltration_rate(psi, eff_theta, eff_sat, K, t)

    assert np.allclose(infil_rate, [1.81596836])

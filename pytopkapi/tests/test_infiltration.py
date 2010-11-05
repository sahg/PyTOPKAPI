"""Tests for the infiltration module

"""
import numpy as np
from scipy.optimize import fsolve

from pytopkapi.infiltration import _green_ampt_cum_eq
from pytopkapi.infiltration import green_ampt_cum_infiltration

def test_basic_green_ampt():
    """Test the Green-Ampt infiltration solution

    This test solves Example 4.3.1 pg 116 from Chow et al. 1988,
    'Applied Hydrology'

    """
    psi = 16.7
    eff_theta = 0.486
    eff_sat = 0.3
    K = 0.65
    dt = 1

    dtheta = (1 - eff_sat)*eff_theta

    F = K*dt # initial guess
    soln, infodict, ierr, mesg = fsolve(_green_ampt_cum_eq, F,
                                        args=(0, psi,
                                              dtheta, K, dt),
                                        full_output=True)
    cum_infil = soln[0]

    assert np.allclose(cum_infil, [3.16721373])

def test_cum_green_ampt_variable_rainfall():
    """Test Green-Ampt with variable rainfall inputs

    This test solves Example 5.4.1 pg 144 from Chow et al., 1988,
    'Applied Hydrology'

    """
    psi = 11.01
    eff_theta = 0.412
    eff_sat = 0.4
    K = 1.09
    dt = 10/60.0
    rain = np.array([0.18,
                     0.21,
                     0.26,
                     0.32,
                     0.37,
                     0.43,
                     0.64,
                     1.14,
                     3.18,
                     1.65,
                     0.81,
                     0.52,
                     0.42,
                     0.36,
                     0.28,
                     0.24,
                     0.19,
                     0.17])

    expected = np.array([0,
                         0.17999999999999999,
                         0.39000000000000001,
                         0.65000000000000002,
                         0.96999999999999997,
                         1.3399999999999999,
                         1.7699999999999998,
                         2.2010980254004107,
                         2.5894006663182325,
                         2.9497199417670243,
                         3.2899523056814881,
                         3.6148926600284406,
                         3.9277082030903774,
                         4.2306187029814302,
                         4.5252501037735513,
                         4.8052501037735516,
                         5.0452501037735518,
                         5.2352501037735522,
                         5.4052501037735521])

    cum_rain = np.concatenate(([0], np.cumsum(rain)))

    F = [0]
    Ft = 0.0
    for r in rain/dt:
        Ft = green_ampt_cum_infiltration(r, psi,
                                           eff_theta, eff_sat, K, dt, Ft)
        F.append(Ft)

    assert(np.allclose(F, expected))

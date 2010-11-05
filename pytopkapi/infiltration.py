"""Infiltration module.

"""
import numpy as np
from scipy.optimize import fsolve

def _green_ampt_cum_eq(F_t1, F_t0, psi, dtheta, K, dt):
    """The Green-Ampt cumulative infiltration equation

    """
    tmp = psi*dtheta

    # np.log(x) computes ln(x)
    return F_t1 - F_t0 - tmp*np.log((F_t1 + tmp)/(F_t0 + tmp)) - K*dt

def _green_ampt_infiltration_rate(F, psi, eff_theta, eff_sat, K):
    """Compute the Green-Ampt infiltration rate

    Compute the infiltration rate using the Green-Ampt cumulative
    infiltration.

    Parameters
    ----------
    F : scalar
        The cumulative infiltration for the time-period.
    psi : scalar
        Soil suction head at wetting front.
    eff_theta : scalar
        Effective porosity.
    eff_sat : scalar
        Effective saturation.
    K : scalar
        Saturated hydraulic conductivity.

    Returns
    -------
    ft : scalar
        Infiltration rate for the given `F`.

    Raises
    ------
    ValueError - If `F` is zero or negative.

    """
    if F <= 0:
        raise ValueError('F must be greater than zero.')

    dtheta = (1 - eff_sat)*eff_theta

    return K*((psi*dtheta)/F + 1)

def green_ampt_cum_infiltration(rain, psi, eff_theta,
                                eff_sat, K, dt, F_t0=0.0):
    """Compute the Green-Ampt cumulative infiltration

    Compute the Green-Ampt cumulative infiltration depth for a given
    rainfall rate and time interval of length `dt`.

    Parameters
    ----------
    rain : scalar
        The constant rainfall rate during the current time interval of
        length `dt`
    psi : scalar
        Soil suction head at wetting front.
    eff_theta : scalar
        Effective porosity.
    eff_sat : scalar
        Effective saturation.
    K : scalar
        Saturated hydraulic conductivity.
    dt : scalar
        Length of the time-step
    F_t0 : scalar
        The cumulative infiltration depth at the beginning of the
        current interval. Default value is zero.

    Returns
    -------
    F_t1 : scalar
        Cumulative infiltration at the end of the current time-step.

    Raises
    ------
    ValueError - If no solution can be found.

    """
    if rain == 0:
        # shortcut: cum infiltration during interval must be zero
        F_t1 = F_t0
    else:
        # compute potential infiltration rate at the start of the
        # interval
        if F_t0 == 0:
            # infinite rate - just make f_t0 > rainfall rate
            f_t0 = rain + 100.0
        else:
            f_t0 = _green_ampt_infiltration_rate(F_t0, psi,
                                                 eff_theta, eff_sat, K)

        if f_t0 <= rain:
            # ponding occurs throughout interval
            dtheta = (1 - eff_sat)*eff_theta

            F = K*dt # initial guess
            soln, infodict, ierr, mesg = fsolve(_green_ampt_cum_eq, F,
                                                args=(F_t0, psi,
                                                      dtheta, K, dt),
                                                full_output=True)
            F_t1 = soln[0]

            if ierr != 1:
                raise ValueError(mesg)
        else:
            # check whether ponding occurs during interval
            Fprime = F_t0 + rain*dt
            fprime = _green_ampt_infiltration_rate(Fprime, psi,
                                                   eff_theta, eff_sat, K)

            if fprime > rain:
                # no ponding infiltrate at rain rate
                F_t1 = Fprime
            else:
                # ponding has occurred
                dtheta = (1 - eff_sat)*eff_theta
                Fp = (K*psi*dtheta)/(rain - K)

                dtprime = (Fp - F_t0)/rain

                F = K*dt # initial guess
                dtp = dt - dtprime
                soln, infodict, ierr, mesg = fsolve(_green_ampt_cum_eq, F,
                                                    args=(Fp, psi,
                                                          dtheta, K, dtp),
                                                    full_output=True)
                F_t1 = soln[0]

                if ierr != 1:
                    raise ValueError(mesg)

    return F_t1

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

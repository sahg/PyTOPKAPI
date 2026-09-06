"""Reference and open-water evaporation formulas for the ET forcing.

Pure functions, no IO: they take numbers (scalars or arrays) and return numbers,
so every one is checkable against the FAO-56 worked examples (Allen et al.,
1998).  The ET build stage calls these per station per day on the daily table
:mod:`met` produces, then disaggregates and interpolates.

Three outputs matter:

* :func:`et0_penman_monteith` -- FAO-56 reference-crop ET (grass), the ``ETr``
  the solver reads.  This is ``ET0`` in FAO notation.
* :func:`et0_hargreaves` -- the temperature-only fallback, for a station with no
  humidity/wind/radiation.  Keeps a first run possible on thin data.
* :func:`e0_open_water` -- Penman (1948) open-water evaporation, the ``ETo`` the
  channel loses at.  Its fallback, when radiation is absent, is ``ETo = ETr``
  (channel evaporation is a sub-percent term in the water balance, so this is
  cheap and defensible -- see the ET design note).

**All three return mm/day.**  The split into mm-per-``Dt`` depths, with a
clear-sky diurnal weight, is ``gauges.disaggregate``'s job downstream, not this
module's.  A value left in mm/day would over-evaporate the catchment by a factor
of ``86400/Dt`` -- the units trap the design note names.

Inputs and their units
-----------------------
======================  ======  =========================================
temperature             degC     ``tmax``, ``tmin``, ``tmean``, ``tdew``
relative humidity       %        ``rh`` (0-100)
wind speed              m/s      ``u`` (adjusted to 2 m by :func:`wind_speed_2m`)
solar radiation         W/m^2    daily-mean (:func:`rs_from_wm2` -> MJ/m^2/day)
elevation               m        station height (for air pressure -> gamma)
latitude                degrees  north positive, south negative
day of year             1-366    (:func:`day_of_year`)
======================  ======  =========================================

Solar is the feed's confirmed unit (W/m^2, daily-mean); :func:`rs_from_wm2`
converts it to the MJ/m^2/day the radiation balance needs.

References: Allen, R.G., Pereira, L.S., Raes, D., Smith, M. (1998), FAO-56;
Hargreaves & Samani (1985); Penman (1948).
"""

from __future__ import annotations

import numpy as np
import pandas as pd

__all__ = [
    "GSC", "SIGMA", "ALBEDO_GRASS", "ALBEDO_WATER", "WM2_TO_MJ_DAY",
    "svp", "svp_slope", "mean_svp",
    "atm_pressure", "psychrometric_constant",
    "ea_from_tdew", "ea_from_rhmean", "ea_from_rhminmax",
    "day_of_year",
    "extraterrestrial_radiation", "daylight_hours", "clear_sky_radiation",
    "rs_from_wm2", "net_shortwave", "net_longwave", "net_radiation",
    "wind_speed_2m",
    "penman_monteith_combine", "et0_penman_monteith",
    "et0_hargreaves", "e0_open_water",
    "reference_et0", "et0_eto_frame",
    # sub-daily / hourly (Allen 2006)
    "SIGMA_HOUR", "CN_DAILY", "CD_DAILY",
    "CN_HOURLY", "CD_HOURLY_DAY", "CD_HOURLY_NIGHT",
    "seasonal_correction", "extraterrestrial_radiation_hourly",
    "daily_cloudiness_factor", "rs_rate_from_wm2", "net_radiation_hourly",
    "et0_penman_monteith_hourly",
    # routing
    "et0_eto_on_clock",
]

#: Solar constant (MJ m^-2 min^-1), FAO-56 Eq. 28.
GSC = 0.0820
#: Stefan-Boltzmann (MJ K^-4 m^-2 day^-1), FAO-56 Eq. 39.
SIGMA = 4.903e-9
#: The same, per hour -- the hourly longwave uses the hour's own T^4.
SIGMA_HOUR = SIGMA / 24.0

#: ASCE/FAO-56 reference-ET numerator/denominator coefficients (Allen 2006).
#: The daily pair is the FAO-56 Eq. 6 constants; the hourly pair (short grass
#: reference) is what makes sub-daily met usable directly instead of a daily
#: total smeared back across the hours.
CN_DAILY, CD_DAILY = 900.0, 0.34
CN_HOURLY = 37.0
CD_HOURLY_DAY, CD_HOURLY_NIGHT = 0.24, 0.96
#: Reference-grass and open-water albedos.
ALBEDO_GRASS = 0.23
ALBEDO_WATER = 0.08
#: Daily-mean W/m^2 -> MJ/m^2/day (86400 s / 1e6).
WM2_TO_MJ_DAY = 0.0864


# ---------------------------------------------------------------------------
# Vapour pressure and its slope
# ---------------------------------------------------------------------------

def svp(t):
    """Saturation vapour pressure e0(T) [kPa] at temperature ``t`` [degC].

    FAO-56 Eq. 11.  ``svp(25) == 3.168``, ``svp(18) == 2.064``.
    """
    t = np.asarray(t, dtype=float)
    return 0.6108 * np.exp(17.27 * t / (t + 237.3))


def svp_slope(t):
    """Slope of the SVP curve Delta [kPa/degC] at ``t`` [degC].  FAO-56 Eq. 13.

    ``svp_slope(16.9) == 0.122``.
    """
    t = np.asarray(t, dtype=float)
    return 4098.0 * svp(t) / (t + 237.3) ** 2


def mean_svp(tmax, tmin):
    """Saturation vapour pressure es [kPa] = mean of e0(Tmax), e0(Tmin).

    FAO-56 Eq. 12 -- the mean of the two, **not** e0 of the mean temperature,
    which would bias es low.  ``mean_svp(25, 18) == 2.616``.
    """
    return (svp(tmax) + svp(tmin)) / 2.0


# ---------------------------------------------------------------------------
# Pressure and the psychrometric constant
# ---------------------------------------------------------------------------

def atm_pressure(elevation):
    """Atmospheric pressure P [kPa] at ``elevation`` [m].  FAO-56 Eq. 7.

    ``atm_pressure(100) == 100.1``; ``atm_pressure(1800) == 81.8``.
    """
    z = np.asarray(elevation, dtype=float)
    return 101.3 * ((293.0 - 0.0065 * z) / 293.0) ** 5.26


def psychrometric_constant(elevation=None, pressure=None):
    """Psychrometric constant gamma [kPa/degC].  FAO-56 Eq. 8.

    Give either ``elevation`` (m) or ``pressure`` (kPa).  At 100 m,
    gamma == 0.0666.
    """
    if pressure is None:
        if elevation is None:
            raise ValueError("give elevation or pressure")
        pressure = atm_pressure(elevation)
    return 0.665e-3 * np.asarray(pressure, dtype=float)


# ---------------------------------------------------------------------------
# Actual vapour pressure ea (three routes, best data first)
# ---------------------------------------------------------------------------

def ea_from_tdew(tdew):
    """Actual vapour pressure ea [kPa] from dew point [degC].  FAO-56 Eq. 14.

    The most reliable route when a dew point (or a humidity that converts to
    one) is measured.
    """
    return svp(tdew)


def ea_from_rhmean(rhmean, tmax, tmin):
    """ea [kPa] from mean relative humidity [%].  FAO-56 Eq. 19.

    ``ea = RHmean/100 * (e0(Tmax) + e0(Tmin)) / 2``.  The route used for the
    eThekwini feed, which reports a mean humidity (``relativeHumidity``).
    """
    rhmean = np.asarray(rhmean, dtype=float)
    return rhmean / 100.0 * mean_svp(tmax, tmin)


def ea_from_rhminmax(rhmax, rhmin, tmax, tmin):
    """ea [kPa] from RHmax and RHmin [%].  FAO-56 Eq. 17 (most accurate)."""
    rhmax = np.asarray(rhmax, dtype=float)
    rhmin = np.asarray(rhmin, dtype=float)
    return (svp(tmin) * rhmax / 100.0 + svp(tmax) * rhmin / 100.0) / 2.0


# ---------------------------------------------------------------------------
# Solar geometry and radiation
# ---------------------------------------------------------------------------

def day_of_year(date):
    """Day of year (1-366) from a date, datetime, or array of them."""
    return pd.DatetimeIndex(pd.to_datetime(np.atleast_1d(date))).dayofyear.to_numpy()


def _sun_geometry(latitude, doy):
    """Return ``(phi, dr, decl, ws)`` for FAO-56 Eqs. 22-25."""
    phi = np.radians(np.asarray(latitude, dtype=float))
    j = np.asarray(doy, dtype=float)
    dr = 1.0 + 0.033 * np.cos(2 * np.pi * j / 365.0)          # Eq. 23
    decl = 0.409 * np.sin(2 * np.pi * j / 365.0 - 1.39)       # Eq. 24
    # clip guards the poles/solstice where |tan phi tan decl| can exceed 1
    ws = np.arccos(np.clip(-np.tan(phi) * np.tan(decl), -1.0, 1.0))   # Eq. 25
    return phi, dr, decl, ws


def extraterrestrial_radiation(latitude, doy):
    """Extraterrestrial radiation Ra [MJ/m^2/day].  FAO-56 Eq. 21.

    Depends only on latitude and day of year -- pure solar geometry, no
    measurement.  ``Ra(-20 deg, 246) == 32.2``; ``Ra(50.8 deg, 187) == 41.1``.
    """
    phi, dr, decl, ws = _sun_geometry(latitude, doy)
    ra = (24 * 60 / np.pi) * GSC * dr * (
        ws * np.sin(phi) * np.sin(decl)
        + np.cos(phi) * np.cos(decl) * np.sin(ws))
    return ra


def daylight_hours(latitude, doy):
    """Daylight hours N.  FAO-56 Eq. 34.  ``N(-20 deg, 246) == 11.7``."""
    _, _, _, ws = _sun_geometry(latitude, doy)
    return 24.0 / np.pi * ws


def clear_sky_radiation(ra, elevation):
    """Clear-sky shortwave Rso [MJ/m^2/day].  FAO-56 Eq. 37."""
    z = np.asarray(elevation, dtype=float)
    return (0.75 + 2e-5 * z) * np.asarray(ra, dtype=float)


def rs_from_wm2(mean_wm2):
    """Daily-mean solar irradiance [W/m^2] -> Rs [MJ/m^2/day].

    ``Rs = mean_wm2 * 86400 / 1e6``.  The feed's confirmed unit is W/m^2, so a
    day's mean of, say, 255.6 W/m^2 is 22.1 MJ/m^2/day.
    """
    return np.asarray(mean_wm2, dtype=float) * WM2_TO_MJ_DAY


def net_shortwave(rs, albedo=ALBEDO_GRASS):
    """Net shortwave Rns [MJ/m^2/day].  FAO-56 Eq. 38."""
    return (1.0 - albedo) * np.asarray(rs, dtype=float)


def net_longwave(tmax, tmin, ea, rs, rso):
    """Net longwave Rnl [MJ/m^2/day] (outgoing positive).  FAO-56 Eq. 39.

    The ``Rs/Rso`` cloudiness factor is clipped to <= 1: a sensor reading
    slightly above the clear-sky estimate must not flip the term negative.
    """
    tmax = np.asarray(tmax, dtype=float)
    tmin = np.asarray(tmin, dtype=float)
    ea = np.asarray(ea, dtype=float)
    tmax_k4 = (tmax + 273.16) ** 4
    tmin_k4 = (tmin + 273.16) ** 4
    cloud = 1.35 * np.clip(np.asarray(rs, dtype=float) / rso, None, 1.0) - 0.35
    return (SIGMA * (tmax_k4 + tmin_k4) / 2.0
            * (0.34 - 0.14 * np.sqrt(np.clip(ea, 0, None))) * cloud)


def net_radiation(rs, tmax, tmin, ea, latitude, elevation, doy,
                  albedo=ALBEDO_GRASS):
    """Net radiation Rn [MJ/m^2/day].  FAO-56 Eq. 40 (Rn = Rns - Rnl).

    With the FAO-56 example inputs (Rs=22.07, Tmax=21.5, Tmin=12.3, ea=1.409,
    lat=50.8, z=100, doy=187) this returns 13.28.
    """
    ra = extraterrestrial_radiation(latitude, doy)
    rso = clear_sky_radiation(ra, elevation)
    rns = net_shortwave(rs, albedo)
    rnl = net_longwave(tmax, tmin, ea, rs, rso)
    return rns - rnl


# ---------------------------------------------------------------------------
# Wind
# ---------------------------------------------------------------------------

def wind_speed_2m(u, measured_height=2.0):
    """Wind speed adjusted to 2 m [m/s].  FAO-56 Eq. 47.

    ``measured_height`` is the anemometer height (m).  The eThekwini height is
    not documented in the feed; leave the default 2.0 (no adjustment) until it
    is confirmed, and set it (often 10.0 for an AWS) once known -- a 10 m speed
    is ~1.2x the 2 m speed, a real bias on the aerodynamic term.
    """
    u = np.asarray(u, dtype=float)
    h = float(measured_height)
    if h == 2.0:
        return u
    return u * 4.87 / np.log(67.8 * h - 5.42)


# ---------------------------------------------------------------------------
# Penman-Monteith reference ET
# ---------------------------------------------------------------------------

def penman_monteith_combine(delta, rn, g, gamma, tmean, u2, es_minus_ea,
                            cn=CN_DAILY, cd=CD_DAILY):
    """The FAO-56 / ASCE combination equation, given its assembled terms.

    The one equation serves both time steps; only the two reference
    coefficients change (Allen et al., 2006):

    ``cn=900, cd=0.34`` (**daily**, the defaults)
        ``rn``/``g`` in MJ/m^2/day, ``es_minus_ea`` from Tmax/Tmin, result
        mm/day.  Checked against the FAO-56 worked example: ``delta=0.122,
        rn=13.28, g=0, gamma=0.0666, tmean=16.9, u2=2.078, es_minus_ea=0.589``
        -> ``3.88`` mm/day.

    ``cn=37, cd=0.24`` day / ``0.96`` night (**hourly**)
        ``rn``/``g`` in MJ/m^2/hour, ``es_minus_ea`` from the hour's own
        temperature, result mm/hour.  ``cd`` (and the ``g`` the caller passes)
        switch on whether the hour is daylit; see
        :func:`et0_penman_monteith_hourly`.
    """
    delta = np.asarray(delta, dtype=float)
    gamma = np.asarray(gamma, dtype=float)
    u2 = np.asarray(u2, dtype=float)
    num = (0.408 * delta * (np.asarray(rn, dtype=float) - g)
           + gamma * (cn / (np.asarray(tmean, dtype=float) + 273.0))
           * u2 * np.asarray(es_minus_ea, dtype=float))
    den = delta + gamma * (1.0 + cd * u2)
    return num / den


def et0_penman_monteith(tmax, tmin, rs, u2, ea, elevation, latitude, doy,
                        tmean=None, g=0.0):
    """FAO-56 reference-crop ET0 [mm/day] from daily met.  FAO-56 Eq. 6.

    ``rs`` is in MJ/m^2/day (use :func:`rs_from_wm2` on the W/m^2 feed); ``u2``
    at 2 m (use :func:`wind_speed_2m`); ``ea`` from :func:`ea_from_rhmean` or
    :func:`ea_from_tdew`.  ``tmean`` defaults to (Tmax+Tmin)/2, the FAO-56
    convention.  ``g`` (soil heat flux) is ~0 for a daily step.

    With the FAO-56 worked-example inputs this returns 3.9 mm/day.
    """
    if tmean is None:
        tmean = (np.asarray(tmax, dtype=float) + np.asarray(tmin, dtype=float)) / 2.0
    delta = svp_slope(tmean)
    gamma = psychrometric_constant(elevation=elevation)
    es = mean_svp(tmax, tmin)
    rn = net_radiation(rs, tmax, tmin, ea, latitude, elevation, doy,
                       albedo=ALBEDO_GRASS)
    return penman_monteith_combine(delta, rn, g, gamma, tmean, u2, es - ea)


def et0_hargreaves(tmax, tmin, doy, latitude, tmean=None):
    """Hargreaves-Samani ET0 [mm/day] -- the temperature-only fallback.

    Hargreaves & Samani (1985): ``ET0 = 0.0023 (Tmean+17.8) sqrt(Tmax-Tmin)
    Ra_mm``, with Ra converted from MJ/m^2/day to mm/day by 0.408.  Needs only
    temperature plus solar geometry (Ra), so it runs on a station with no
    humidity, wind, or radiation.
    """
    tmax = np.asarray(tmax, dtype=float)
    tmin = np.asarray(tmin, dtype=float)
    if tmean is None:
        tmean = (tmax + tmin) / 2.0
    ra_mm = 0.408 * extraterrestrial_radiation(latitude, doy)
    dtr = np.clip(tmax - tmin, 0.0, None)      # guard a bad Tmax < Tmin day
    return 0.0023 * (tmean + 17.8) * np.sqrt(dtr) * ra_mm


# ---------------------------------------------------------------------------
# Open-water evaporation (ETo)
# ---------------------------------------------------------------------------

def e0_open_water(tmax, tmin, rs, u2, ea, elevation, latitude, doy, tmean=None):
    """Penman (1948) open-water evaporation E0 [mm/day] -- the channel's ``ETo``.

    Same energy + drying-power split as Penman-Monteith but for a free water
    surface: open-water albedo (0.08, so a little more net radiation than
    grass) and Penman's classic wind function in the aerodynamic term,
    ``Ea = 0.26 (1 + 0.54 u2)(es - ea)`` [mm/day].  Combined as
    ``E0 = (Delta*Rn_mm + gamma*Ea) / (Delta + gamma)``, with
    ``Rn_mm = 0.408*Rn``.

    Open water evaporates a little faster than grass, so E0 > ET0 under the same
    weather -- a useful sanity check.  When radiation is unavailable, the design
    note's fallback is ``ETo = ETr`` (channel evaporation is a sub-percent term),
    handled by the caller, not here.
    """
    tmax = np.asarray(tmax, dtype=float)
    tmin = np.asarray(tmin, dtype=float)
    if tmean is None:
        tmean = (tmax + tmin) / 2.0
    delta = svp_slope(tmean)
    gamma = psychrometric_constant(elevation=elevation)
    es = mean_svp(tmax, tmin)
    rn = net_radiation(rs, tmax, tmin, ea, latitude, elevation, doy,
                       albedo=ALBEDO_WATER)
    rn_mm = 0.408 * rn
    ea_term = 0.26 * (1.0 + 0.54 * np.asarray(u2, dtype=float)) * (es - ea)
    return (delta * rn_mm + gamma * ea_term) / (delta + gamma)


# ---------------------------------------------------------------------------
# Sub-daily radiation and reference ET (FAO-56 hourly / Allen 2006)
#
# When the met feed is already sub-daily -- the eThekwini network reports every
# half hour -- computing ET at the model step from those readings is strictly
# better than collapsing to a daily total and smearing it back across the hours
# by a clear-sky shape: the measured radiation, temperature, humidity and wind
# carry the real diurnal signal, clouds and all. That is the "use it explicitly"
# path. Disaggregation is reserved for genuinely daily input (see the routing in
# :func:`et0_eto_on_clock`).
# ---------------------------------------------------------------------------

def seasonal_correction(doy):
    """Seasonal correction for solar time Sc [hour].  FAO-56 Eqs. 32-33."""
    b = 2 * np.pi * (np.asarray(doy, dtype=float) - 81) / 364.0
    return 0.1645 * np.sin(2 * b) - 0.1255 * np.cos(b) - 0.025 * np.sin(b)


def extraterrestrial_radiation_hourly(latitude, longitude, doy, clock_hour,
                                       period_hours=1.0, tz_meridian=None):
    """Extraterrestrial radiation Ra for one sub-daily period [MJ/m^2/period].

    FAO-56 Eq. 28.  ``clock_hour`` is the **midpoint** of the period in local
    standard clock time; ``tz_meridian`` is the longitude (deg, east positive)
    of that clock's standard meridian -- 30 for SAST (UTC+2).  Pass
    ``tz_meridian=None`` to treat the clock as already solar (no longitude/
    equation-of-time shift), which is the right default for a naive clock.

    This doubles as the clear-sky **shape** for the daily-input disaggregation
    path: night periods return ~0, midday the most, from pure solar geometry.
    """
    phi = np.radians(np.asarray(latitude, dtype=float))
    j = np.asarray(doy, dtype=float)
    dr = 1.0 + 0.033 * np.cos(2 * np.pi * j / 365.0)
    decl = 0.409 * np.sin(2 * np.pi * j / 365.0 - 1.39)

    t = np.asarray(clock_hour, dtype=float)
    if tz_meridian is None:
        w = (np.pi / 12.0) * (t - 12.0)                       # clock == solar
    else:
        # FAO Eq. 31 in east-positive longitude: +(Lm - Lz)/15 hours.
        lm = np.asarray(longitude, dtype=float)
        shift = (lm - float(tz_meridian)) / 15.0 + seasonal_correction(j)
        w = (np.pi / 12.0) * ((t + shift) - 12.0)

    t1 = np.pi * float(period_hours) / 24.0
    w1, w2 = w - t1, w + t1
    ra = (12 * 60 / np.pi) * GSC * dr * (
        (w2 - w1) * np.sin(phi) * np.sin(decl)
        + np.cos(phi) * np.cos(decl) * (np.sin(w2) - np.sin(w1)))
    return np.clip(ra, 0.0, None)                             # night -> 0


def daily_cloudiness_factor(rs_day, ra_day, elevation):
    """Longwave cloudiness factor ``1.35 Rs/Rso - 0.35`` from **daily** totals.

    The hourly longwave needs the ``Rs/Rso`` ratio, but at night Rso is ~0 and
    the ratio is undefined.  FAO-56 handles this by carrying a daytime value
    through the night; taking the ratio from the day's totals does the same
    thing more simply and is stable, since cloudiness varies slowly.  The caller
    computes this once per station-day and applies it to every hour of the day.
    """
    rso_day = clear_sky_radiation(ra_day, elevation)
    ratio = np.clip(np.asarray(rs_day, dtype=float) / rso_day, None, 1.0)
    return np.clip(1.35 * ratio - 0.35, 0.05, 1.0)


def rs_rate_from_wm2(mean_wm2):
    """Period-mean solar irradiance [W/m^2] -> Rs **rate** [MJ/m^2/hour].

    ``Rs_rate = mean_wm2 * 3600 / 1e6``.  The hourly ET stack works in rates
    (per hour), so the model step ``Dt`` scales the depth once at the end,
    letting the very same code serve an hourly or a half-hourly ``Dt``.
    """
    return np.asarray(mean_wm2, dtype=float) * 3600.0 / 1e6


def net_radiation_hourly(rs, t, ea, fcd, albedo=ALBEDO_GRASS):
    """Net-radiation **rate** Rn [MJ/m^2/hour] for a sub-daily period.

    ``rs`` is the period-mean shortwave as a rate MJ/m^2/hour (use
    :func:`rs_rate_from_wm2` on the W/m^2 feed); ``fcd`` the cloudiness factor
    from :func:`daily_cloudiness_factor`.  The longwave uses the period's own
    temperature ``t`` and the hourly Stefan-Boltzmann constant.
    """
    rns = (1.0 - albedo) * np.asarray(rs, dtype=float)
    tk4 = (np.asarray(t, dtype=float) + 273.16) ** 4
    rnl = (SIGMA_HOUR * tk4
           * (0.34 - 0.14 * np.sqrt(np.clip(np.asarray(ea, dtype=float), 0, None)))
           * np.asarray(fcd, dtype=float))
    return rns - rnl


def et0_penman_monteith_hourly(t, rs, u2, ea, elevation, fcd, *,
                               period_hours=1.0):
    """FAO-56 hourly reference ET0 [mm/period] from one period's met.

    ``t`` the period-mean temperature [degC], ``rs`` the period-mean shortwave
    as a **rate** [MJ/m^2/hour] (:func:`rs_rate_from_wm2`), ``u2`` wind at 2 m
    [m/s], ``ea`` actual vapour pressure [kPa], ``fcd`` the daily cloudiness
    factor.  Day vs night (the ``Cd`` and soil-heat-flux switch) is decided from
    the sign of Rn -- positive by day, negative at night -- exactly as ASCE
    prescribes.  The hourly rate is scaled by ``period_hours`` (the model ``Dt``
    in hours), so a run of periods sums to the day's ET with no separate step.
    """
    t = np.asarray(t, dtype=float)
    rn = net_radiation_hourly(rs, t, ea, fcd, albedo=ALBEDO_GRASS)
    day = rn > 0
    g = np.where(day, 0.1 * rn, 0.5 * rn)
    cd = np.where(day, CD_HOURLY_DAY, CD_HOURLY_NIGHT)
    delta = svp_slope(t)
    gamma = psychrometric_constant(elevation=elevation)
    es = svp(t)                       # hourly: es is e0 at the hour's own T
    et = penman_monteith_combine(delta, rn, g, gamma, t, u2, es - ea,
                                 cn=CN_HOURLY, cd=cd)
    # A cold, dry, windy night can push the combination slightly below zero;
    # ET does not run backwards, so floor it. (~0 already, per FAO Ex. 19.)
    return np.clip(et * float(period_hours), 0.0, None)


# ---------------------------------------------------------------------------
# Convenience dispatcher
# ---------------------------------------------------------------------------

def reference_et0(row, latitude, measured_wind_height=2.0):
    """ET0 [mm/day] for one daily record, choosing the method from the data.

    ``row`` is a mapping/Series with the :mod:`met` daily-table fields
    (``tmax, tmin, tmean, date`` always; ``solar, wind, rh, tdew`` where the
    station has them).  Full Penman-Monteith is used when radiation **and** a
    humidity route (``rh`` or ``tdew``) are present; otherwise Hargreaves.

    Returns ``(et0, method)`` where method is ``"penman_monteith"`` or
    ``"hargreaves"``.  ``et0`` is NaN when even Tmax/Tmin are missing.
    """
    tmax, tmin = row.get("tmax"), row.get("tmin")
    tmean = row.get("tmean")
    doy = day_of_year(row["date"])[0]

    if pd.isna(tmax) or pd.isna(tmin):
        return float("nan"), "none"

    solar, wind = row.get("solar"), row.get("wind")
    rh, tdew = row.get("rh"), row.get("tdew")
    have_humidity = not pd.isna(rh) or not pd.isna(tdew)
    have_energy = not pd.isna(solar) and not pd.isna(wind)

    if have_energy and have_humidity:
        ea = (ea_from_tdew(tdew) if not pd.isna(tdew)
              else ea_from_rhmean(rh, tmax, tmin))
        rs = rs_from_wm2(solar)
        u2 = wind_speed_2m(wind, measured_wind_height)
        et0 = et0_penman_monteith(tmax, tmin, rs, u2, ea,
                                  row["elevation_m"], latitude, doy, tmean=tmean)
        return float(et0), "penman_monteith"

    et0 = et0_hargreaves(tmax, tmin, doy, latitude, tmean=tmean)
    return float(et0), "hargreaves"


# ---------------------------------------------------------------------------
# CLI -- daily ET0/E0 per station from a met daily table
# ---------------------------------------------------------------------------

def et0_eto_frame(daily, manifest, *, wind_height=2.0, method="auto"):
    """Compute daily ``et0`` (ETr) and ``eto`` (E0 open water) per station-day.

    ``daily`` is the :func:`met.daily_table` output (a DataFrame or its CSV path);
    ``manifest`` is the weather manifest (DataFrame or path) -- it supplies each
    station's ``elevation_m`` and, via its coordinates, latitude.  ``eto`` is the
    open-water Penman evaporation where the full energy balance is available and
    falls back to ``eto = et0`` otherwise (the design-note channel-evaporation
    shortcut).  Returns a frame ``station_id, date, et0, eto, method``.
    """
    from pyproj import Transformer
    from . import met as _met

    if not isinstance(daily, pd.DataFrame):
        daily = pd.read_csv(daily, parse_dates=["date"])
    daily = daily.copy()
    daily["station_id"] = daily["station_id"].astype(str).str.strip()

    man = manifest if isinstance(manifest, pd.DataFrame) else _met.read_manifest(manifest)
    # latitude (degrees) per station: reproject the model-CRS coordinates back
    # to lon/lat, so it is correct whatever CRS the manifest declared.
    tx = Transformer.from_crs(str(man["crs"].iloc[0]), "EPSG:4326", always_xy=True)
    _, lat = tx.transform(man["x"].to_numpy(float), man["y"].to_numpy(float))
    lat_of = dict(zip(man.index, lat))
    elev_of = man["elevation_m"].to_dict()

    rows = []
    for rec in daily.to_dict("records"):
        sid = rec["station_id"]
        if sid not in lat_of:
            continue
        rec["elevation_m"] = elev_of.get(sid, np.nan)
        latitude = lat_of[sid]

        forced_harg = method == "hargreaves"
        et0, used = reference_et0(rec, latitude, measured_wind_height=wind_height)
        if forced_harg and not pd.isna(rec.get("tmax")) and not pd.isna(rec.get("tmin")):
            doy = day_of_year(rec["date"])[0]
            et0 = float(et0_hargreaves(rec["tmax"], rec["tmin"], doy, latitude,
                                       tmean=rec.get("tmean")))
            used = "hargreaves"

        # open-water ETo: full Penman where energy + humidity exist, else = ETr
        eto = et0
        if used == "penman_monteith":
            ea = (ea_from_tdew(rec["tdew"]) if not pd.isna(rec.get("tdew"))
                  else ea_from_rhmean(rec["rh"], rec["tmax"], rec["tmin"]))
            eto = float(e0_open_water(
                rec["tmax"], rec["tmin"], rs_from_wm2(rec["solar"]),
                wind_speed_2m(rec["wind"], wind_height), ea,
                rec["elevation_m"], latitude, day_of_year(rec["date"])[0],
                tmean=rec.get("tmean")))

        rows.append({"station_id": sid, "date": rec["date"],
                     "et0": et0, "eto": eto, "method": used})

    return pd.DataFrame(rows, columns=["station_id", "date", "et0", "eto", "method"])


def _station_lonlat_elev(manifest):
    """Return ``{sid: (lon, lat, elevation_m)}`` from the weather manifest."""
    from pyproj import Transformer
    tx = Transformer.from_crs(str(manifest["crs"].iloc[0]), "EPSG:4326",
                              always_xy=True)
    lon, lat = tx.transform(manifest["x"].to_numpy(float),
                            manifest["y"].to_numpy(float))
    elev = manifest["elevation_m"].to_numpy(float)
    return {sid: (float(lo), float(la), float(el))
            for sid, lo, la, el in zip(manifest.index, lon, lat, elev)}


def _clearsky_shape(timeline, latitude, longitude, tz_meridian):
    """Clear-sky Ra at each step's midpoint -- the diurnal disaggregation shape."""
    ends = timeline.times
    mids = ends - timeline.dt / 2
    hours = mids.hour + mids.minute / 60.0
    ra = extraterrestrial_radiation_hourly(
        latitude, longitude, mids.dayofyear.to_numpy(),
        hours.to_numpy(), period_hours=timeline.dt_seconds / 3600.0,
        tz_meridian=tz_meridian)
    return pd.Series(ra, index=ends)


def et0_eto_on_clock(clean, manifest, timeline, *, method="auto",
                     wind_height=2.0, tz_meridian=None,
                     min_coverage=None):
    """Per-station ET0/ETo on the model clock, routed by data resolution.

    This is the resolution decision the rainfall side makes per gauge, applied
    to ET.  For each station, comparing its native sampling step to ``Dt``:

    * **native <= Dt and Dt is sub-daily** -> compute **hourly Penman-Monteith**
      at ``Dt`` from the step-aggregated met (:func:`et0_penman_monteith_hourly`).
      The sub-daily feed is used as measured; nothing is disaggregated.
    * **native > Dt and Dt is sub-daily** -> compute **daily** ET0
      (:func:`et0_eto_frame`) and split it across the hours with the clear-sky
      solar shape.  This is the *only* path that disaggregates.
    * **Dt is daily** -> daily ET0 mapped to each daily step, no disaggregation.

    Open-water ``eto`` on the hourly path is the hourly ``et0`` scaled by the
    day's open-water ratio ``E0/ET0`` (from the daily Penman open-water), so it
    keeps the sub-percent channel term consistent with the daily path and falls
    back to ``eto = et0`` on temperature-only days.

    Parameters
    ----------
    clean : DataFrame
        Cleaned long measurements from :func:`met.clean_measurements`.
    manifest : DataFrame
        From :func:`met.read_manifest`; its index is the station-column order.
    timeline : Timeline
        The one clock; ``tz`` set here is honoured throughout.
    tz_meridian : float, optional
        Standard-meridian longitude of the clock (30 for SAST).  Only shifts the
        clear-sky *shape* by a few minutes on the disaggregation path; leave
        ``None`` for a naive clock.

    Returns
    -------
    et0_readings, eto_readings : ndarray, shape (n_t, n_stations), mm per step
    available : ndarray of bool, shape (n_t, n_stations)
    report : DataFrame
        One row per station: the ``path`` taken and its native step.
    """
    from . import met as _met

    mc = _met.DEFAULT_MIN_COVERAGE if min_coverage is None else min_coverage
    tz = timeline.times.tz
    dt_hours = timeline.dt_seconds / 3600.0
    dt_is_subdaily = timeline.dt < pd.Timedelta("1D")
    station_ids = [str(s) for s in manifest.index]
    geo = _station_lonlat_elev(manifest)

    steps = _met.step_table(clean, timeline, min_coverage=mc)
    daily = _met.daily_table(clean, tz=tz, min_coverage=mc)
    daily_ee = et0_eto_frame(daily, manifest, wind_height=wind_height,
                             method=method)

    n_t = timeline.n_t
    et0 = np.zeros((n_t, len(station_ids)))
    eto = np.zeros((n_t, len(station_ids)))
    avail = np.zeros((n_t, len(station_ids)), dtype=bool)
    report = []

    for col, sid in enumerate(station_ids):
        native = _met.station_native_step(clean, sid)
        lon, lat, elev = geo.get(sid, (np.nan, np.nan, np.nan))
        subdaily_met = native is not None and native <= timeline.dt

        if dt_is_subdaily and subdaily_met:
            path = "subdaily"
            e0, eo, ok = _hourly_path(
                steps[steps["station_id"] == sid], daily_ee, sid, timeline,
                lat, elev, wind_height, dt_hours)
        else:
            path = "daily->disaggregated" if dt_is_subdaily else "daily"
            e0, eo, ok = _daily_path(
                daily_ee, sid, timeline, lat, lon, tz_meridian, dt_is_subdaily)

        et0[:, col], eto[:, col], avail[:, col] = e0, eo, ok
        report.append({"station_id": sid, "path": path,
                       "native_step": None if native is None else str(native)})

    return et0, eto, avail, pd.DataFrame(report)


def _daily_open_water_ratio(daily_ee, sid):
    """Per-date ``E0/ET0`` for one station; 1.0 where ET0 is 0/NaN/Hargreaves."""
    de = daily_ee[daily_ee["station_id"] == sid]
    r = {}
    for rec in de.to_dict("records"):
        e0 = rec["et0"]
        r[pd.Timestamp(rec["date"]).date()] = (
            rec["eto"] / e0 if pd.notna(e0) and e0 > 0 else 1.0)
    return r


def _hourly_path(st, daily_ee, sid, timeline, lat, elev, wind_height, dt_hours):
    """Hourly Penman-Monteith at Dt from step-aggregated met."""
    ends = timeline.times
    st = st.set_index("datetime").reindex(ends)
    temp = st["temp"].to_numpy(float)
    rh = st["rh"].to_numpy(float)
    wind = st["wind"].to_numpy(float)
    solar = st["solar"].to_numpy(float)

    ok = np.isfinite(temp) & np.isfinite(rh) & np.isfinite(wind) & np.isfinite(solar)

    # cloudiness factor per day, from the daily solar total vs clear-sky Ra
    doy = ends.dayofyear.to_numpy()
    ra_day = extraterrestrial_radiation(lat, doy)
    de = daily_ee[daily_ee["station_id"] == sid].set_index("date")
    # daily mean solar (W/m^2) -> MJ/m^2/day; pull from the daily table via met
    # is cleaner, but daily_ee lacks it, so recompute the day means from steps:
    day_key = pd.Index(ends.date)
    rs_day_mm = np.full(len(ends), np.nan)
    tmp = pd.Series(np.where(ok, solar, np.nan), index=ends)
    day_mean_wm2 = tmp.groupby(ends.date).transform("mean").to_numpy()
    rs_day_mm = rs_from_wm2(day_mean_wm2)
    fcd = daily_cloudiness_factor(rs_day_mm, ra_day, elev)

    et0 = np.zeros(len(ends))
    with np.errstate(invalid="ignore"):
        es = svp(temp)
        ea = rh / 100.0 * es
        rs_rate = rs_rate_from_wm2(solar)
        u2 = wind_speed_2m(wind, wind_height)
        vals = et0_penman_monteith_hourly(temp, rs_rate, u2, ea, elev, fcd,
                                           period_hours=dt_hours)
    et0[ok] = np.nan_to_num(vals[ok], nan=0.0)

    # open-water: scale the hour by the day's E0/ET0 ratio
    ratio_of = _daily_open_water_ratio(daily_ee, sid)
    r = np.array([ratio_of.get(d, 1.0) for d in ends.date])
    eto = et0 * r
    return et0, eto, ok


def _daily_path(daily_ee, sid, timeline, lat, lon, tz_meridian, dt_is_subdaily):
    """Daily ET0, then either disaggregate to sub-daily Dt or map to daily Dt."""
    from .gauges import disaggregate

    ends = timeline.times
    de = daily_ee[daily_ee["station_id"] == sid].copy()
    de = de[pd.notna(de["et0"])]
    if de.empty:
        return np.zeros(len(ends)), np.zeros(len(ends)), np.zeros(len(ends), bool)

    tz = ends.tz
    # daily total, interval-ending at the end of its calendar day
    day_end = pd.DatetimeIndex(
        [pd.Timestamp(d) + pd.Timedelta("1D") for d in de["date"]])
    if tz is not None:
        day_end = day_end.tz_localize(tz)
    et0_day = pd.Series(de["et0"].to_numpy(float), index=day_end)
    eto_day = pd.Series(de["eto"].to_numpy(float), index=day_end)

    if dt_is_subdaily:
        shape = _clearsky_shape(timeline, lat, lon, tz_meridian)
        et0_s = disaggregate(et0_day, timeline, shape=shape, step="1D")
        eto_s = disaggregate(eto_day, timeline, shape=shape, step="1D")
        e0 = et0_s.to_numpy(float)
        eo = eto_s.to_numpy(float)
    else:
        # daily Dt: each step carries the total for the day it covers
        covered = pd.Index((ends - timeline.dt).date)
        by_date_et0 = {d.date(): v for d, v in zip(day_end - pd.Timedelta("1D"),
                                                   et0_day.to_numpy())}
        by_date_eto = {d.date(): v for d, v in zip(day_end - pd.Timedelta("1D"),
                                                   eto_day.to_numpy())}
        e0 = np.array([by_date_et0.get(d, np.nan) for d in covered])
        eo = np.array([by_date_eto.get(d, np.nan) for d in covered])

    ok = np.isfinite(e0)
    return np.nan_to_num(e0, nan=0.0), np.nan_to_num(eo, nan=0.0), ok


def _cli(argv=None):
    import argparse

    p = argparse.ArgumentParser(
        description="Daily reference ET0 (ETr) and open-water E0 (ETo) per "
                    "station from a `met` daily table.",
    )
    p.add_argument("--daily", required=True,
                   help="weather_daily.csv from `met` (station_id, date, tmax, …)")
    p.add_argument("--manifest", required=True,
                   help="weather_manifest.csv (elevation_m + coordinates)")
    p.add_argument("--out", default=None, help="write per-station-day ET0/ETo CSV")
    p.add_argument("--wind-height", type=float, default=2.0,
                   help="anemometer height (m) for the 2 m adjustment "
                        "(default 2.0; set 10 for a typical AWS once confirmed)")
    p.add_argument("--method", choices=["auto", "hargreaves"], default="auto",
                   help="'auto' uses Penman-Monteith where radiation+humidity "
                        "allow, else Hargreaves; 'hargreaves' forces the "
                        "temperature-only method everywhere")
    args = p.parse_args(argv)

    out = et0_eto_frame(args.daily, args.manifest,
                        wind_height=args.wind_height, method=args.method)
    good = out[out["et0"].notna()]

    print("=== method used (station-days) ===")
    print(out["method"].value_counts().to_string())
    print("\n=== ET0 (mm/day) per station ===")
    summ = (good.groupby("station_id")["et0"]
            .agg(["count", "min", "mean", "max"]).round(2))
    print(summ.to_string())
    print("\n=== head ===")
    print(out.head(10).to_string(index=False))

    if args.out:
        out.to_csv(args.out, index=False)
        print(f"\n{len(out)} station-days -> {args.out}")


if __name__ == "__main__":
    _cli()

"""Temporal gap-filling for the ET demand series (tiers 1-2, no reanalysis yet).

Tier 1 -- lean on the neighbours -- already lives in the spatial step:
``etfields.build_and_write_etfields`` zeroes an absent station's weight column
for the steps it is out and renormalises the survivors, so a lone dropout is
already covered.  This module is **tier 2**: filling a station's *own* series in
time, so that a step where *every* station is down -- the case that trips the
writer's all-empty guard -- still has values to spread.

Two fills, chosen by how long the hole is:

* **short gap -> linear interpolation.**  ET0 is smooth hour to hour, so a hole
  of a few steps is well drawn by a straight line between the values that
  bracket it.
* **long gap -> the station's mean diurnal cycle.**  A straight line across a
  whole night would ramp through hours ET never reaches; instead each missing
  step takes the station's average value for that *hour of day* (and, by
  default, that *calendar month*), which keeps the day/night rhythm intact.

Both are built per station from that station's own available record.  A station
with no record at all in the window is left untouched for the spatial step to
drop.  ET0 is non-negative, so fills are floored at zero.

This is deliberately the light, no-download tier.  The heavier reanalysis
backfill (ERA5-Land, bias-corrected) is a later, separate piece.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

__all__ = [
    "DEFAULT_MAX_INTERP_GAP",
    "diurnal_climatology",
    "fill_series_in_time",
    "fill_readings_in_time",
]

#: Linear interpolation is used only across holes at most this many steps long;
#: longer holes fall to the diurnal climatology.  Three hours is short enough
#: that a line does not span a full night at an hourly Dt.
DEFAULT_MAX_INTERP_GAP = 3


def diurnal_climatology(values, available, times, *, by_month=True):
    """Return a lookup ``clim(i)`` giving the climatological value for step ``i``.

    The climatology is the mean of the **available** steps grouped by
    ``(month, hour-of-day)`` (or just ``hour-of-day`` when ``by_month`` is
    False).  Lookups fall back ``(month, hour) -> (hour) -> NaN``, so a month
    with no reading at some hour still borrows the all-month value for that hour.
    """
    v = np.asarray(values, dtype=float)
    ok = np.asarray(available, dtype=bool) & np.isfinite(v)
    hour = np.asarray(times.hour)
    month = np.asarray(times.month)

    if not ok.any():
        return lambda i: np.nan

    obs = pd.DataFrame({"v": v[ok], "h": hour[ok], "m": month[ok]})
    hh = obs.groupby("h")["v"].mean()
    mh = obs.groupby(["m", "h"])["v"].mean() if by_month else None

    def clim(i):
        if by_month:
            key = (month[i], hour[i])
            if key in mh.index:
                return float(mh.loc[key])
        h = hour[i]
        if h in hh.index:
            return float(hh.loc[h])
        return np.nan

    return clim


def _gap_runs(ok):
    """Yield ``(start, stop)`` half-open index ranges of contiguous False runs."""
    n = len(ok)
    i = 0
    while i < n:
        if ok[i]:
            i += 1
            continue
        j = i
        while j < n and not ok[j]:
            j += 1
        yield i, j
        i = j


def fill_series_in_time(values, available, times, *,
                        max_interp_gap=DEFAULT_MAX_INTERP_GAP, by_month=True):
    """Fill one station's series in time.  Returns ``(filled, available)``.

    Short holes (<= ``max_interp_gap`` steps) bracketed by observations on both
    sides are linearly interpolated; every other hole is filled from the diurnal
    climatology.  ``available`` is returned with the filled steps flipped True;
    a step the climatology cannot reach (an hour never observed) is left False
    for the spatial step.
    """
    v = np.asarray(values, dtype=float).copy()
    ok = np.asarray(available, dtype=bool).copy()
    if not ok.any():
        return v, ok                       # nothing to build from

    clim = diurnal_climatology(v, ok, times, by_month=by_month)
    n = len(v)

    for i, j in _gap_runs(ok):
        left = i - 1 if i - 1 >= 0 else None       # bracket indices (observed)
        right = j if j < n else None
        length = j - i
        if length <= max_interp_gap and left is not None and right is not None:
            span = right - left
            for k in range(i, j):
                w = (k - left) / span
                v[k] = (1.0 - w) * v[left] + w * v[right]
                ok[k] = True
        else:
            for k in range(i, j):
                c = clim(k)
                if np.isfinite(c):
                    v[k] = c
                    ok[k] = True

    np.clip(v, 0.0, None, out=v)            # ET does not run backwards
    return v, ok


def fill_readings_in_time(readings, available, timeline, *,
                          max_interp_gap=DEFAULT_MAX_INTERP_GAP, by_month=True):
    """Fill every station column in time.

    Returns ``(filled_readings, filled_available, n_filled)`` where ``n_filled``
    is how many station-steps were turned from gap into value -- the number the
    caller reports.
    """
    r = np.asarray(readings, dtype=float).copy()
    a = np.asarray(available, dtype=bool).copy()
    times = timeline.times
    n_gap_before = int((~a).sum())
    for c in range(r.shape[1]):
        r[:, c], a[:, c] = fill_series_in_time(
            r[:, c], a[:, c], times,
            max_interp_gap=max_interp_gap, by_month=by_month,
        )
    n_filled = n_gap_before - int((~a).sum())
    return r, a, n_filled

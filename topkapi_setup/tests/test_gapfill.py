"""Tests for the tier-2 temporal gap-fill (forcing/gapfill.py).

Pin the two behaviours that matter: a short hole is drawn straight across, a
long hole takes the station's own diurnal shape (not a line ramping through the
night), and a column with nothing in it is left for the spatial step.
"""

import numpy as np
import pandas as pd
import pytest

from topkapi_setup.forcing import gapfill as gf
from topkapi_setup.forcing.gauges import Timeline


def _times(n, start="2025-01-02 00:00"):
    return pd.date_range(start, periods=n, freq="h")


def test_short_gap_is_linearly_interpolated():
    t = _times(6)
    v = np.array([1.0, np.nan, np.nan, 4.0, 5.0, 6.0])
    ok = np.array([True, False, False, True, True, True])
    filled, av = gf.fill_series_in_time(v, ok, t, max_interp_gap=3)
    assert av.all()
    assert filled[1] == pytest.approx(2.0)      # 1 -> 4 over three steps
    assert filled[2] == pytest.approx(3.0)


def test_long_gap_uses_diurnal_climatology_not_a_line():
    # three days hourly; a diurnal bump (0 at night, peak midday)
    t = _times(72)
    h = t.hour.to_numpy()
    shape = np.clip(np.sin((h - 6) / 12 * np.pi), 0, None)   # day bump
    v = shape.copy()
    ok = np.ones(72, bool)
    ok[24:48] = False                                        # blank the middle day
    filled, av = gf.fill_series_in_time(v, ok, t, max_interp_gap=3)
    assert av.all()
    # the filled midday hour ~ the mean of the other days' midday, not a ramp
    midday = np.where(h == 12)[0]
    assert filled[24 + 12] == pytest.approx(v[midday[0]], abs=1e-6)
    # a filled night hour stays ~0 (a straight line would have lifted it)
    assert filled[24 + 3] < 0.05


def test_empty_column_left_for_the_spatial_step():
    t = _times(5)
    v = np.full(5, np.nan)
    ok = np.zeros(5, bool)
    filled, av = gf.fill_series_in_time(v, ok, t)
    assert not av.any()                          # nothing to build from


def test_fill_readings_reports_count_and_stays_nonnegative():
    t = _times(12)
    tl = Timeline(str(t[0]), str(t[-1]), 3600)
    r = np.tile(np.linspace(0, 2, 12)[:, None], (1, 2))
    a = np.ones((12, 2), bool)
    a[5:8, 0] = False                            # a 3-step hole in station 0
    a[2, 1] = False                              # a 1-step hole in station 1
    filled, av, n = gf.fill_readings_in_time(r, a, tl, max_interp_gap=3)
    assert n == 4 and av.all()
    assert (filled >= 0).all()

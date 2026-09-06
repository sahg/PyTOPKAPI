"""Third-party compatibility shims.

pysheds (latest release 0.5, Aug 2025) calls ``numpy.in1d``. NumPy deprecated
``in1d`` in 2.0 and *removed* it in 2.4.0 (superseded by ``numpy.isin``), so on
NumPy >= 2.4 pysheds raises ``AttributeError`` on import of the terrain stage
until the alias is restored. This shim restores it; it is a no-op on NumPy
versions that still ship ``in1d`` (< 2.4).

Import this module *before* importing pysheds.

This is a bridge, not a fix -- it exists only because no released pysheds yet
supports NumPy >= 2.4. Durable options, in increasing order of effort, for when
you want to delete it:

  1. Pin a pysheds release that supports NumPy >= 2.4, once one ships
     (watch https://github.com/mdbartos/pysheds/releases).
  2. Vendor the one-line ``in1d`` -> ``isin`` change into a local pysheds.
  3. Move terrain analysis to ``richdem`` (the alternative noted in
     conda-env.yaml).

Revisit at the M0 environment refresh.
"""

import numpy as np

if not hasattr(np, "in1d"):
    # numpy.isin is the documented drop-in replacement; its signature matches
    # pysheds' usage (a 1-D membership test).
    np.in1d = np.isin  # type: ignore[attr-defined]

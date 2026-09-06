"""Forcing builders (M3).

One array is the target: rainfall depth of shape ``(n_t, n_cells)``, on the
canonical cell order, written to ``rainfields.h5``.  Everything in this package
is machinery for filling it.

Modules land in build order:

===================  ============================================  ========
``interpolate.py``   the ``W``-builders (mean/thiessen/idw/kriging)  present
``gauges.py``        manifest + measurements IO, clock, gaps         present
``rainfields.py``    the ``(n_t, n_cells)`` writer                   present
``build.py``         end-to-end CLI (+ ``forcing_manifest.json``)    present
``gridded.py``       CHIRPS/IMERG resample onto grid                 to do
===================  ============================================  ========

Run the whole stage with ``python -m topkapi_setup.forcing`` (see ``build.py``).
"""

from topkapi_setup.forcing.gauges import (  # noqa: F401
    Timeline,
    aggregate,
    align_to_clock,
    coverage,
    disaggregate,
    gauge_xy,
    read_manifest,
    read_measurements,
)
from topkapi_setup.forcing.rainfields import (  # noqa: F401
    build_and_write_rainfields,
    build_rainfield,
    check_cell_order,
    read_rainfields,
    write_rainfields,
)
from topkapi_setup.forcing.etfields import (  # noqa: F401
    build_and_write_etfields,
    build_etfields,
    read_etfields,
    write_etfields,
)
from topkapi_setup.forcing.interpolate import (  # noqa: F401
    DEFAULT_BUFFER_M,
    DEFAULT_IDW_POWER,
    DEFAULT_METHOD,
    METHODS,
    build_weights,
    catchment_cell_xy,
    renormalise,
    select_gauges,
)
from topkapi_setup.forcing.build import (  # noqa: F401
    ForcingResult,
    build_rainfields,
)

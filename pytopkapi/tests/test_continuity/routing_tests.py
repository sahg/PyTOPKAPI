"""Regression tests for inter-cell routing.

These tests guard against a class of bug where the parallel execution
path silently drops the flow routed between cells, so that each cell is
solved as though it were an isolated headwater. Mass balance (continuity)
still holds in that broken state -- water that fails to route simply
stays in the upstream cell stores -- so the continuity tests do not catch
it. The check here is different: the parallel and serial solvers must
produce the *same* outlet hydrograph.

The 4-cell catchment (Parak, 2006) is used as the fixture because three
of its four cells have upstream contributors, so the outlet flow is
dominated by routed water rather than by local rainfall or initial
conditions.
"""

import os
import contextlib

from configparser import ConfigParser

import h5py
import numpy as np

import pytopkapi

# Directory holding this test module and its ini/fixture files. Resolving
# paths against __file__ keeps the test independent of the caller's
# working directory, so it behaves the same under pytest, nose, or when
# run directly.
TEST_DIR = os.path.dirname(os.path.abspath(__file__))


@contextlib.contextmanager
def _working_dir(path):
    """Temporarily change the working directory.

    The ini files reference their inputs by relative path, so the model
    must run with the test directory as the working directory.
    """
    original = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(original)


def _run_and_read_channel_flow(ini_fname, parallel_exec):
    """Run a simulation and return the channel outflow array.

    The result file named in the ini is removed before returning so the
    two execution modes cannot read each other's output.
    """
    with _working_dir(TEST_DIR):
        pytopkapi.run(ini_fname,
                      parallel_exec=parallel_exec,
                      nworkers=1,
                      quiet=True)

        config = ConfigParser()
        config.read(ini_fname)
        result_fname = config.get('output_files', 'file_out')

        with h5py.File(result_fname, 'r') as h5file:
            Qc_out = h5file['/Channel/Qc_out'][...]

        os.remove(result_fname)

    return Qc_out


def test_parallel_matches_serial_4cell():
    """Parallel and serial solvers agree on the 4-cell catchment.

    A tolerance of 1e-6 m3/s is far tighter than the ~0.3 m3/s outlet
    discrepancy produced when routing is dropped, so this fails loudly
    on the routing bug while tolerating ordinary floating-point noise.
    """
    ini_fname = '4cells.ini'
    outlet_cell = 3

    old_settings = np.seterr(all='ignore')
    try:
        Qc_serial = _run_and_read_channel_flow(ini_fname, parallel_exec=False)
        Qc_parallel = _run_and_read_channel_flow(ini_fname, parallel_exec=True)
    finally:
        np.seterr(**old_settings)

    # Guard against a trivial pass: the outlet must actually carry routed
    # flow, well above its initial condition. If this ever fails, the
    # fixture has stopped exercising routing and the equality check below
    # would be meaningless.
    peak = Qc_serial[:, outlet_cell].max()
    initial = Qc_serial[0, outlet_cell]
    assert peak > 10 * initial, \
        'fixture no longer routes flow to the outlet (peak=%g)' % peak

    # The two execution modes must produce the same hydrograph everywhere.
    max_abs_diff = np.abs(Qc_serial - Qc_parallel).max()
    assert max_abs_diff < 1e-6, \
        'parallel and serial channel flow disagree by %g m3/s' % max_abs_diff


if __name__ == '__main__':
    test_parallel_matches_serial_4cell()
    print('routing regression test PASSED')

"""Wrapper script to set-up proper environment to run GRASS

Before importing the grass package, the correct environment variables
must be set. This script sets up the environment for a Linux system,
then calls the processing script in a subprocess.

"""
import os
import subprocess

def prepend_env_var(key, val):
    """Prepend `val` to an environment variable

    A new variable is created if `key` doesn't exist, otherwise `val`
    is prepended to the existing environment variable.

    """
    if key in os.environ.keys():
        old_val = os.environ[key]

        os.environ[key] = ':'.join([val, old_val])
    else:
        os.environ[key] = val

def find_gisbase():
    """Find the value of the GISBASE environ variable on Linux

    A bit hackish, but seems to work for now..

    """
    gscript = subprocess.check_output(['which', 'grass']).strip('\n')

    fp = open(gscript, 'r')
    for line in fp:
        if 'GISBASE=' in line:
            gis_base = line.split('=')[1]
    fp.close()

    return gis_base.strip() # remove line-endings

gis_base = find_gisbase()
home = os.environ['HOME']

prepend_env_var('GISBASE', '%s' % gis_base)
prepend_env_var('PATH', '%s/bin:%s/scripts' % (gis_base, gis_base))
prepend_env_var('LD_LIBRARY_PATH', '%s/lib' % gis_base)
prepend_env_var('PYTHONPATH', '%s/etc/python' % gis_base)
prepend_env_var('GISRC', '%s/.grassrc6' % home)
prepend_env_var('GIS_LOCK', '%s' % os.getpid())

cmd = ['python', 'process-catchment.py']
subprocess.call(cmd)

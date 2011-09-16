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

# To Do: Find out how to discover the GRASS installation directory and
# avoid hard-coding this
gis_base = '/usr/lib/grass64'
home = os.environ['HOME']

prepend_env_var('PATH', '%s/bin:%s/scripts' % (gis_base, gis_base))
prepend_env_var('PYTHONPATH', '%s/etc/python' % gis_base)
prepend_env_var('GISRC', '%s/.grassrc6' % home)

cmd = ['python', 'process-catchment.py']
subprocess.call(cmd)

import os
import subprocess

# Before importing the grass package, set environment variables for a
# Linux system
def prepend_env_var(key, val):
    """Prepend `val` to an environment variable

    A new variable is created if `key` doesn't exist.

    """
    if key in os.environ.keys():
        print('Key %s exists:' % key)
        print(os.environ[key])

        old_val = os.environ[key]

        os.environ[key] = ':'.join([val, old_val])
    else:
        os.environ[key] = val

    print('Edited %s:' % key)
    print(os.environ[key])

gis_base = '/usr/lib/grass64'
home = os.environ['HOME']

prepend_env_var('PATH', '%s/bin:%s/scripts' % (gis_base, gis_base))
prepend_env_var('PYTHONPATH', '%s/etc/python' % gis_base)
prepend_env_var('GISRC', '%s/.grassrc6' % home)

cmd = ['python', 'process-catchment.py']
subprocess.call(cmd)

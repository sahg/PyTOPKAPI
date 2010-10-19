.. _install:

===================================
PyTOPKAPI installation instructions
===================================

Download
--------

PyTOPKAPI source releases and a 32-bit Windows installer can be
downloaded from http://github.com/sahg/PyTOPKAPI/downloads

Prerequisites
-------------

Before using the PyTOPKAPI package, you will need a working Python
installation on your computer. Python can be obtained from
http://www.python.org if it is not already present on your system. In
addition, the model depends on a number of 3rd party Python packages,
these are: `NumPy <http://www.numpy.org>`_, `SciPy
<http://www.scipy.org>`_ and `PyTables
<http://www.pytables.org>`_. All of the 3rd party packages must be
installed and working on your system for PyTOPKAPI to work correctly.

Windows Binary Installer
------------------------

Simply double-click the installer and follow the instructions.

Source Install
--------------

The PyTOPKAPI package currently consists of pure Python code and does
not require any extension modules to be built, so installation from
source is straightforward on all popular operating systems. Windows
users who have installed PyTOPKAPI using the binary installer do not
need to perform a source install.

Once the prerequisite Python packages have been successfully
installed, PyTOPKAPI can be installed by issuing the following command
at the system command line::

    $ python setup.py install

If this fails, you'll need to check that the `python` executable is on
your system path.

Next steps
----------

Once installed a basic test of the installation can be made in the
Python interpreter::

    >>> import pytopkapi
    >>> print pytopkapi.__version__

New users should then follow the basic :ref:`tutorial <tutorial>`.

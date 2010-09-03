.. _install:

===================================
PyTOPKAPI installation instructions
===================================

Download
--------

Source releases and 32-bit Windows installers of the PyTOPKAPI code
can be downloaded from http://github.com/sahg/PyTOPKAPI/downloads

Install
-------

The PyTOPKAPI package currently consists of pure Python code and does
not require any extension modules to be built. However the model
depends on a number of 3rd party Python packages, these are: NumPy
(http://www.numpy.org), SciPy (http://www.scipy.org) and PyTables
(http://www.pytables.org).

Once the pre-requisite Python packages have been successfully
installed, PyTOPKAPI can be installed at the command line for any
operating system::

    $ python setup.py install

The installation can be tested in the Python interpreter::

    >>> import pytopkapi

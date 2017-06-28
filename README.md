![PyTOPKAPI logo](http://sahg.github.io/PyTOPKAPI/static/logo.png)

PyTOPKAPI is a BSD licensed Python library implementing the TOPKAPI
Hydrological model (Liu and Todini, 2002). The model is a
physically-based and fully distributed hydrological model, which has
already been successfully applied in several countries around the
world (e.g. Liu and Todini, 2002; Bartholomes and Todini, 2005; Liu et
al., 2005; Martina et al., 2006; Vischel et al., 2008a/b).

Install
-------

Download the latest release
[![DOI](https://zenodo.org/badge/689231.svg)](https://zenodo.org/badge/latestdoi/689231)
Alternatively available on
[Github](https://github.com/sahg/PyTOPKAPI/releases)

The PyTOPKAPI package currently consists of pure Python code and does
not require any extension modules to be built. However, the model
depends on a number of 3rd party Python packages with compiled code
linking to external libraries. These are most easily managed using the
conda package manager from
[Anaconda](https://docs.continuum.io/anaconda/).

Once Anaconda is installed, a conda environment with the pre-requisite
Python packages can be created and activated [this should be done from
the root of the unpacked PyTOPKAPI archive]:

    $ conda env create -f conda-env.yaml
    $ source activate pytopkapi-0.4.x

PyTOPKAPI can be installed into the conda environment using pip [note
the 'dot' at the end, which refers to the current directory]:

    $ pip install --upgrade .

The package tests can be run as follows:

    $ nosetests

An example simulation can be run from the *example_simulation*
directory:

    $ cd example_simulation
    $ python run_example.py

Documentation
-------------

The most recent documentation for PyTOPKAPI can be found at
http://sahg.github.io/PyTOPKAPI/


Citation
--------

Latest release:

[![DOI](https://zenodo.org/badge/689231.svg)](https://zenodo.org/badge/latestdoi/689231)

Selected publications:

Sinclair S. and Pegram G.G.S., (2010), "A comparison of ASCAT and
modeled soil moisture over South Africa, using TOPKAPI in land surface
mode", Hydrol. Earth Syst. Sci., 14, 613-626. [Full Text via
CrossRef](http://dx.doi.org/10.5194/hess-14-613-2010)

Sinclair S. and Pegram G.G.S., (2013), "A sensitivity assessment of
the TOPKAPI model with an added infiltration module", J. Hydrol., 479,
100-112. [Full Text via
CrossRef](http://dx.doi.org/10.1016/j.jhydrol.2012.11.061)

Vischel T., Pegram G.G.S., Sinclair S. and Parak M., (2008a),
"Implementation of the TOPKAPI model in South Africa: Initial results
from the Liebenbergsvlei catchment", Water SA, 34(3), 1-12. [Full
Text](http://hdl.handle.net/10520/EJC116536)

Vischel T., Pegram G.G.S., Sinclair S., Wagner W. and Bartsch A.,
(2008b), "Comparison of soil moisture fields estimated by catchment
modelling and remote sensing: a case study in South Africa",
Hydrol. Earth Syst. Sci., 12, 751-767. [Full Text via
CrossRef](http://dx.doi.org/10.5194/hess-12-751-2008)
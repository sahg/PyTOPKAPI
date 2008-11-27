Introduction
============

This section introduces the TOPKAPI model as it has been implemented and
presents the general layout of the documentation.

Generalities
------------

The TOPKAPI_SA  package
The file TOPKAPI_SA.zip contains:

* Software_package: a directory containing the Python installation executable and the Python packages required by the TOPKAPI model, as well as some other useful software utilities.
* TOPKAPI_python_package.exe: the package containing the TOPKAPI program.
* TOPKAPI_Example: a directory containing an example to run the model.

Conditions of use of the TOPKAPI package
The use of the TOPKAPI package is free, however it is the property of the developers of the model: Théo Vischel, Geoff Pegram, Scott Sinclair. Under this condition:
Users are required to inform the authors about their use and application of the TOPKAPI package.
The developers must be associated or acknowledged in any publications that would use the codes of the TOPKAPI package.
Comments and proposition of improvements of the model are welcomed, as well as students and researchers who would like to join the developer team!

Objective of the manual
This manual is aiming at giving to the future users of TOPKAPI the keys to handle the TOPKAPI package. It is not a description of the theoretical aspects of the TOPKAPI model, nor a detailed description of the modelling methodology, although some methodological aspects are addressed in the manual. For an all-comprehensive description of the theory and methodology associated with the TOPKAPI model, readers will refer to the following references: Liu and Todini (2002); Pegram et al. (2007); Vischel et al. (2007a) and Vischel et al. (2007b, submitted to Water SA).
There are three main steps in the use of the TOPKAPI model. The first step is the creation of the parameter files of the models according mainly to the catchment data, as well as the creation of the forcing variables (rain and reference crop evapotranspiration fields). The second step is to simply run the model with the created inputs. The third step is the analysis of the simulations by plotting graphics. 
Nota Bene
The extraction of catchment data (mainly by using GIS) as well as the creation of forcing fields is left to the initiative of the users who are assumed to make use of their own expertise, software and numerical tools to do it. 
Numerically speaking, efforts have been made by the developers in coding and organising properly the TOPKAPI model itself, so that is can be easily run once the parameter files are created. 
Scientifically speaking, it is clear that the most crucial step is the creation of the parameter files, and particularly the way they are used to link the data available on the catchment, to the effective physical parameters required by TOPKAPI. This step requires knowledge of hydrology and must be conducted under the expertise of the user. However, in order to facilitate the process of using TOPKAPI, a module called parameter_utils was included in the TOPKAPI package. This module is only a suggestion from the developers to help future users in the creation of the parameter files starting from the catchment data. However, since this module was created for the specific case study on the Liebenbergsvlei catchment (4625 km2, South Africa) using the personal choices of the developers, it is to be used with care and requires from the user to carefully read and sometimes modify the Python codes if different applications or scientific choices would have to be made.
The analysis of the results of the TOPKAPI model is reduced to a simple analysis of the simulations, so that the user can rapidly evaluate the quality of the run.

Overview of the manual
The manual is organised as follows. The second section describes the installation of the model. In the third section are described the TOPKAPI Python package and the example of application used for the tutorial that is presented in Sections 4, 5 & 6. This tutorial mainly explains how to (i) run the model, (ii) analyse the simulation results and (iii) create and modify the parameter files. 

2.Installation of the TOPKAPI model
TOPKAPI was coded in the Python programming language. Python is a package of free software (see http://www.python.org for detailed information) that was chosen, primarily since it is suited to the management of large data arrays and because it can be extended and interfaced with many other existing languages (FORTRAN, C, C++, etc.).
The minimal version of Python to be used for installing the TOPKAPI_python_package is Python 2.5. The following packages are required: Numpy (at least version 1.0.3), Scipy (at least version 0.5.2), Tables (at least version 2.0), Matplotlib (at least version 0.90.1)
For a convenient command shell the authors recommend the use of IPython (at least version 0.8.1) that requires pyreadline (at least version 1.4.4).
Python and required package executables are available in the directory software_package and must be installed in the following order:
For Python:
Execute (double click) python-2.5.1.msi
For the required packages:
Execute (double click) numpy-1.0.3.1.win32-py2.5.exe
Execute scipy-0.5.2.1.win32-py2.5.exe
Execute tables-2.0.win32-py2.5.exe
Execute matplotlib-0.90.1.win32-py2.5.exe
Execute TOPKAPI_python_package.exe
For the IPython command shell
Execute pyreadline-1.4.4.win32.exe
Execute ipython-0.8.1.win32.exe
Some specific formats called HDF5 formats are used in the TOPKAPI package for their convenience in compressing and archiving the data (see http://hdf.ncsa.uiuc.edu/products/hdf5/index.html for more information). Several programs exist to navigate into the HDF5 files through user-friendly interfaces. 
The authors recommend the use of HDF_Explorer, which is provided in the Software_package and can be installed by: 
Executing HDF_ExploreSetup.exe.

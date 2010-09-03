============
Introduction
============

This manual is aiming at giving to the future users of TOPKAPI the
keys to handle the TOPKAPI package. It is not a description of the
theoretical aspects of the TOPKAPI model, nor a detailed description
of the modelling methodology, although some methodological aspects are
addressed in the manual. For an all-comprehensive description of the
theory and methodology associated with the TOPKAPI model, readers will
refer to the following references: Liu and Todini (2002); Pegram et
al. (2007); Vischel et al. (2007a) and Vischel et al. (2007b,
submitted to Water SA).

There are three main steps in the use of the TOPKAPI model. The first
step is the creation of the parameter files of the models according
mainly to the catchment data, as well as the creation of the forcing
variables (rain and reference crop evapotranspiration fields). The
second step is to simply run the model with the created inputs. The
third step is the analysis of the simulations by plotting graphics.

Nota Bene
~~~~~~~~~

* The extraction of catchment data (mainly by using GIS) as well as
  the creation of forcing fields is left to the initiative of the
  users who are assumed to make use of their own expertise, software
  and numerical tools to do it.

* Numerically speaking, efforts have been made by the developers in
  coding and organising properly the TOPKAPI model itself, so that is
  can be easily run once the parameter files are created.

* Scientifically speaking, it is clear that the most crucial step is
  the creation of the parameter files, and particularly the way they
  are used to link the data available on the catchment, to the
  effective physical parameters required by TOPKAPI. This step
  requires knowledge of hydrology and must be conducted under the
  expertise of the user. However, in order to facilitate the process
  of using TOPKAPI, a module called parameter_utils was included in
  the TOPKAPI package. This module is only a suggestion from the
  developers to help future users in the creation of the parameter
  files starting from the catchment data. However, since this module
  was created for the specific case study on the Liebenbergsvlei
  catchment (4625 km2, South Africa) using the personal choices of the
  developers, it is to be used with care and requires from the user to
  carefully read and sometimes modify the Python codes if different
  applications or scientific choices would have to be made.

* The analysis of the results of the TOPKAPI model is reduced to a
  simple analysis of the simulations, so that the user can rapidly
  evaluate the quality of the run.

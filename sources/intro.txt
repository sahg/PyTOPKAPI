============
Introduction
============

This documentation aims at providing sufficient information to allow
users of the PyTOPKAPI package to set-up and run the model for their
own purposes. It is not intended to be a definitive description of the
theoretical aspects of the TOPKAPI model, nor a detailed description
of the modelling methodology, although some methodological aspects are
addressed. For a comprehensive description of the theory and
methodology associated with the TOPKAPI model, readers should refer to
the following references: Liu and Todini (2002), Pegram et al. (2006),
Vischel et al. (2008a) and Vischel et al. (2008b).

There are three basic steps in the use of the PyTOPKAPI package. The
first step is the creation of one or more model parameter files
derived from the available catchment data, as well as the associated
forcing variables (rainfall and evapotranspiration fields). The second
step is to simply run the model using the designated parameter file
and forcing variables. The third step is the analysis of the resulting
simulations. The PyTOPKAPI package provides tools to assist with each
of these tasks, but the user does need to exercise experience and
judgement:

* The extraction of catchment data (mainly using GIS) as well as the
  creation of forcing fields is left to the initiative of the users
  who are assumed to make use of their own expertise, software and
  numerical tools.

* Numerically speaking, efforts have been made by the developers in
  coding and organising the PyTOPKAPI package, so that it can be
  easily run once the parameter files have been created.

* Scientifically speaking, it is clear that the most crucial step is
  the creation of the parameter files, and particularly the way they
  are used to link the data available on the catchment, to the
  effective physical parameters required by the TOPKAPI model. This
  step requires knowledge of hydrology and must be conducted under the
  expertise of the user. However, in order to facilitate the process
  of using the PyTOPKAPI package, a module called `parameter_utils` is
  included in PyTOPKAPI. This module is only a suggestion from the
  developers to help future users in the creation of the parameter
  files starting from the catchment data. However, since this module
  was created based on the developers needs, it should be used with
  care and requires the user to carefully read and sometimes modify
  the Python codes if different applications or scientific choices are
  made.

* The analysis of the results of the PyTOPKAPI package is reduced to a
  simple analysis of the simulations, so that the user can rapidly
  evaluate the quality of the run.

References
----------

**Liu Z. and Todini E., 2002**, Towards a comprehensive physically-based
rainfall-runoff model, Hydrol. Earth Syst. Sci., 6(5), 859 – 881.

**Pegram G., Sinclair S., Parak M., Sakulski D. and Nxumalo N.,
2006**, National Flood Nowcasting System: Towards an integrated
mitigation strategy, WRC Report No. 1429/1/06, Water Research
Commission, Pretoria, South Africa.

**Vischel T., Pegram G.G.S., Sinclair S., Wagner W. and Bartsch A.,
2008a**, Comparison of soil moisture fields estimated by catchment
modelling and remote sensing: a case study in South Africa,
Hydrol. Earth Syst. Sci., 12, 751-767.

**Vischel T., Pegram G.G.S., Sinclair S.  and Parak M., 2008b**,
Implementation of the TOPKAPI model in South Africa: Initial results
from the Liebenbergsvlei catchment, Water SA, 34(3), 1-12.

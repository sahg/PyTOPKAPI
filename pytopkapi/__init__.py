"""
PyTOPKAPI is a BSD licensed Python package implementing the TOPKAPI
Hydrological model (Liu and Todini, 2002). The model is a
physically-based and fully distributed hydrological model, which has
already been successfully applied in several countries around the
world (Liu and Todini, 2002; Bartholomes and Todini, 2005; Liu et al.,
2005; Martina et al., 2006; Vischel et al., 2008, Sinclair and Pegram,
2010).

"""

try:
    from __dev_version import version as __version__
    from __dev_version import git_revision as __git_revision__
except ImportError:
    from __version import version as __version__
    from __version import git_revision as __git_revision__

import model
from model import *
import results_analysis
from results_analysis import *

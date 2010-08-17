"""Package providing an implementation of the TOPKAPI model and some utilities.

The interface isn't stable yet so be prepared to update your code
on a regular basis...

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

# This test works on a Windows XP system running:
# Python 2.5.1
# numpy.__version__ '1.0.3'
# scipy.__version__ '0.5.2'
# tables.__version__ '2.0'

import TOPKAPI # import the TOPKAPI package

TOPKAPI.model.run('run_the_model/TOPKAPI.ini') # run the model using the configuration in TOPKAPI.ini

# plot the simulation results using the config in plot_Qsim_Qobs_Rain.ini
from TOPKAPI.results_analysis import plot_Qsim_Qobs_Rain

plot_Qsim_Qobs_Rain.run('analyse_the_results/plot_Qsim_Qobs_Rain.ini')

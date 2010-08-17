import pytopkapi

# run a model simulation using the configuration in TOPKAPI.ini
pytopkapi.run('run_the_model/TOPKAPI.ini')

# plot the simulation results (rainfall-runoff graphics) using the
# config in plot_Qsim_Qobs_Rain.ini
from pytopkapi.results_analysis import plot_Qsim_Qobs_Rain

plot_Qsim_Qobs_Rain.run('analyse_the_results/plot_Qsim_Qobs_Rain.ini')


# plot the simulation results (soil moisture maps) using the config in
# plot_soil_moisture_maps.ini
from pytopkapi.results_analysis import plot_soil_moisture_maps

plot_soil_moisture_maps.run('analyse_the_results/plot_soil_moisture_maps.ini')

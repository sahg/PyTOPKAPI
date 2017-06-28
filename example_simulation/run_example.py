import pytopkapi
from pytopkapi.results_analysis import plot_Qsim_Qobs_Rain, \
                                       plot_soil_moisture_maps

if __name__ == '__main__':
    # Run a model simulation using the configuration in
    # model-simulation.ini
    pytopkapi.run('model-simulation.ini')

    # Plot the simulation results (rainfall-runoff graphics) using the
    # config in plot-flow-precip.ini
    plot_Qsim_Qobs_Rain.run('plot-flow-precip.ini')


    # Plot the simulation results (soil moisture maps) using the config in
    # plot-soil-moisture-maps.ini
    plot_soil_moisture_maps.run('plot-soil-moisture-maps.ini')

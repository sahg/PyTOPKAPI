"""Core logic of PyTOPKAPI.

The `run` function in this module contains the logic to run a TOPKAPI
simulation based on the parameters specified in an INI file.

"""

#General module importation
import os.path
from configparser import SafeConfigParser

import h5py
import numpy as np
from tqdm import tqdm

#Personnal module importation
import pytopkapi
from . import utils as ut
from . import pretreatment as pm
from . import fluxes as fl
from . import ode as om
from . import evap as em
from .infiltration import green_ampt_cum_infiltration

no_data = np.nan

def run(ini_file='TOPKAPI.ini',
        verbose=False, quiet=False, parallel_exec=True):
    """Run the model.

    Parameters
    ----------
    ini_file : str
       The name of the PyTOPKAPI initialization file. This file
       describes the locations of the parameter files and model setup
       options. Default is to use a file named `TOPKAPI.ini` in the
       current directory.
    verbose : bool
        Prints runtime information [default False - don't display
        runtime info]. Is independent of the `quiet` keyword argument.
    quiet : bool
        Toggles whether to display an informational banner at runtime
        [default False - display banner]. Is independent of the
        `verbose` keyword argument.

    """

    ##================================##
    ##  Read the input file (*.ini)   ##
    ##================================##
    config = SafeConfigParser()
    config.read(ini_file)

    ##~~~~~~ Numerical_options ~~~~~~##
    solve_s = config.getint('numerical_options', 'solve_s')
    solve_o = config.getint('numerical_options', 'solve_o')
    solve_c = config.getint('numerical_options', 'solve_c')
    only_channel_output = config.getboolean('numerical_options',
                                            'only_channel_output')

    ##~~~~~~~~~~~ input files ~~~~~~~~~~~##
    #Param
    file_global_param = config.get('input_files', 'file_global_param')
    file_cell_param = config.get('input_files', 'file_cell_param')
    #Rain
    file_rain = config.get('input_files', 'file_rain')
    #ETP
    file_ET = config.get('input_files', 'file_ET')

    #~~~~~~~~~~~ Group (simulated event) ~~~~~~~~~~~##
    group_name = config.get('groups', 'group_name')

    ##~~~~~~ Calibration ~~~~~~##
    fac_L = config.getfloat('calib_params', 'fac_L')
    fac_Ks = config.getfloat('calib_params', 'fac_Ks')
    fac_n_o = config.getfloat('calib_params', 'fac_n_o')
    fac_n_c = config.getfloat('calib_params', 'fac_n_c')

    ##~~~~~~ External flows ~~~~~~##
    external_flow = config.getboolean('external_flow', 'external_flow')
    if external_flow:
        file_Qexternal_flow = config.get('external_flow',
                                         'file_Qexternal_flow')
        Xexternal_flow = config.getfloat('external_flow', 'Xexternal_flow')
        Yexternal_flow = config.getfloat('external_flow', 'Yexternal_flow')

    ##~~~~~~~~~~~ output files ~~~~~~~~~~##
    file_out = config.get('output_files', 'file_out')
    ut.check_file_exist(file_out) #create path_out if it doesn't exist
    if os.path.exists(file_out):
        first_run = False
    else:
        first_run = True

    append_output = config.getboolean('output_files', 'append_output')
    if append_output is True:
        fmode = 'a'
    else:
        fmode = 'w'

    ##============================##
    ##   Read the forcing data    ##
    ##============================##
    if verbose:
        print('Read the forcing data')

    #~~~~Rainfall
    h5_rain = h5py.File(file_rain)
    dset_name = '/{}/rainfall'.format(group_name)
    rainfall_forcing = h5_rain[dset_name][...]
    h5_rain.close()

    #~~~~ETr - Reference crop ET
    h5_ET = h5py.File(file_ET)
    dset_name = '/{}/ETr'.format(group_name)
    ETr_forcing = h5_ET[dset_name][...]

    #~~~~ETo - Open water potential evap.
    dset_name = '/{}/ETo'.format(group_name)
    ET0_forcing = h5_ET[dset_name][...]
    h5_ET.close()

    #~~~~external_flow flows
    if external_flow:
        external_flow_records = np.loadtxt(file_Qexternal_flow)[:, 5]
    else:
        external_flow_records = None


    ##============================##
    ## Pretreatment of input data ##
    ##============================##
    if verbose:
        print('Pretreatment of input data')

    #~~~~Read Global parameters file
    X, Dt, alpha_s, \
    alpha_o, alpha_c, \
    A_thres, W_min, W_max = pm.read_global_parameters(file_global_param)

    #~~~~Read Cell parameters file
    ar_cell_label, ar_coorx, \
    ar_coory, channel_flag, \
    Xc, ar_dam, \
    ar_tan_beta, ar_tan_beta_channel, \
    ar_L, Ks, \
    ar_theta_r, ar_theta_s, \
    ar_n_o, ar_n_c, \
    ar_cell_down, ar_pVs_t0, \
    ar_Vo_t0, ar_Qc_t0, \
    Kc, psi_b, lamda = pm.read_cell_parameters(file_cell_param)

    #~~~~Number of cell in the catchment
    nb_cell = len(ar_cell_label)

    #~~~~Computation of cell order
    node_hierarchy = pm.compute_node_hierarchy(ar_cell_label, ar_cell_down)
    ar_label_sort = pm.sort_cell(ar_cell_label, ar_cell_down)

    #~~~~Computation of upcells
    li_cell_up = pm.direct_up_cell(ar_cell_label, ar_cell_down, ar_label_sort)

    #~~~~Computation of drained area
    ar_A_drained = pm.drained_area(ar_label_sort, li_cell_up, X)

    #~~~~Apply calibration factors to the parameter values
    ar_L = ar_L*fac_L
    Ks = Ks*fac_Ks
    ar_n_o = ar_n_o*fac_n_o
    ar_n_c = ar_n_c*fac_n_c

    if verbose:
        print('Max L=', max(ar_L))
        print('Max Ks=', max(Ks))
        print('Max n_o=', max(ar_n_o))
        print('Max n_c=', max(ar_n_c))

    #~~~~Computation of model parameters from physical parameters
    Vsm, b_s, b_o, \
    W, b_c = pm.compute_cell_param(X, Xc, Dt, alpha_s,
                                         alpha_o, alpha_c, nb_cell,
                                         A_thres, W_max, W_min,
                                         channel_flag, ar_tan_beta,
                                         ar_tan_beta_channel, ar_L,
                                         Ks, ar_theta_r, ar_theta_s,
                                         ar_n_o, ar_n_c, ar_A_drained)

    #~~~~Look for the cell of external_flow tunnel
    if external_flow:
        cell_external_flow = ut.find_cell_coordinates(ar_cell_label,
                                                      Xexternal_flow,
                                                      Yexternal_flow,
                                                      ar_coorx,
                                                      ar_coory,
                                                      channel_flag)

        if verbose:
            print('external flows will be taken into account for cell no',\
                  cell_external_flow, ' coordinates ('\
                  ,Xexternal_flow,',',Yexternal_flow,')')
    else:
        cell_external_flow = None

    #~~~~Number of simulation time steps
    nb_time_step = rainfall_forcing.shape[0]


    ##=============================##
    ##  Variable array definition  ##
    ##=============================##

    ## Initialisation of the reservoirs
    #Matrix of soil,overland and channel store at the begining of the time step
    if append_output and not first_run:
        if verbose:
            print('Initialize from simulation file')

        h5file_in = h5py.File(file_out)

        Vs_t0 = h5file_in['/Soil/V_s'][-1, :]
        Vc_t0 = h5file_in['/Channel/V_c'][-1, :]
        Vo_t0 = h5file_in['/Overland/V_o'][-1, :]

        h5file_in.close()
    else:
        if verbose:
            print('Initialize from parameters')
        Vs_t0 = fl.initial_volume_soil(ar_pVs_t0, Vsm)
        Vo_t0 = ar_Vo_t0
        Vc_t0 = fl.initial_volume_channel(ar_Qc_t0, W, X, ar_n_c)

    ##=============================##
    ## HDF5 output file definition ##
    ##=============================##
    h5file, dset_Vs, dset_Vo, dset_Vc,     \
    dset_Qs_out, dset_Qo_out, dset_Qc_out, \
    dset_Q_down, dset_ET_out, dset_Ec_out  \
                                    = ut.open_simulation_file(file_out, fmode,
                                                   Vs_t0, Vo_t0, Vc_t0, no_data,
                                                   nb_cell, nb_time_step,
                                                   append_output, first_run)

    eff_theta = ar_theta_s - ar_theta_r

    ##===========================##
    ##     Core of the Model     ##
    ##===========================##
    if not quiet:
        ut.show_banner(ini_file, nb_cell, nb_time_step)
        progress_desc = 'Simulation'
    else:
        progress_desc = 'PyTOPKAPI v{}'.format(pytopkapi.__version__)

    # prepare parameter dict
    exec_params = {'nb_cell' : nb_cell,
                   'nb_time_step' : nb_time_step,
                   'progress_desc' : progress_desc,
                   'Dt' : Dt,
                   'rainfall_forcing' : rainfall_forcing,
                   'ETr_forcing' : ETr_forcing,
                   'ET0_forcing' : ET0_forcing,
                   'psi_b' : psi_b,
                   'lamda' : lamda,
                   'eff_theta' : eff_theta,
                   'Ks' : Ks,
                   'X' : X,
                   'b_s' : b_s,
                   'b_o' : b_o,
                   'b_c' : b_c,
                   'alpha_s' : alpha_s,
                   'alpha_o' : alpha_o,
                   'alpha_c' : alpha_c,
                   'Vs_t0' : Vs_t0,
                   'Vo_t0' : Vo_t0,
                   'Vc_t0' : Vc_t0,
                   'Vsm' : Vsm,
                   'dset_Vs' : dset_Vs,
                   'dset_Vo' : dset_Vo,
                   'dset_Vc' : dset_Vc,
                   'dset_Qs_out' : dset_Qs_out,
                   'dset_Qo_out' : dset_Qo_out,
                   'dset_Qc_out' : dset_Qc_out,
                   'dset_Q_down' : dset_Q_down,
                   'dset_ET_out' : dset_ET_out,
                   'dset_Ec_out' : dset_Ec_out,
                   'solve_s' : solve_s,
                   'solve_o' : solve_o,
                   'solve_c' : solve_c,
                   'channel_flag' : channel_flag,
                   'W' : W,
                   'Xc' : Xc,
                   'Kc' : Kc,
                   'cell_external_flow' : cell_external_flow,
                   'external_flow_records' : external_flow_records,
                   'node_hierarchy' : node_hierarchy,
                   'li_cell_up' : li_cell_up}

    if not parallel_exec:
        # Serial execution. Solve by timestep in a single process.
        # Outer loop timesteps - inner loop cells
        _serial_execute(exec_params)
    else:
        # Parallel execution. Solve by cell using multiple processes.
        # Outer loop cells - inner loop timesteps
        _parallel_execute(exec_params)

    h5file.close()

def _solve_cell(params):
    """Core calculations for a model cell.

    """
    Dt = params['Dt']
    rain_depth = params['rain_depth']
    psi = params['psi']
    eff_theta = params['eff_theta']
    eff_sat = params['eff_sat']
    Ks = params['Ks']
    X = params['X']
    soil_upstream_inflow = params['soil_upstream_inflow']
    b_s = params['b_s']
    alpha_s = params['alpha_s']
    Vs0 = params['Vs0']
    solve_s = params['solve_s']
    Vsm = params['Vsm']
    b_o = params['b_o']
    alpha_o = params['alpha_o']
    Vo0 = params['Vo0']
    solve_o = params['solve_o']
    channel_flag = params['channel_flag']
    W = params['W']
    Xc = params['Xc']
    channel_upstream_inflow = params['chan_up_inflow']
    Kc = params['Kc']
    ETr = params['ETr']
    b_c = params['b_c']
    alpha_c = params['alpha_c']
    Vc0 = params['Vc0']
    solve_c = params['solve_c']
    ET0 = params['ET0']
    external_flow_flag = params['external_flow_flag']
    external_flow = params['external_flow']

    ## ======================== ##
    ## ===== INFILTRATION ===== ##
    ## ======================== ##
    rain_rate = rain_depth/Dt

    infiltration_depth = green_ampt_cum_infiltration(rain_rate, psi,
                                                     eff_theta, eff_sat, Ks, Dt)

    ## ====================== ##
    ## ===== SOIL STORE ===== ##
    ## ====================== ##
    #~~~~ Computation of soil input
    a_s = fl.input_soil(infiltration_depth, Dt, X, soil_upstream_inflow)

    #~~~~ Resolution of the equation dV/dt=a_s-b_s*V^alpha_s
    # Calculate the volume in the soil store at the end of the
    # current time-step.
    Vs_prim = om.solve_storage_eq(a_s, b_s, alpha_s, Vs0, Dt, solve_s)

    #~~~~ Computation of soil outflow and overland input
    Qs_out, Vs1 = fl.output_soil(Vs0, Vs_prim, Vsm, a_s, b_s, alpha_s, Dt)

    ## ========================== ##
    ## ===== OVERLAND STORE ===== ##
    ## ========================== ##
    #~~~~ Computation of overland input
    rain_excess = rain_depth - infiltration_depth
    # convert mm to m^3/s
    rain_excess = max(0, (rain_excess*(10**-3)/Dt)*X**2)

    a_o = max(0, a_s - ((Vs1-Vs0)/Dt + Qs_out) + rain_excess)

    #~~~~ Resolution of the equation dV/dt=a_o-b_o*V^alpha_o

    Vo1 = om.solve_storage_eq(a_o, b_o, alpha_o, Vo0, Dt, solve_o)

    #~~~~ Computation of overland outflows
    Qo_out = fl.Qout_computing(Vo0, Vo1, a_o, Dt)

    ## ============================= ##
    ## ===== FLOW PARTITIONING ===== ##
    ## ============================= ##

    Q_down, Q_to_channel  = fl.flow_partitioning(channel_flag,
                                                 Qs_out, Qo_out, W, X, Xc)

    ## ======================== ##
    ## ===== CHANNEL STORE ==== ##
    ## ======================== ##
    if channel_flag == 1:
        #~~~~ Computation of channel input
        a_c = fl.input_channel(channel_upstream_inflow, Q_to_channel)

        if external_flow_flag:
            a_c = a_c + external_flow

        #~~~~ Resolution of the equation dV/dt=a_c-b_c*V^alpha_c

        Vc1 = om.solve_storage_eq(a_c, b_c, alpha_c, Vc0, Dt, solve_c)

        #~~~~ Computation of channel outflows
        Qc_out = fl.Qout_computing(Vc0, Vc1, a_c, Dt)

    else:
        a_c = 0.
        Vc1 = 0.
        Qc_out = 0.

    ## ============================== ##
    ## ===== EVAPOTRANSPIRATION ===== ##
    ## ============================== ##
    #~~~~~ From soil
    ETa, Vs1, Vo1 = em.evapot_soil_overland(Vo1, Vs1, Vsm, Kc, ETr, X)

    #~~~~~ Evaporation from channel
    if channel_flag == 1:
        ET_channel, Vc1 = em.evapor_channel(Vc1, ET0, W, Xc)
    else:
        ET_channel = 0

    return Qs_out, Qo_out, Qc_out, Q_down, Vs1, Vo1, Vc1, ETa, ET_channel

def _serial_execute(model_params):
    """Execute serially on a single processor core.

    """
    nb_cell = model_params['nb_cell']
    nb_time_step = model_params['nb_time_step']
    progress_desc = model_params['progress_desc']
    Vsm = model_params['Vsm']
    psi_b = model_params['psi_b']
    lamda = model_params['lamda']
    node_hierarchy = model_params['node_hierarchy']
    li_cell_up = model_params['li_cell_up']
    cell_external_flow = model_params['cell_external_flow']
    external_flow_records = model_params['external_flow_records']
    Dt = model_params['Dt']
    rainfall_forcing = model_params['rainfall_forcing']
    ETr_forcing = model_params['ETr_forcing']
    ET0_forcing = model_params['ET0_forcing']
    eff_theta = model_params['eff_theta']
    Ks = model_params['Ks']
    W = model_params['W']
    X = model_params['X']
    Xc = model_params['Xc']
    Kc = model_params['Kc']
    b_s = model_params['b_s']
    b_o = model_params['b_o']
    b_c = model_params['b_c']
    alpha_s = model_params['alpha_s']
    alpha_o = model_params['alpha_o']
    alpha_c = model_params['alpha_c']
    solve_s = model_params['solve_s']
    solve_o = model_params['solve_o']
    solve_c = model_params['solve_c']
    channel_flag = model_params['channel_flag']
    dset_Vs = model_params['dset_Vs']
    dset_Vo = model_params['dset_Vo']
    dset_Vc = model_params['dset_Vc']
    dset_Qs_out = model_params['dset_Qs_out']
    dset_Qo_out = model_params['dset_Qo_out']
    dset_Qc_out = model_params['dset_Qc_out']
    dset_Q_down = model_params['dset_Q_down']
    dset_ET_out = model_params['dset_ET_out']
    dset_Ec_out = model_params['dset_Ec_out']

    # Initialize start and end of timestep soil, overland and channel stores
    Vs0 = np.array(model_params['Vs_t0'])
    Vo0 = np.array(model_params['Vo_t0'])
    Vc0 = np.array(model_params['Vc_t0'])

    Vs1 = np.ones(nb_cell)*no_data
    Vo1 = np.ones(nb_cell)*no_data
    Vc1 = np.ones(nb_cell)*no_data

    # Initialize outflows during the time step
    Qs_out = np.ones(nb_cell)*no_data
    Qo_out = np.ones(nb_cell)*no_data
    Qc_out = np.zeros(nb_cell)

    Q_down = np.ones(nb_cell)*no_data

    # Initialize evaporation during the timestep
    ETa = np.zeros(nb_cell)
    ET_channel = np.zeros(nb_cell)

    ## Loop on time
    for t in tqdm(range(nb_time_step), ascii=True, desc=progress_desc):
        eff_sat = Vs0/Vsm

        # estimate soil suction head using Brookes and Corey (1964)
        psi = psi_b/np.power(eff_sat, 1.0/lamda)

        ## Loop on cell hierarchy
        for lvl in range(len(node_hierarchy.keys())):
            for cell in node_hierarchy[lvl]:

                soil_upstream_inflow = Q_down[li_cell_up[cell]]
                channel_upstream_inflow = Qc_out[li_cell_up[cell]]

                if cell == cell_external_flow:
                    cell_params = {'Dt' : Dt,
                                  'rain_depth' : rainfall_forcing[t, cell],
                                  'psi' : psi[cell],
                                  'eff_theta' : eff_theta[cell],
                                  'eff_sat' : eff_sat[cell],
                                  'Ks' : Ks[cell],
                                  'X' : X,
                                  'soil_upstream_inflow' : soil_upstream_inflow,
                                  'b_s' : b_s[cell],
                                  'alpha_s' : alpha_s,
                                  'Vs0' : Vs0[cell],
                                  'solve_s' : solve_s,
                                  'Vsm' : Vsm[cell],
                                  'b_o' : b_o[cell],
                                  'alpha_o' : alpha_o,
                                  'Vo0' : Vo0[cell],
                                  'solve_o' : solve_o,
                                  'channel_flag' : channel_flag[cell],
                                  'W' : W[cell],
                                  'Xc' : Xc[cell],
                                  'chan_up_inflow' : channel_upstream_inflow,
                                  'Kc' : Kc[cell],
                                  'ETr' : ETr_forcing[t, cell],
                                  'b_c' : b_c[cell],
                                  'alpha_c' : alpha_c,
                                  'Vc0' : Vc0[cell],
                                  'solve_c' : solve_c,
                                  'ET0' : ET0_forcing[t, cell],
                                  'external_flow_flag' : True,
                                  'external_flow' : external_flow_records[t]}
                else:
                    cell_params = {'Dt' : Dt,
                                  'rain_depth' : rainfall_forcing[t, cell],
                                  'psi' : psi[cell],
                                  'eff_theta' : eff_theta[cell],
                                  'eff_sat' : eff_sat[cell],
                                  'Ks' : Ks[cell],
                                  'X' : X,
                                  'soil_upstream_inflow' : soil_upstream_inflow,
                                  'b_s' : b_s[cell],
                                  'alpha_s' : alpha_s,
                                  'Vs0' : Vs0[cell],
                                  'solve_s' : solve_s,
                                  'Vsm' : Vsm[cell],
                                  'b_o' : b_o[cell],
                                  'alpha_o' : alpha_o,
                                  'Vo0' : Vo0[cell],
                                  'solve_o' : solve_o,
                                  'channel_flag' : channel_flag[cell],
                                  'W' : W[cell],
                                  'Xc' : Xc[cell],
                                  'chan_up_inflow' : channel_upstream_inflow,
                                  'Kc' : Kc[cell],
                                  'ETr' : ETr_forcing[t, cell],
                                  'b_c' : b_c[cell],
                                  'alpha_c' : alpha_c,
                                  'Vc0' : Vc0[cell],
                                  'solve_c' : solve_c,
                                  'ET0' : ET0_forcing[t, cell],
                                  'external_flow_flag' : False,
                                  'external_flow' : None}

                Qs_out[cell], Qo_out[cell], Qc_out[cell], Q_down[cell], \
                Vs1[cell], Vo1[cell], Vc1[cell], \
                ETa[cell], ET_channel[cell] = \
                                           _solve_cell(cell_params)

        ####===================================####
        #### Affectation of new vector values  ####
        ####===================================####
        Vs0 = np.array(Vs1)
        Vo0 = np.array(Vo1)
        Vc0 = np.array(Vc1)

        ####===================================####
        #### Results writing at each time step ####
        ####===================================####
        dset_Vs[t+1] = Vs1
        dset_Vo[t+1] = Vo1
        dset_Vc[t+1] = Vc1

        dset_Qs_out[t+1] = Qs_out
        dset_Qo_out[t+1] = Qo_out
        dset_Qc_out[t+1] = Qc_out

        dset_Q_down[t+1] = Q_down

        dset_ET_out[t+1] = ETa
        dset_Ec_out[t+1] = ET_channel*1e-3 * W * Xc

def _parallel_execute(model_params):
    nb_cell = model_params['nb_cell']
    nb_time_step = model_params['nb_time_step']
    progress_desc = model_params['progress_desc']

    node_hierarchy = model_params['node_hierarchy']
    li_cell_up = model_params['li_cell_up']

    with tqdm(total=nb_cell, ascii=True, desc=progress_desc) as pbar:
        ## Loop on cell hierarchy
        for lvl in range(len(node_hierarchy.keys())):
            for cell in node_hierarchy[lvl]:

                if cell == model_params['cell_external_flow']:
                    external_flow_flag = True
                else:
                    external_flow_flag = False

                if len(li_cell_up[cell]) > 0:
                    soil_upstream_inflow =    \
                               model_params['dset_Q_down'][1:, li_cell_up[cell]]
                    channel_upstream_inflow = \
                               model_params['dset_Qc_out'][1:, li_cell_up[cell]]
                else:
                    soil_upstream_inflow = [np.array([])
                                            for i in range(nb_time_step)]
                    channel_upstream_inflow = [np.array([])
                                               for i in range(nb_time_step)]

                ts_params = {
                 'nb_time_step' : model_params['nb_time_step'],
                 'Vs_t0' : model_params['Vs_t0'][cell],
                 'Vo_t0' : model_params['Vo_t0'][cell],
                 'Vc_t0' : model_params['Vc_t0'][cell],
                 'psi_b' : model_params['psi_b'][cell],
                 'lamda' : model_params['lamda'][cell],
                 'external_flow_flag' : external_flow_flag,
                 'rainfall_forcing' : model_params['rainfall_forcing'][:, cell],
                 'ETr_forcing' : model_params['ETr_forcing'][:, cell],
                 'ET0_forcing' : model_params['ET0_forcing'][:, cell],
                 'soil_upstream_inflow' : soil_upstream_inflow,
                 'channel_upstream_inflow' : channel_upstream_inflow,
                 'eff_theta' : model_params['eff_theta'][cell],
                 'X' : model_params['X'],
                 'W' : model_params['W'][cell],
                 'Dt' : model_params['Dt'],
                 'Xc' : model_params['Xc'][cell],
                 'Kc' : model_params['Kc'][cell],
                 'Ks' : model_params['Ks'][cell],
                 'b_s' : model_params['b_s'][cell],
                 'b_o' : model_params['b_o'][cell],
                 'b_c' : model_params['b_c'][cell],
                 'alpha_s' : model_params['alpha_s'],
                 'alpha_o' : model_params['alpha_o'],
                 'alpha_c' : model_params['alpha_c'],
                 'solve_s' : model_params['solve_s'],
                 'solve_o' : model_params['solve_o'],
                 'solve_c' : model_params['solve_c'],
                 'Vsm' : model_params['Vsm'][cell],
                 'channel_flag' : model_params['channel_flag'][cell],
                 'external_flow_records' : model_params['external_flow_records']
                    }

                ts_result = _solve_cell_timeseries(ts_params)

                # Write results to disk
                model_params['dset_Vs'][1:, cell] = ts_result['Vs1']
                model_params['dset_Vo'][1:, cell] = ts_result['Vo1']
                model_params['dset_Vc'][1:, cell] = ts_result['Vc1']

                model_params['dset_Qs_out'][1:, cell] = ts_result['Qs_out']
                model_params['dset_Qo_out'][1:, cell] = ts_result['Qo_out']
                model_params['dset_Qc_out'][1:, cell] = ts_result['Qc_out']

                model_params['dset_Q_down'][1:, cell] = ts_result['Q_down']

                model_params['dset_ET_out'][1:, cell] = ts_result['ETa']
                Ec = ts_result['ET_channel']*1e-3
                Ec = Ec * model_params['W'][cell]
                Ec = Ec * model_params['Xc'][cell]
                model_params['dset_Ec_out'][1:, cell] = Ec

                pbar.update()

def _solve_cell_timeseries(tseries_params):
    """Solve the full simulation timeseries for a single cell.

    """
    # Initial values for cell store volumes
    Vs0 = tseries_params['Vs_t0']
    Vo0 = tseries_params['Vo_t0']
    Vc0 = tseries_params['Vc_t0']

    # Storage for end of timestep soil, overland and channel store volumes
    Vs1 = np.ones(tseries_params['nb_time_step'])*no_data
    Vo1 = np.ones(tseries_params['nb_time_step'])*no_data
    Vc1 = np.ones(tseries_params['nb_time_step'])*no_data

    # Storage for outflows during the time step
    Qs_out = np.ones(tseries_params['nb_time_step'])*no_data
    Qo_out = np.ones(tseries_params['nb_time_step'])*no_data
    Qc_out = np.zeros(tseries_params['nb_time_step'])

    Q_down = np.ones(tseries_params['nb_time_step'])*no_data

    # Storage for evaporation during the timestep
    ETa = np.zeros(tseries_params['nb_time_step'])
    ET_channel = np.zeros(tseries_params['nb_time_step'])

    ## Loop on time
    for t in range(tseries_params['nb_time_step']):
        eff_sat = Vs0/tseries_params['Vsm']

        # estimate soil suction head using Brookes and Corey (1964)
        psi = tseries_params['psi_b']/np.power(eff_sat,
                                               1.0/tseries_params['lamda'])

        if tseries_params['external_flow_flag']:
            # This cell receives external flow inputs
            cell_params = {
             'rain_depth' : tseries_params['rainfall_forcing'][t],
             'ETr' : tseries_params['ETr_forcing'][t],
             'ET0' : tseries_params['ET0_forcing'][t],
             'soil_upstream_inflow' : tseries_params['soil_upstream_inflow'][t],
             'chan_up_inflow' : tseries_params['channel_upstream_inflow'][t],
             'psi' : psi,
             'eff_sat' : eff_sat,
             'eff_theta' : tseries_params['eff_theta'],
             'X' : tseries_params['X'],
             'W' : tseries_params['W'],
             'Dt' : tseries_params['Dt'],
             'Xc' : tseries_params['Xc'],
             'Kc' : tseries_params['Kc'],
             'Ks' : tseries_params['Ks'],
             'b_s' : tseries_params['b_s'],
             'b_o' : tseries_params['b_o'],
             'b_c' : tseries_params['b_c'],
             'alpha_s' : tseries_params['alpha_s'],
             'alpha_o' : tseries_params['alpha_o'],
             'alpha_c' : tseries_params['alpha_c'],
             'solve_s' : tseries_params['solve_s'],
             'solve_o' : tseries_params['solve_o'],
             'solve_c' : tseries_params['solve_c'],
             'Vs0' : Vs0,
             'Vo0' : Vo0,
             'Vc0' : Vc0,
             'Vsm' : tseries_params['Vsm'],
             'channel_flag' : tseries_params['channel_flag'],
             'external_flow_flag' : tseries_params['external_flow_flag'],
             'external_flow' : tseries_params['external_flow_records'][t]
                            }
        else:
            cell_params = {
             'rain_depth' : tseries_params['rainfall_forcing'][t],
             'ETr' : tseries_params['ETr_forcing'][t],
             'ET0' : tseries_params['ET0_forcing'][t],
             'soil_upstream_inflow' : tseries_params['soil_upstream_inflow'][t],
             'chan_up_inflow' : tseries_params['channel_upstream_inflow'][t],
             'psi' : psi,
             'eff_sat' : eff_sat,
             'eff_theta' : tseries_params['eff_theta'],
             'X' : tseries_params['X'],
             'W' : tseries_params['W'],
             'Dt' : tseries_params['Dt'],
             'Xc' : tseries_params['Xc'],
             'Kc' : tseries_params['Kc'],
             'Ks' : tseries_params['Ks'],
             'b_s' : tseries_params['b_s'],
             'b_o' : tseries_params['b_o'],
             'b_c' : tseries_params['b_c'],
             'alpha_s' : tseries_params['alpha_s'],
             'alpha_o' : tseries_params['alpha_o'],
             'alpha_c' : tseries_params['alpha_c'],
             'solve_s' : tseries_params['solve_s'],
             'solve_o' : tseries_params['solve_o'],
             'solve_c' : tseries_params['solve_c'],
             'Vs0' : Vs0,
             'Vo0' : Vo0,
             'Vc0' : Vc0,
             'Vsm' : tseries_params['Vsm'],
             'channel_flag' : tseries_params['channel_flag'],
             'external_flow_flag' : tseries_params['external_flow_flag'],
             'external_flow' : None
                            }

        Qs_out[t], Qo_out[t], Qc_out[t], Q_down[t], \
        Vs1[t], Vo1[t], Vc1[t], ETa[t], ET_channel[t] = _solve_cell(cell_params)

        # Update storages for next time-step
        Vs0 = Vs1[t]
        Vo0 = Vo1[t]
        Vc0 = Vc1[t]

    result = {
        'Vs1' : Vs1,
        'Vo1' : Vo1,
        'Vc1' : Vc1,
        'Qs_out' : Qs_out,
        'Qo_out' : Qo_out,
        'Qc_out' : Qc_out,
        'Q_down' : Q_down,
        'ETa' : ETa,
        'ET_channel' : ET_channel
              }

    return result

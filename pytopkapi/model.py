"""Core logic of PyTOPKAPI.

The `run` function in this module contains the logic to run a TOPKAPI
simulation based on the parameters specified in an INI file.

"""

#General module importation
import os.path
from configparser import SafeConfigParser

import h5py
import numpy as np
import tables as h5

#Personnal module importation
import pytopkapi
from . import utils as ut
from . import pretreatment as pm
from . import fluxes as fl
from . import ode as om
from . import evap as em
from .infiltration import green_ampt_cum_infiltration

def run(ini_file='TOPKAPI.ini'):
    """Run the model with the set-up defined by `ini_file`.

    """

    ##================================##
    ##  Read the input file (*.ini)   ##
    ##================================##
    config = SafeConfigParser()
    config.read(ini_file)
    print('Read the file ',ini_file)

    ##~~~~~~ Numerical_options ~~~~~~##
    solve_s = config.getfloat('numerical_options', 'solve_s')
    solve_o = config.getfloat('numerical_options', 'solve_o')
    solve_c = config.getfloat('numerical_options', 'solve_c')
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
    print('Read the forcing data')

    #~~~~Rainfall
    h5file_in = h5.open_file(file_rain,mode='r')
    group = '/'+group_name+'/'
    node = h5file_in.get_node(group+'rainfall')
    rainfall_forcing = node.read()
    h5file_in.close()

    #~~~~ETr - Reference crop ET
    h5file_in = h5.open_file(file_ET,mode='r')
    group = '/'+group_name+'/'
    node = h5file_in.get_node(group+'ETr')
    ETr_forcing = node.read()
    h5file_in.close()

    #~~~~ETo - Open water potential evap.
    h5file_in = h5.open_file(file_ET,mode='r')
    group = '/'+group_name+'/'
    node = h5file_in.get_node(group+'ETo')
    ET0_forcing = node.read()
    h5file_in.close()

    #~~~~external_flow flows
    if external_flow:
        external_flow_records = np.loadtxt(file_Qexternal_flow)[:, 5]


    ##============================##
    ## Pretreatment of input data ##
    ##============================##
    print('Pretreatment of input data')

    #~~~~Read Global parameters file
    X, Dt, alpha_s, \
    alpha_o, alpha_c, \
    A_thres, W_min, W_max = pm.read_global_parameters(file_global_param)

    #~~~~Read Cell parameters file
    ar_cell_label, ar_coorx, \
    ar_coory, ar_lambda, \
    Xc, ar_dam, \
    ar_tan_beta, ar_tan_beta_channel, \
    ar_L0, Ks0, \
    ar_theta_r, ar_theta_s, \
    ar_n_o0, ar_n_c0, \
    ar_cell_down, ar_pVs_t0, \
    ar_Vo_t0, ar_Qc_t0, \
    kc, psi_b, lamda = pm.read_cell_parameters(file_cell_param)

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
    ar_L = ar_L0*fac_L
    Ks = Ks0*fac_Ks
    ar_n_o = ar_n_o0*fac_n_o
    ar_n_c = ar_n_c0*fac_n_c

    print('Max L=', max(ar_L))
    print('Max Ks=', max(Ks))
    print('Max n_o=', max(ar_n_o))
    print('Max n_c=', max(ar_n_c))

    #~~~~Computation of model parameters from physical parameters
    Vsm, b_s, b_o, \
    W, b_c = pm.compute_cell_param(X, Xc, Dt, alpha_s,
                                         alpha_o, alpha_c, nb_cell,
                                         A_thres, W_max, W_min,
                                         ar_lambda, ar_tan_beta,
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
                                                      ar_lambda)

        print('external flows will be taken into account for cell no',\
            cell_external_flow, ' coordinates ('\
            ,Xexternal_flow,',',Yexternal_flow,')')

    #~~~~Number of simulation time steps
    nb_time_step = len(rainfall_forcing[:,0])


    ##=============================##
    ##  Variable array definition  ##
    ##=============================##

    ## Initialisation of the reservoirs
    #Matrix of soil,overland and channel store at the begining of the time step
    if append_output and not first_run:
        print('Initialize from file')

        h5file_in = h5py.File(file_out)

        Vs0 = h5file_in['/Soil/V_s'][-1, :]
        Vc0 = h5file_in['/Channel/V_c'][-1, :]
        Vo0 = h5file_in['/Overland/V_o'][-1, :]

        h5file_in.close()
    else:
        print('Initialize from parms')
        Vs0 = fl.initial_volume_soil(ar_pVs_t0, Vsm)
        Vo0 = ar_Vo_t0
        Vc0 = fl.initial_volume_channel(ar_Qc_t0, W, X, ar_n_c)

    ## Computed variables
    #Matrix of soil,overland and channel store at the end of the time step
    ar_Vs1 = np.ones(nb_cell)*-99.9
    ar_Vo1 = np.ones(nb_cell)*-99.9
    ar_Vc1 = np.ones(nb_cell)*-99.9

    #Matrix of outflows between two time steps
    ar_Qs_out = np.ones(nb_cell)*-99.9
    ar_Qo_out = np.ones(nb_cell)*-99.9
    ar_Qc_out = np.zeros(nb_cell)

    ## Intermediate variables
    ar_Q_to_next_cell = np.ones(nb_cell)*-99.9
    ar_Q_to_channel = np.ones(nb_cell)*-99.9
    ar_Qc_cell_up = np.zeros(nb_cell)
    ar_ETa = np.zeros(nb_cell)
    ar_ET_channel = np.zeros(nb_cell)


    ##=============================##
    ## HDF5 output file definition ##
    ##=============================##
    h5file = h5.open_file(file_out, mode=fmode, title='TOPKAPI_out')

    root = h5file.get_node('/')
    root._v_attrs.pytopkapi_version = pytopkapi.__version__
    root._v_attrs.pytopkapi_git_revision = pytopkapi.__git_revision__

    atom = h5.Float32Atom()
    h5filter = h5.Filters(9)# maximum compression

    # create file structure as necessary
    grp_name = '/Soil'
    if grp_name not in h5file:
        h5file.create_group('/', 'Soil', 'Soil arrays')
    if grp_name+'/Qs_out' not in h5file:
        array_Qs_out = h5file.create_earray(grp_name, 'Qs_out',
                                           atom, shape=(0,nb_cell),
                                           title='m3/s', filters=h5filter,
                                           expectedrows=nb_time_step)
    else:
        array_Qs_out = h5file.get_node(grp_name+'/Qs_out')
    if grp_name+'/V_s' not in h5file:
        array_Vs = h5file.create_earray(grp_name, 'V_s',
                                       atom, shape=(0, nb_cell),
                                       title='m3', filters=h5filter,
                                       expectedrows=nb_time_step+1)
    else:
        array_Vs = h5file.get_node(grp_name+'/V_s')

    grp_name = '/Overland'
    if grp_name not in h5file:
        h5file.create_group('/', 'Overland', 'Overland arrays')
    if grp_name+'/Qo_out' not in h5file:
        array_Qo_out = h5file.create_earray(grp_name, 'Qo_out',
                                           atom, shape=(0,nb_cell),
                                           title='m3/s', filters=h5filter,
                                           expectedrows=nb_time_step)
    else:
        array_Qo_out = h5file.get_node(grp_name+'/Qo_out')
    if grp_name+'/V_o' not in h5file:
        array_Vo = h5file.create_earray(grp_name, 'V_o',
                                       atom, shape=(0,nb_cell),
                                       title='m3', filters=h5filter,
                                       expectedrows=nb_time_step+1)
    else:
        array_Vo = h5file.get_node(grp_name+'/V_o')

    grp_name = '/Channel'
    if grp_name not in h5file:
        h5file.create_group('/', 'Channel', 'Channel arrays')
    if grp_name+'/Qc_out' not in h5file:
        array_Qc_out = h5file.create_earray(grp_name, 'Qc_out',
                                           atom, shape=(0,nb_cell),
                                           title='m3/s', filters=h5filter,
                                           expectedrows=nb_time_step)
    else:
        array_Qc_out = h5file.get_node(grp_name+'/Qc_out')
    if grp_name+'/V_c' not in h5file:
        array_Vc = h5file.create_earray(grp_name, 'V_c',
                                       atom, shape=(0,nb_cell),
                                       title='m3', filters=h5filter,
                                       expectedrows=nb_time_step)
    else:
        array_Vc = h5file.get_node(grp_name+'/V_c')
    if grp_name+'/Ec_out' not in h5file:
        array_Ec_out = h5file.create_earray(grp_name, 'Ec_out',
                                           atom, shape=(0,nb_cell),
                                           title='m3', filters=h5filter,
                                           expectedrows=nb_time_step)
    else:
        array_Ec_out = h5file.get_node(grp_name+'/Ec_out')

    if '/ET_out' not in h5file:
        array_ET_out = h5file.create_earray('/', 'ET_out',
                                           atom, shape=(0,nb_cell),
                                           title='mm', filters=h5filter,
                                           expectedrows=nb_time_step)
    else:
        array_ET_out = h5file.get_node('/ET_out')

    if '/Q_down' not in h5file:
        array_Q_down = h5file.create_earray('/', 'Q_down',
                                           atom, shape=(0,nb_cell),
                                           title='m3/s', filters=h5filter,
                                           expectedrows=nb_time_step)
    else:
        array_Q_down = h5file.get_node('/Q_down')

    if append_output is False or first_run is True:
        #Write the initial values into the output file
        array_Vs.append(Vs0.reshape((1,nb_cell)))
        array_Vo.append(Vo0.reshape((1,nb_cell)))
        array_Vc.append(Vc0.reshape((1,nb_cell)))

        array_Qs_out.append(ar_Qs_out.reshape((1,nb_cell)))
        array_Qo_out.append(ar_Qo_out.reshape((1,nb_cell)))
        array_Qc_out.append(ar_Qc_out.reshape((1,nb_cell)))

        array_Q_down.append(ar_Q_to_next_cell.reshape((1,nb_cell)))

        array_ET_out.append(ar_ETa.reshape((1,nb_cell)))

        E_vol = ar_ET_channel*1e-3 * W * Xc
        array_Ec_out.append(E_vol.reshape((1,nb_cell)))

    eff_theta = ar_theta_s - ar_theta_r

    ##===========================##
    ##     Core of the Model     ##
    ##===========================##
    print('** NB_CELL=',nb_cell)
    print('** NB_TIME_STEP=',nb_time_step)
    print('--> SIMULATIONS <--')

    ## Loop on time
    for t in range(nb_time_step):
        print(t+1, '/', nb_time_step)

        eff_sat = Vs0/Vsm

        # estimate soil suction head using Brookes and Corey (1964)
        psi = psi_b/np.power(eff_sat, 1.0/lamda)

        ## Loop on cell hierarchy
        for lvl in range(len(node_hierarchy.keys())):
            for cell1 in node_hierarchy[lvl]:
                cell=np.where(ar_cell_label==cell1)[0][0]


                if external_flow:
                    _solve_cell(cell,
                                Dt, rainfall_forcing[t, cell], psi[cell],
                                eff_theta[cell], eff_sat[cell],Ks[cell], X,
                                ar_Q_to_next_cell, li_cell_up, b_s[cell],
                                alpha_s, Vs0[cell], solve_s, Vsm[cell],
                                ar_Qs_out, ar_Vs1, b_o[cell], alpha_o,
                                Vo0[cell], solve_o, ar_Vo1, ar_Qo_out,
                                ar_lambda, W[cell], Xc[cell], ar_Q_to_channel,
                                ar_Qc_out, ar_Qc_cell_up, ar_cell_label,
                                ar_Vc1, kc[cell], ETr_forcing[t, cell], ar_ETa,
                                ar_cell_down, b_c[cell], alpha_c, Vc0[cell],
                                solve_c, ET0_forcing[t, cell], ar_ET_channel,
                                external_flow,
                                cell_external_flow, external_flow_records[t])
                else:
                    _solve_cell(cell,
                                Dt, rainfall_forcing[t, cell], psi[cell],
                                eff_theta[cell], eff_sat[cell],Ks[cell], X,
                                ar_Q_to_next_cell, li_cell_up, b_s[cell],
                                alpha_s, Vs0[cell], solve_s, Vsm[cell],
                                ar_Qs_out, ar_Vs1, b_o[cell], alpha_o,
                                Vo0[cell], solve_o, ar_Vo1, ar_Qo_out,
                                ar_lambda, W[cell], Xc[cell], ar_Q_to_channel,
                                ar_Qc_out, ar_Qc_cell_up, ar_cell_label,
                                ar_Vc1, kc[cell], ETr_forcing[t, cell], ar_ETa,
                                ar_cell_down, b_c[cell], alpha_c, Vc0[cell],
                                solve_c, ET0_forcing[t, cell], ar_ET_channel,
                                external_flow)

        ####===================================####
        #### Affectation of new vector values  ####
        ####===================================####
        Vs0 = np.array(ar_Vs1)
        Vo0 = np.array(ar_Vo1)
        Vc0 = np.array(ar_Vc1)

        ####===================================####
        #### Results writing at each time step ####
        ####===================================####
        array_Vs.append(ar_Vs1.reshape((1,nb_cell)))
        array_Vo.append(ar_Vo1.reshape((1,nb_cell)))
        array_Vc.append(ar_Vc1.reshape((1,nb_cell)))

        array_Qs_out.append(ar_Qs_out.reshape((1,nb_cell)))
        array_Qo_out.append(ar_Qo_out.reshape((1,nb_cell)))
        array_Qc_out.append(ar_Qc_out.reshape((1,nb_cell)))

        array_Q_down.append(ar_Q_to_next_cell.reshape((1,nb_cell)))

        array_ET_out.append(ar_ETa.reshape((1,nb_cell)))

        E_vol = ar_ET_channel*1e-3 * W * Xc
        array_Ec_out.append(E_vol.reshape((1,nb_cell)))

    h5file.close()

    print(' ')
    print('***** THE END *****')

def _solve_cell(cell,
                Dt, rain_depth, psi, eff_theta, eff_sat, Ks, X,
                ar_Q_to_next_cell, li_cell_up, b_s, alpha_s, Vs0,
                solve_s, Vsm, ar_Qs_out, ar_Vs1, b_o, alpha_o,
                Vo0, solve_o, ar_Vo1, ar_Qo_out, ar_lambda, W, Xc,
                ar_Q_to_channel, ar_Qc_out,
                ar_Qc_cell_up, ar_cell_label, ar_Vc1, kc, ETr, ar_ETa,
                ar_cell_down,b_c, alpha_c, Vc0, solve_c, ET0,
                ar_ET_channel, external_flow_flag, cell_external_flow=None,
                external_flow=None):
    """Core calculations for a model cell.

    """
    ## ======================== ##
    ## ===== INTERCEPTION ===== ##
    ## ======================== ##
    ## No interception for the moment

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
    a_s = fl.input_soil(infiltration_depth, Dt,
                                 X, ar_Q_to_next_cell, li_cell_up[cell])

    #~~~~ Resolution of the equation dV/dt=a_s-b_s*V^alpha_s
    # Calculate the volume in the soil store at the end of the
    # current time-step.

    Vs_prim = om.solve_storage_eq(a_s, b_s,
                                  alpha_s, Vs0, Dt, solve_s)

    #~~~~ Computation of soil outflow and overland input
    ar_Qs_out[cell], ar_Vs1[cell] = fl.output_soil(Vs0, Vs_prim,
                                                   Vsm, a_s,
                                                   b_s, alpha_s, Dt)

    if ar_Qs_out[cell] < 0:
        print('Problem Soil:output greater than input....')
        print('n=', n, 'label=', cell)
        stop

    ## ========================== ##
    ## ===== OVERLAND STORE ===== ##
    ## ========================== ##
    #~~~~ Computation of overland input
    rain_excess = rain_depth - infiltration_depth
    # convert mm to m^3/s
    rain_excess = max(0, (rain_excess*(10**-3)/Dt)*X**2)

    a_o = max(0,
                       a_s \
                       - ((ar_Vs1[cell]-Vs0)/Dt \
                       + ar_Qs_out[cell]) \
                       + rain_excess)

    #~~~~ Resolution of the equation dV/dt=a_o-b_o*V^alpha_o

    ar_Vo1[cell] = om.solve_storage_eq(a_o, b_o,
                                       alpha_o, Vo0, Dt, solve_o)

    #~~~~ Computation of overland outflows
    ar_Qo_out[cell] = fl.Qout_computing(Vo0,
                                        ar_Vo1[cell], a_o, Dt)

    if ar_Qo_out[cell] < 0:
        print('Problem Overland:output greater than input....')
        print('n=', n, 'label=', cell)
        stop

    ## ============================= ##
    ## ===== FLOW PARTITIONING ===== ##
    ## ============================= ##

    ar_Q_to_next_cell[cell], \
    ar_Q_to_channel[cell]  = fl.flow_partitioning(ar_lambda[cell],
                                                     ar_Qs_out[cell],
                                                     ar_Qo_out[cell],
                                                     W, X, Xc)

    ## ======================== ##
    ## ===== CHANNEL STORE ==== ##
    ## ======================== ##
    if ar_lambda[cell] == 1:
        if ar_cell_down[cell] >= 0 \
        and ar_lambda[ar_cell_down[cell]] == 0:

            print('Problem: the present cell has a channel but not the cell down...')
            Stop

        #~~~~ Computation of channel input
        a_c, \
        ar_Qc_cell_up[cell] = fl.input_channel(ar_Qc_out,
                                               ar_Q_to_channel[cell],
                                               li_cell_up[cell])

        #TO DO: Handle external flows properly. Vars not passed into
        #this function currently.
        if external_flow_flag \
          and cell == np.where(ar_cell_label==cell_external_flow)[0][0]:
            a_c = a_c + external_flow

        #~~~~ Resolution of the equation dV/dt=a_c-b_c*V^alpha_c

        ar_Vc1[cell] = om.solve_storage_eq(a_c, b_c,
                                           alpha_c, Vc0, Dt, solve_c)

        #~~~~ Computation of channel outflows
        ar_Qc_out[cell] = fl.Qout_computing(Vc0,
                                            ar_Vc1[cell], a_c, Dt)

        if ar_Qc_out[cell] < 0:
            print('Problem Channel: output greater than input....')
            stop
        if str(ar_Qc_out[cell]).count('N') > 0:
            print(ar_Qc_out[cell])
            print('Problem Channel: Non authorized operand....')
            stop

    else:
        a_c = 0.
        ar_Vc1[cell] = 0.
        ar_Qc_out[cell] = 0.


    ## ============================== ##
    ## ===== EVAPOTRANSPIRATION ===== ##
    ## ============================== ##
    #~~~~~ From soil
    ar_ETa[cell], \
    ar_Vs1[cell], \
    ar_Vo1[cell] = em.evapot_soil_overland(ar_Vo1[cell], ar_Vs1[cell],
                                           Vsm, kc,
                                           ETr, X)

    #~~~~~ Evaporation from channel
    if ar_lambda[cell] == 1:
        ar_ET_channel[cell], \
        ar_Vc1[cell] = em.evapor_channel(ar_Vc1[cell],
                                         ET0,
                                         W, Xc)

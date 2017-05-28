"""

Functions required to compute the intrinsic TOPKAPI parameters from the
physical parameters

"""

import numpy as np
import networkx as nx

def read_global_parameters(file_name):
    """Read global model parameters from file.

    Read the file that specifies the parameters of the model, which
    are globally applied to all cells in the model.

    Parameters
    ----------
    file_name : string
        Full name of the parameter file (incl. path)

    Returns
    -------
    X : scalar
        The lateral dimension of the grid-cell (:math:`m`)
    Dt : scalar
        The length of the current time-step in seconds
    alpha_s : scalar
        The dimensionless pore-size distribution parameter for the
        soil store
    alpha_o : scalar
        Power co-efficient for Mannings Equation applied to the
        overland store the value is typically 5/3
    alpha_c : scalar
        Power co-efficient for Mannings Equation applied to the
        channel store the value is typically 5/3
    A_thres : scalar
        The minimum area of upstream contributing cells required
        before a cell is considered to initiate a river channel
        (:math:`m^2`)
    W_min : scalar
        The minimum width of a channel (:math:`m`)
    W_max : scalar
        The maximum width of a channel (:math:`m`)

    """
    file_read=open(file_name,'r')
    ## reading lines
    tab_read=file_read.readlines()
    for line in tab_read[1:]:
        ##split the blanks
        donnees=line.split()
        ##append the values to variable lists
        X=float(donnees[0])
        Dt=float(donnees[1])
        alpha_s=float(donnees[2])
        alpha_o=float(donnees[3])
        alpha_c=float(donnees[4])
        A_thres=float(donnees[5])
        W_min=float(donnees[6])
        W_max=float(donnees[7])
        ##end of iteration the file is closed
        file_read.close()

    return X,Dt,alpha_s,alpha_o,alpha_c,A_thres,W_min,W_max

def read_cell_parameters(file_name):
    """Read the spatially variable cell parameters from file.

    Read the file containing the physical parameters of each
    cell. This information governs the distributed behaviour of the
    model.

    Parameters
    ----------
    file_name : string
        Full name of the parameter file (incl. path)

    Returns
    -------
    ar_cell_label : (N,) int array
        Numbers labelling each cell
    ar_coorx : (N,) float array
        The x co-ordinate of the centre of each cell (:math:`m`). This
        is the Longitude expressed in metres using an appropriate map
        projection
    ar_coory : (N,) float array
        The y co-ordinate of the centre of each cell (:math:`m`). This
        is the Latitude expressed in metres using an appropriate map
        projection
    ar_lambda : (N,) int array
        Switch indicating whether the current cell contains a
        channel. A value of `1` indicates a channel cell, `0`
        indicates no channel
    ar_Xc : (N,) float array
        The length of the channel in a cell, this can be different
        from the lateral dimension of the grid cell if the channel
        runs along the cell diagonal (:math:`m`)
    ar_dam : (N,) int array
        Switch indicating whether the current cell contains a dam. A
        value of `1` indicates a dam cell, `0` indicates no dam. The
        switch currently has no influence in the model
    ar_tan_beta : (N,) float array
        The tangent of the surface slope for each cell. The surface
        slope affects processes in the overland and soil stores
    ar_tan_beta_channel : (N,) float array
        The tangent of the channel slope for each cell. This is
        conceptually different from the surface slope and affects the
        channel store
    ar_L : (N,) float array
        The depth of the soil store in each cell (:math:`m`)
    ar_Ks : (N,) float array
        The saturated hydraulic conductivity of each cell
        (:math:`mm/s`)
    ar_theta_r : (N,) float array
        The residual soil moisture content for each cell
    ar_theta_s : (N,) float array
        The saturated soil moisture content for each cell
    ar_n_o : (N,) float array
        Manning's roughness coefficient for overland flows in each
        cell
    ar_n_c : (N,) float array
        Manning's roughness coefficient for channel flows in each cell
    ar_cell_down : (N,) int array
        The label (from `ar_cell_label`) associated with the cell
        downstream of the current cell
    ar_pVs_t0 : (N,) float array
        The initial saturation of each soil store (%)
    ar_Vo_t0 : (N,) float array
        The initial volume of water in each overland store
        (:math:`m^3`)
    ar_Qc_t0 : (N,) float array
        The initial channel discharge for each cell, if applicable
        (:math:`m^3/s`)
    ar_kc : (N,) float array
        The crop co-efficient for each cell

    """
    tab_read = np.loadtxt(file_name)
    tab_read = np.atleast_2d(tab_read)

    ar_cell_label=np.array(tab_read[:,0], int)
    ar_coorx=tab_read[:,1]
    ar_coory=tab_read[:,2]
    ar_lambda=np.array(tab_read[:,3], int)
    ar_Xc=np.array(tab_read[:,4])
    ar_dam=np.array(tab_read[:,5], int)
    ar_tan_beta=tab_read[:,6]
    ar_tan_beta_channel=tab_read[:,7]
    ar_L=tab_read[:,8]
    ar_Ks=tab_read[:,9]
    ar_theta_r=tab_read[:,10]
    ar_theta_s=tab_read[:,11]
    ar_n_o=tab_read[:,12]
    ar_n_c=tab_read[:,13]
    ar_cell_down=np.array(tab_read[:,14], int)
    ar_pVs_t0=tab_read[:,15]
    ar_Vo_t0=tab_read[:,16]
    ar_Qc_t0=tab_read[:,17]
    ar_kc=tab_read[:,18]
    psi_b=tab_read[:,19]
    lamda=tab_read[:,20]

    return ar_cell_label, ar_coorx,ar_coory, ar_lambda, \
           ar_Xc,ar_dam, ar_tan_beta, ar_tan_beta_channel, \
           ar_L,ar_Ks, ar_theta_r, ar_theta_s, \
           ar_n_o, ar_n_c, ar_cell_down, \
           ar_pVs_t0, ar_Vo_t0, ar_Qc_t0, ar_kc, psi_b, lamda

def sort_cell(ar_cell_label, ar_cell_down):
    """Determine a suitable computation order for the cells.

    Sort the cells into a valid computation order based on the
    distance to the outlet of the catchment. This is necessary as each
    cell depends on the contribution from upstream cells.

    Parameters
    ----------
    ar_cell_label : (N,) int array
        Numbers labelling each cell
    ar_cell_down : (N,) int array
        The label (from `ar_cell_label`) associated with the cell
        downstream of the current cell

    Returns
    -------
    ar_label_sort : (N,) int array
        Sorted array of cell labels

    """
    network_dag = _generate_network_dag(ar_cell_label, ar_cell_down)

    ar_label_sort = np.array(nx.topological_sort(network_dag))

    return ar_label_sort

def _generate_network_dag(nodes, downstream_nodes):
    """DAG from list of nodes and downstream nodes.

    Generates a networkx Directed Acyclic Graph (DAG) from a list of
    cell nodes and their corresponding downstream neighbours.

    Parameters
    ----------
    nodes : (N,) int array
        Integer IDs labelling each cell
    downstream_nodes : (N,) int array
        The ID associated with the cell downstream of the current cell

    Returns
    -------
    dag : Networkx DiGraph
        A networkx DiGraph instance describing the catchment topology.

    """
    edges = []
    for k in nodes:
        if downstream_nodes[k] >= 0:
            edges.append((k, downstream_nodes[k]))

    dag = nx.DiGraph()
    dag.add_nodes_from(nodes)
    dag.add_edges_from(edges)

    return dag

def direct_up_cell(ar_cell_label, ar_cell_down, ar_label_sort):
    """Calculate the upstream cells for each cell.

    This function calculates the immediate upstream cells contributing
    flow to each cell in the catchment.

    Parameters
    ----------
    ar_cell_label : (N,) int array
        Numbers labelling each cell
    ar_cell_down : (N,) int array
        The label (from `ar_cell_label`) associated with the cell
        downstream of the current cell
    ar_label_sort : (N,) int array
        Sorted array of cell labels

    Returns
    -------
    li_cell_up : List of (M,) int arrays
        A list of arrays containing the upstream cell labels for each
        cell in a catchment

    """
    nb_cell=len(ar_label_sort)
    li_cell_up=[]
    for ncell in range(len(ar_label_sort)):
        cell_label=ar_cell_label[ncell]
        up_cell=ar_cell_label[np.where(ar_cell_down==cell_label)]
        li_cell_up.append(up_cell)

    return li_cell_up

def drained_area(ar_label_sort, li_cell_up, X):
    """Compute the drained area for each cell.

    This function calculates the total area drained for each cell in
    the catchment, as the sum of it's area and the upstream area.

    Parameters
    ----------
    ar_label_sort : (N,) int array
        Sorted array of cell labels
    li_cell_up : List of (M,) int arrays
        A list of arrays containing the upstream cell labels for each
        cell in a catchment
    X : scalar
        The lateral dimension of the grid-cell (:math:`m`)

    Returns
    -------
    ar_A_drained : (N,) float array
        The total drained area associated with each cell (:math:`m^2`)

    """
    A_cell=X**2
    nb_cell=len(ar_label_sort)
    ar_A_drained=np.ones(nb_cell)*-99.9
    for cell in ar_label_sort:
        up_cell=li_cell_up[cell]
        ar_A_drained[cell]=sum(ar_A_drained[up_cell])+A_cell

    return ar_A_drained

def compute_cell_param(X, ar_Xc, Dt, alpha_s, alpha_o,
                       alpha_c, nb_cell, A_thres, W_max,
                       W_min, ar_lambda, ar_tan_beta,
                       ar_tan_beta_channel, ar_L, ar_Ks,
                       ar_theta_r, ar_theta_s, ar_n_o,
                       ar_n_c, ar_A_drained):
    """Compute model parameters from physical parameters.

    This function uses the physically based parameters and constraints
    for each model cell to compute, the saturated soil moisture
    volume, channel width and constant terms for the differential
    equations of the soil, overland and channel stores.

    Parameters
    ----------
    X : scalar
        The lateral dimension of the grid-cell (:math:`m`)
    ar_Xc : (N,) float array
        The length of the channel in a cell, this can be different
        from the lateral dimension of the grid cell if the channel
        runs along the cell diagonal (:math:`m`)
    Dt : scalar
        The length of the current time-step in seconds
    alpha_s : scalar
        The dimensionless pore-size distribution parameter for the
        soil store
    alpha_o : scalar
        Power co-efficient for Mannings Equation applied to the
        overland store the value is typically 5/3
    alpha_c : scalar
        Power co-efficient for Mannings Equation applied to the
        channel store the value is typically 5/3
    nb_cell : scalar
        The number of cells in the catchment
    A_thres : scalar
        The minimum area of upstream contributing cells required
        before a cell is considered to initiate a river channel
        (:math:`m^2`)
    W_min : scalar
        The minimum width of a channel (:math:`m`)
    W_max : scalar
        The maximum width of a channel (:math:`m`)
    ar_lambda : (N,) int array
        Switch indicating whether the current cell contains a
        channel. A value of `1` indicates a channel cell, `0`
        indicates no channel
    ar_tan_beta : (N,) float array
        The tangent of the surface slope for each cell. The surface
        slope affects processes in the overland and soil stores
    ar_tan_beta_channel : (N,) float array
        The tangent of the channel slope for each cell. This is
        conceptually different from the surface slope and affects the
        channel store
    ar_L : (N,) float array
        The depth of the soil store in each cell (:math:`m`)
    ar_Ks : (N,) float array
        The saturated hydraulic conductivity of each cell
        (:math:`mm/s`)
    ar_theta_r : (N,) float array
        The residual soil moisture content for each cell
    ar_theta_s : (N,) float array
        The saturated soil moisture content for each cell
    ar_n_o : (N,) float array
        Manning's roughness coefficient for overland flows in each
        cell
    ar_n_c : (N,) float array
        Manning's roughness coefficient for channel flows in each cell
    ar_A_drained : (N,) float array
        The total drained area associated with each cell (:math:`m^2`)

    Returns
    -------
    ar_Vsm : (N,) float array
        The saturated moisture volume of the soil store for each cell
        (:math:`m^3`)
    ar_b_s : (N,) float array
        The constant term of the non differential equation for each
        soil store --> :math:`dV_s/dt = a_s - b_s V_s^{\\alpha_s}`
    ar_b_o : (N,) float array
        The constant term of the non differential equation for each
        overland store --> :math:`dV_o/dt = a_o - b_o V_o^{\\alpha_o}`
    ar_W : (N,) float array
        The channel width for each cell (:math:`m`)
    ar_b_c : (N,) float array
        The constant term of the non differential equation for each
        channel store --> :math:`dV_c/dt = a_c - b_c V_c^{\\alpha_c}`

    """
    ##Soil parameters
    ar_Vsm=(ar_theta_s-ar_theta_r)*(X**2)*ar_L
    ar_Cs=(ar_L*ar_Ks*ar_tan_beta)/(((ar_theta_s-ar_theta_r)*ar_L)**alpha_s)
    ar_b_s=ar_Cs*X/(X**(2*alpha_s))
    ##Overland parameters
    ar_Co=(1/ar_n_o)*(ar_tan_beta)**0.5
    ar_b_o=ar_Co*X/(X**(2*alpha_o))
    ##Channel parameters
    A_total=nb_cell*X**2
    ar_W = W_max + ((W_max-W_min)/(A_total**0.5-A_thres**0.5)) \
           * (ar_A_drained**0.5-A_total**0.5)

    ar_Cc=(1/ar_n_c)*(ar_tan_beta_channel)**0.5

    ar_b_c=ar_Cc*ar_W/((ar_Xc*ar_W)**(alpha_c))
    ar_W[ar_lambda==0]=-99.9
    ar_b_c[ar_lambda==0]=-99.9

    return ar_Vsm, ar_b_s, ar_b_o, ar_W, ar_b_c

#``````````````````````````````````````````
def read_column_input(file_name,nb_cell):
    """ read_column_input

        Read the file containing data in column format:
        Cell1  Cell2  Cell3 ...
        1.3    4.3     5.2  ...
        2.3    5.6     4.2  ...

        Return a matrix mat_out(nrow, ncol)
    """
    file_read=open(file_name,'r')
    tab_read=file_read.readlines()

    nb_time_step=len(tab_read)-1
    a=np.zeros(nb_cell*nb_time_step)
    mat_out=a.reshape(nb_time_step,nb_cell)

    i=-1
    for line in tab_read[1:]:
        i=i+1
        donnees=line.split()
        mat_out[i,]=[float(elem) for elem in donnees]
    ##end of iteration the file is closed
    file_read.close()

    return mat_out

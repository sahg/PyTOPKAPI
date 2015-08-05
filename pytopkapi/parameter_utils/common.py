import numpy as np
import pandas as pd

def parameter_file_to_dataframe(param_fname):
    column_labels = ['cell_label',
                     'X',
                     'Y',
                     'lambda',
                     'Xc',
                     'dam',
                     'tan_beta',
                     'tan_beta_channel',
                     'L',
                     'Ks',
                     'theta_r',
                     'theta_s',
                     'n_o',
                     'n_c',
                     'cell_down',
                     'pVs_t0',
                     'ar_Vo_t0',
                     'Qc_t0',
                     'Kc',
                     'psi_b',
                     'lambda']

    params = np.loadtxt(param_fname)

    params = pd.DataFrame(params[:, 1:],
                          index=params[:, 0], columns=column_labels[1:])

    return params

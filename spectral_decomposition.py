#################
#  Calculate spectral decomposition given connectivity matrices
#  Connectivity matrices from https://doi.org/10.1371/journal.pone.0014832
#  Pablo F. Damasceno Aug 2019
#  pablo.damasceno@ucsf.edu
################

from scipy.io import loadmat
import numpy as np
import os
# import matplotlib.pyplot as plt

def get_connectivity_matrices(filepath):
    '''
    read connectivity matrices from .mat file
    
    Arguments:
        filepath {str} -- full path to .mat file containing connectivity matrices
    '''
    from pathlib import Path
    my_file = Path(filepath)
    if my_file.is_file():
        healthy_14    = loadmat(my_file)['S'] #load data for all 14 patients
        conn_matrices = np.transpose(healthy_14)
        return(conn_matrices)
    else:
        raise(FileNotFoundError('Is path to file correct??'))
    
def get_path():
    '''
    returns path to connectivity_matrices.mat
    '''

    here_dir          = os.path.dirname(os.path.realpath('__file__'))
    conn_matrices_dir = os.path.join(here_dir, 'data', 'connectivity_matrices.mat')
    return(conn_matrices_dir)

def calculate_laplacian(connec_matrix):
    '''[Calculate Laplacian matrix given the Connectivity Matrix C: 
    L = D - C, where D is the Identity Matrix (diagonal)]
    
    Arguments:
        connectivity_matrix {np.array} -- [weighted matrix of connections between regions in ATLAS]
    '''
    degrees     = [np.sum(connec_matrix[i]) for i in range(len(connec_matrix[0]))]
    diag_matrix = np.diag(degrees)
    lap_matrix  = diag_matrix - connec_matrix
    return(lap_matrix)

    if __name__ == 'main':
        np.sum(laplacian[0])


# get array of connectivity matrices for all 14 patients
conn_matrices_dir = get_path()
conn_matrices     = get_connectivity_matrices(conn_matrices_dir)

laplacian = calculate_laplacian(conn_matrices[0])

np.sum(laplacian[0])
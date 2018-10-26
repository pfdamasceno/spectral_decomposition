###########################################################################
#  Calculate spectral cluster of brain regions based on their connections
#  connectivity matrices from https://doi.org/10.1371/journal.pone.0014832
#  prunning code from https://doi.org/10.1371/journal.pone.0035029.g002
###########################################################################

from scipy.io import loadmat
import numpy as np
import os
import sklearn.cluster
import matplotlib.pyplot as plt
import numpy as np
from nilearn import plotting
import pandas as pd

def get_path(filename):
    '''
    returns path to file: './data/filename'
    '''

    here_dir = os.path.dirname(os.path.realpath('__file__'))
    file_dir = os.path.join(here_dir, 'data', filename)
    return(file_dir)

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
        print(my_file)
        raise(FileNotFoundError('Is path to file correct??'))

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

def prune_connectivity_matrices(conn_matrices, p = 0.0001):
    '''[prune connectivity matrix based on PLoS paper
    doi:10.1371/journal.pone.0035029.g002]

    Arguments:
        conn_matrices {np.array} -- [array of weighted connectivity matrices]
    '''

    mean_matrix = np.mean(conn_matrices, axis=0)
    shifted_connectivity_matrices  = conn_matrices - mean_matrix


    def prune_matrix(matrix, z_matrix, p):
        '''[]

        Arguments:
            matrix {np.array} -- [matrix to be pruned]
            z_matrix {np.array} -- [matrix of z_scores for each node]
            p {double} -- [z-score below which connections will be pruned]
        '''

        indices_for_pruning = np.argwhere(np.abs(z_matrix) < p)
        for [i,j] in indices_for_pruning:
            matrix[int(i)][int(j)] = 0.000
        return(matrix)

    num_patients = range(len(conn_matrices))

    for i in range(10):

        variance_matrix = np.var(shifted_connectivity_matrices, axis=0)
        z_matrices      = [np.nan_to_num(shifted_connectivity_matrices[i]/variance_matrix) for i in num_patients]
        pruned_matrices = np.array([prune_matrix(conn_matrices[i], z_matrices[i], p) for i in num_patients])

        distance = np.sum(pruned_matrices**2 - conn_matrices**2)

        mean_matrix = np.mean(pruned_matrices, axis=0)
        shifted_connectivity_matrices  = pruned_matrices - mean_matrix

    return(pruned_matrices)

def plot_glass_brains(color = color):
    coordinates_path = get_path('coordinates.csv')
    df = pd.read_csv(coordinates_path)
    color = color

    coords = [(df['X'][i],df['Y'][i],df['Z'][i]) for i in range(90)]
    connec = np.array([[0]*90]*90)

    plotting.plot_connectome(connec, coords, node_color=color, display_mode='lyrz')

if __name__ == 'main':
    conn_matrices_dir = get_path('connectivity_matrices.mat')
    conn_matrices     = get_connectivity_matrices(conn_matrices_dir)

    pruned = prune_connectivity_matrices(conn_matrices)
    shifted_conn_matrices = pruned + 0.0001
    inv_conn_matrices = 1/ shifted_conn_matrices

    kernelized = np.exp(- inv_conn_matrices ** 2 / (2. * 0.003 ** 2))

    spectral = sklearn.cluster.SpectralClustering(n_clusters=2, eigen_solver='arpack')
    spectral.fit(kernelized[0])

    color = spectral.labels_

    plot_glass_brains(color)

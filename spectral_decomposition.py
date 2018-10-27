###########################################################################
#  Calculate spectral cluster of brain regions based on their connections
#  connectivity matrices from https://doi.org/10.1371/journal.pone.0014832
#  prunning code from https://doi.org/10.1371/journal.pone.0035029.g002
#
# TODO:
# 1. get submatrix of nodes based on color
# 2. iterate over submatrices to calculate decomposition of decompositions
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
    """Find filename in ./data/ directory.

    Args:
        filename (str): file we're looking for in the ./data/ directory.

    Returns:
        str: path to file "filename" in ./data/ dir.

    """

    here_dir = os.path.dirname(os.path.realpath('__file__'))
    file_dir = os.path.join(here_dir, 'data', filename)
    return(file_dir)

def get_connectivity_matrices(filepath):
    """Read connectivity matrices from .mat file.

    Args:
        filepath (str): full path to .mat file containing connectivity matrices.

    Returns:
        arr: Array of connectivity matrices, one per patient.

    """

    from pathlib import Path
    my_file = Path(filepath)

    if my_file.is_file():
        healthy_14    = loadmat(my_file)['S'] #load data for all 14 patients
        conn_matrices = np.transpose(healthy_14)
        return(conn_matrices)

    else:
        print(my_file)
        raise(FileNotFoundError('Is path to file correct??'))

def prune_connectivity_matrices(conn_matrices, p = 0.0001):
    """prune connectivity matrix based on PLoS paper
                doi:10.1371/journal.pone.0035029.g002.

    Args:
        conn_matrices (arr): Array of connectivity matrices (90 x 90 x num_patients).
        p (double): Prunning parameter below which connections will be ignored.

    Returns:
        arr: Array of pruned connectivity matrices.

    """

    mean_matrix = np.mean(conn_matrices, axis=0)
    shifted_connectivity_matrices  = conn_matrices - mean_matrix


    def prune_matrix(matrix, z_matrix, p):
        """Prune connections whose z-score are below cutoff parameter `p`.

        Args:
            matrix (arr): Connectivity matrix to be prunned.
            z_matrix (arr): Matrix of z-scores based on all patients.
            p (double): Prunning parameter below which connections will be ignored.

        Returns:
            arr: Prunned matrix.

        """

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

def calculate_laplacian(connec_matrix):
    """Calculate Laplacian matrix given the Connectivity Matrix C:
            L = D - C, where D is the Identity Matrix (diagonal).

    Args:
        connec_matrix (arr): Array of connectivity matrices (90 x 90 x num_patients).

    Returns:
        arr: Laplacian matrix.

    """

    degrees     = [np.sum(connec_matrix[i]) for i in range(len(connec_matrix[0]))]
    diag_matrix = np.diag(degrees)
    lap_matrix  = diag_matrix - connec_matrix

    return(lap_matrix)

def get_submatrix(matrix, index_list):
    """Remove `index_list` rows and columns from `matrix`.

    Args:
        matrix (arr): Matrix to calculate submatrix from.
        index_list (list): list of indices of rows / columns to remove from Matrix.

    Returns:
        arr: original `Matrix` with rows and columns removed.

    """

    submatrix = np.delete(matrix, index_list,0)
    submatrix = np.delete(submatrix, index_list,1)
    return(submatrix)

def plot_glass_brains(color):
    """Plot a glass brain for a 90 regions ATLAS with nodes colored by `color`.

    Args:
        color (list): Color indices. e.g. [0,1,1,0] will color nodes [1,2] differently.

    Returns:
        matplotlib.plot: A plot object.

    """

    coordinates_path = get_path('coordinates.csv')
    df = pd.read_csv(coordinates_path)
    color = color

    coords = [(df['X'][i],df['Y'][i],df['Z'][i]) for i in range(90)]
    connec = np.array([[0]*90]*90)

    plotting.plot_connectome(connec, coords, node_color=color, display_mode='lyrz')

conn_matrices_dir = get_path('connectivity_matrices.mat')
conn_matrices     = get_connectivity_matrices(conn_matrices_dir)

pruned = prune_connectivity_matrices(conn_matrices)
shifted_conn_matrices = pruned + 0.0001
inv_conn_matrices = 1/ shifted_conn_matrices

kernelized = np.exp(- inv_conn_matrices ** 2 / (2. * 0.005 ** 2))

spectral = sklearn.cluster.SpectralClustering(n_clusters=2, eigen_solver='arpack')
spectral.fit(kernelized[0])

color_step1 = spectral.labels_
color_step1

plot_glass_brains(color_step1)

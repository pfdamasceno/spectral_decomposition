#################################
#  Calculate spectral decomposition for the brain given connectivity matrices
#  Connectivity matrices from https://doi.org/10.1371/journal.pone.0014832
#  Pablo F. Damasceno Aug 2019
#  pablo.damasceno@ucsf.edu
################################

from scipy.io import loadmat
import numpy as np
import os
import sklearn.cluster
import matplotlib.pyplot as plt
import numpy as np
from nilearn import plotting
import pandas as pd

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

def get_path(filename):
    '''
    returns path to connectivity_matrices.mat
    '''

    here_dir = os.path.dirname(os.path.realpath('__file__'))
    file_dir = os.path.join(here_dir, 'data',filename)
    return(file_dir)

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

# get array of connectivity matrices for all 14 patients
conn_matrices_dir = get_path('connectivity_matrices.mat')
conn_matrices     = get_connectivity_matrices(conn_matrices_dir)

laplacian = calculate_laplacian(conn_matrices[0])

spectral = sklearn.cluster.SpectralClustering(n_clusters=2)
spectral.fit(conn_matrices[0])
color = spectral.labels_

def plot_glass_brains(color = color):
    coordinates_path = get_path('coordinates.csv')

    color = color

    coords = [(df['X'][i],df['Y'][i],df['Z'][i]) for i in range(90)]
    connec = np.array([[0]*90]*90)

    plotting.plot_connectome(connec, coords, node_color=color, display_mode='lyrz')


plot_glass_brains(color)

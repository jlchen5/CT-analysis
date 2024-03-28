# -*- coding: utf-8 -*-
"""
Created on Fri Nov 19 10:10:23 2021

@author: scalvinib
"""


import numpy as np
import matplotlib.pyplot as plt
import PIL.Image
import networkx as nx
import pandas as pd
import string

from functions.genome_topology import open_pdb
from functions.genome_topology import select_chrom
from functions.genome_topology import geom_distance
from functions.genome_topology import make_graph
from functions.genome_topology import fractal_dimension
from functions.genome_topology import get_matrix
from functions.genome_topology import normalize_psc

def load_pdb(file_path):
    """
    Load a PDB file and extract atom coordinates.

    Parameters:
    file_path (str): Path to the PDB file.

    Returns:
    int: The number of atoms.
    list of tuples: List where each tuple represents the coordinates (x, y, z) of an atom.
    """
    coordinates = []  # List to hold the coordinates
    with open(file_path, 'r') as pdb_file:
        for line in pdb_file:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                # Extract the x, y, and z coordinates for each atom
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                coordinates.append((x, y, z))  # Add the coordinates as a tuple to the list
    
    n = len(coordinates)  # Number of atoms
    return n, coordinates  # Return the number of atoms and the list of coordinates


def main(parameters, path_data,path_results):
    letters=list(string.ascii_lowercase)
    chr_vec=['chr {}'.format(letter) for letter in letters[:n_all_chr]]
    cell=path_data[-5:] 
    
    for n_chr, chrom in enumerate(chr_vec):
        print('analyzing the chrom',chrom)
        
        # n, coord= select_chrom(n_chr, path_data) 
        n, coord_list = load_pdb(path_data)
        coord = np.array(coord_list)  # Convert list of tuples into NumPy array       
        iterations=np.int(n/parameters['resolution'])
        start_iteration= np.int(parameters['init']/parameters['resolution'])
        
        Parallel=np.zeros(iterations)
        Series=np.zeros(iterations)
        Cross=np.zeros(iterations)
        Dim_fractal=np.zeros(iterations)
        r2_fractalfit=np.zeros(iterations)
        N_contacts=np.zeros(iterations)
        clustering= np.zeros(iterations)
        
        for t in range(start_iteration,iterations):
            
            n_atoms=np.int(resolution*(t+1))
            print('analyzing the n_tom',n_atoms)
            coord_cut=coord[0: n_atoms]  
            dist, N_contacts[t], index=geom_distance(coord_cut, 
            parameters['cutoff'], parameters['neighbors'])
             
            try:
                
                mat, stats = get_matrix(index,chrom)
                Parallel[t], Series[t], Cross[t]=normalize_psc(stats,N_contacts[t])
                Dim_fractal[t], r2_fractalfit[t]=fractal_dimension(mat, plot_fig=0)
                G=make_graph(index)
                clustering[t]= nx.average_clustering(G)
            except:
                print('WARNING: NOT ENOUGH CONTACTS FOR ANALYSIS')
        
        topology_parameters = {'Parallel (%)':Parallel, 'Series (%)':Series, 
        'Cross (%)':Cross, 'N contacts': N_contacts,
        'Fractal dimension':Dim_fractal, 'r squared': r2_fractalfit,
        'Clustering': clustering}
        topology_parameters= pd.DataFrame(topology_parameters)
        topology_parameters.to_csv('{}/Top_parameters_{}_{}.csv'.format(
            path_results, cell, chrom))

     
        
        
r_cutoff=1.0
neighbours=1
n_all_chr=20
resolution=5
start_from= 0
parameters={'cutoff':r_cutoff,'neighbors': neighbours,'N chromosomes': 
              n_all_chr, 'resolution': resolution, 'init': start_from}

path_data= 'data/pdb/Cell_ID-01/chr_a.pdb'
path_results= 'results/cumulative analysis/{}'.format(path_data[-20:-10])      
        
main(parameters, path_data,path_results)
# -*- coding: utf-8 -*-
"""
Optimized Script for Analyzing a Single Chromosome in Genome Topology
"""

import numpy as np
import pandas as pd
import os
import networkx as nx

from functools import lru_cache

# Assuming the 'functions.genome_topology' module is properly defined and accessible
from functions.genome_topology import select_chrom, geom_distance, make_graph, fractal_dimension, get_matrix, normalize_psc

cache = {}

def cached_geom_distance(coord_cut, cutoff, neighbors):
    # Convert array to a hashable key
    key = (tuple(coord_cut.flatten()), cutoff, neighbors)
    
    if key in cache:
        return cache[key]
    else:
        result = geom_distance(coord_cut, cutoff, neighbors)
        cache[key] = result
        return result


def process_chromosome(parameters, path_data, resolution, chrom="a"):
    path_results, cell = parameters['path_results'], parameters['cell']
    
    print(f"Processing the chromosome {chrom}")
    
    # Assuming `select_chrom` can directly accept the chromosome identifier (e.g., "a") as an argument
    n, coord = select_chrom(chrom, path_data)
    n, coord = select_chrom(chrom, path_data)
    print(f"Number of coordinates: {n}, Sample coordinates: {coord[:5]}")  # Print a sample of coordinates
    iterations = int(n / parameters['resolution'])
    start_iteration = int(parameters['init'] / parameters['resolution'])
    print(f"Iterations: {iterations}, Start from: {start_iteration}")

    
    Parallel = np.zeros(iterations)
    Series = np.zeros(iterations)
    Cross = np.zeros(iterations)
    Dim_fractal = np.zeros(iterations)
    r2_fractalfit = np.zeros(iterations)
    N_contacts = np.zeros(iterations)
    clustering = np.zeros(iterations)
    
    for t in range(start_iteration, iterations):
        n_atoms = int(resolution * (t + 1))
        coord_cut = coord[:n_atoms]
        dist, N_contacts[t], index = cached_geom_distance(coord_cut, parameters['cutoff'], parameters['neighbors'])
        
        try:
            mat, stats = get_matrix(index, chrom)
            Parallel[t], Series[t], Cross[t] = normalize_psc(stats, N_contacts[t])
            Dim_fractal[t], r2_fractalfit[t] = fractal_dimension(mat, plot_fig=0)
            G = make_graph(index)
            clustering[t] = nx.average_clustering(G)
        except Exception as e:
            print(f'WARNING: NOT ENOUGH CONTACTS FOR ANALYSIS ON {chrom}. Error: {e}')
    
    topology_parameters = {'Parallel (%)': Parallel, 'Series (%)': Series, 'Cross (%)': Cross, 'N contacts': N_contacts, 'Fractal dimension': Dim_fractal, 'r squared': r2_fractalfit, 'Clustering': clustering}
    topology_parameters = pd.DataFrame(topology_parameters)
    topology_parameters.to_csv(f'{path_results}/Top_parameters_{cell}_{chrom}.csv')
    
    return f"Processed {chrom}"

if __name__ == "__main__":
    r_cutoff = 1.0
    neighbours = 1
    resolution = 5
    start_from = 0
    parameters = {'cutoff': r_cutoff, 'neighbors': neighbours, 'resolution': resolution, 'init': start_from}

    path_data = '/Users/jialechen/Desktop/PhD/CT/Pang_2022_GenomeBiol_3D/data/pdb/Cell_ID-03'
    path_results = '/Users/jialechen/Desktop/PhD/CT/Pang_2022_GenomeBiol_3D/results/cumulative analysis/'
    parameters['path_results'] = path_results
    parameters['cell'] = path_data.split('/')[-1]

    result = process_chromosome(parameters, path_data, parameters['resolution'], chrom="a")
    print(result)

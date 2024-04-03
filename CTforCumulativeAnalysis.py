import numpy as np
import pandas as pd
import os
import networkx as nx

# Assuming the 'functions.genome_topology' module is properly defined and accessible
from functions.genome_topology import select_chrom, geom_distance, make_graph, fractal_dimension, get_matrix, normalize_psc

cache = {}

def cached_geom_distance(coord_cut, cutoff, neighbors):
    key = (tuple(coord_cut.flatten()), cutoff, neighbors)
    
    if key in cache:
        return cache[key]
    else:
        result = geom_distance(coord_cut, cutoff, neighbors)
        cache[key] = result
        return result

def process_chromosome(parameters, path_data, resolution, chrom="a"):
    path_results, cell = parameters['path_results'], parameters['cell']
    
    print(f"Processing the chromosome {chrom} for {cell}")
    
    n, coord = select_chrom(chrom, path_data)
    print(f"Number of coordinates: {n}, Sample coordinates: {coord[:5]}")
    iterations = int(n / parameters['resolution'])
    start_iteration = int(parameters['init'] / parameters['resolution'])
    print(f"Iterations: {iterations}, Start from: {start_iteration}")

    # Initialize arrays to store the results
    Parallel, Series, Cross, Dim_fractal, r2_fractalfit, N_contacts, clustering = [np.zeros(iterations) for _ in range(7)]
    
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
    
    topology_parameters = pd.DataFrame({
        'Parallel (%)': Parallel, 
        'Series (%)': Series, 
        'Cross (%)': Cross, 
        'N contacts': N_contacts, 
        'Fractal dimension': Dim_fractal, 
        'r squared': r2_fractalfit, 
        'Clustering': clustering
    })
    topology_parameters.to_csv(f'{path_results}/Top_parameters_{cell}_{chrom}.csv')
    
    return f"Processed {chrom} for {cell}"

def analyze_cells(start_id=1, end_id=24, chrom="a"):
    r_cutoff = 1.0
    neighbours = 1
    resolution = 5
    start_from = 0
    base_path_data = '/Users/jialechen/Desktop/PhD/CT/Pang_2022_GenomeBiol_3D/data/pdb/Cell_ID-'
    path_results = '/Users/jialechen/Desktop/PhD/CT/Pang_2022_GenomeBiol_3D/results/cumulative analysis/'

    for cell_id in range(start_id, end_id + 1):
        formatted_cell_id = f"{cell_id:02}"  # Ensure the cell ID is properly formatted with leading zeros
        path_data = f"{base_path_data}{formatted_cell_id}"
        parameters = {
            'cutoff': r_cutoff, 
            'neighbors': neighbours, 
            'resolution': resolution, 
            'init': start_from, 
            'path_results': path_results, 
            'cell': f"Cell_ID-{formatted_cell_id}"
        }
        result = process_chromosome(parameters, path_data, resolution, chrom=chrom)
        print(result)

if __name__ == "__main__":
    analyze_cells()  # Call the function to start processing

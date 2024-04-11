import numpy as np
import matplotlib.pyplot as plt
import PIL.Image
import networkx as nx
import pandas as pd
import string
import concurrent.futures
import os

from functions.genome_topology import open_pdb
from functions.genome_topology import select_chrom
from functions.genome_topology import geom_distance
from functions.genome_topology import make_graph
from functions.genome_topology import fractal_dimension
from functions.genome_topology import get_matrix
from functions.genome_topology import normalize_psc

def main(parameters, path_data, path_results):
    letters = list(string.ascii_lowercase)
    chr_vec = ['chr {}'.format(letter) for letter in letters[:parameters['N chromosomes']]]
    cell = path_data.split('/')[-1]  # Corrected to remove trailing slash issue
    print("Analyzing cell:", cell)
    
    for n_chr, chrom in enumerate(chr_vec):
        print("Analyzing:", chrom)
        file_name = f"chrom{chrom[-1]}.pdb"  # Assumes file names like 'chroma.pdb', 'chromb.pdb', etc.
        file_path = os.path.join(path_data, file_name)

        if not os.path.exists(file_path):
            print(f"File not found: {file_path}, skipping...")
            continue  # Skip this chromosome if the file does not exist

        n, coord = select_chrom(n_chr, path_data)        
        iterations = int(n / parameters['resolution'])
        start_iteration = int(parameters['init'] / parameters['resolution'])
        
        Parallel = np.zeros(iterations)
        Series = np.zeros(iterations)
        Cross = np.zeros(iterations)
        Dim_fractal = np.zeros(iterations)
        r2_fractalfit = np.zeros(iterations)
        N_contacts = np.zeros(iterations)
        clustering = np.zeros(iterations)
        
        for t in range(start_iteration, iterations):
            n_atoms = int(parameters['resolution'] * (t + 1))
            print("Analysis the N_atoms:{}".format(n_atoms))
            coord_cut = coord[:n_atoms]
            dist, N_contacts[t], index = geom_distance(coord_cut, parameters['cutoff'], parameters['neighbors'])
            
            try:
                mat, stats = get_matrix(index, chrom)
                Parallel[t], Series[t], Cross[t] = normalize_psc(stats, N_contacts[t])
                Dim_fractal[t], r2_fractalfit[t] = fractal_dimension(mat, plot_fig=0)
                G = make_graph(index)
                clustering[t] = nx.average_clustering(G)
            except Exception as e:
                print('WARNING: NOT ENOUGH CONTACTS FOR ANALYSIS ON {}. Error: {}'.format(chrom, e))
        
        # Define the complete file path
        csv_file_path = os.path.join(path_results, f'Top_parameters_{cell}_{chrom}.csv')
        
        # Ensure the results directory exists
        os.makedirs(os.path.dirname(csv_file_path), exist_ok=True)

        # Create DataFrame
        topology_parameters = pd.DataFrame({
            'Parallel (%)': Parallel, 
            'Series (%)': Series, 
            'Cross (%)': Cross, 
            'N contacts': N_contacts,
            'Fractal dimension': Dim_fractal, 
            'r squared': r2_fractalfit,
            'Clustering': clustering
        })

        # Save DataFrame to CSV using the correct path
        topology_parameters.to_csv(csv_file_path)
        print(f"Results for {chrom} saved successfully at {csv_file_path}.")

# Parameters
r_cutoff = 1.0
neighbours = 1
n_all_chr = 23
resolution = 20
start_from = 0
parameters = {
    'cutoff': r_cutoff,
    'neighbors': neighbours,
    'N chromosomes': n_all_chr, 
    'resolution': resolution, 
    'init': start_from
}

# Base directory paths
base_data_path = '/home/cheetah/jiale/pdb/'
base_results_path = '/home/cheetah/jiale/results/cumulative analysis/'

# Thread pool execution
with concurrent.futures.ThreadPoolExecutor(max_workers=100) as executor:
    futures = []
    for i in range(1, 25):  # 25 is exclusive, up to 24
        cell_id = f"Cell_ID-{i:02}"  # Formats number with leading zeros
        path_data = os.path.join(base_data_path, cell_id)
        path_results = os.path.join(base_results_path, cell_id)
         
        
        # Submit the task to the thread pool
        futures.append(executor.submit(main, parameters, path_data, path_results))

    # Wait for all threads to complete
    for future in concurrent.futures.as_completed(futures):
        print(future.result()) 

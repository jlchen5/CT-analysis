import numpy as np
import os
import pandas as pd
import string
import networkx as nx
import concurrent.futures

from functions.genome_topology import open_pdb, select_chrom, geom_distance, make_graph, fractal_dimension, get_matrix, normalize_psc

def analyze_chromosome(parameters, path_data, path_results, cell, chrom):
    print(f"Analyzing cell: {cell}, chromosome: {chrom}")
    file_name = f"chrom{chrom[-1]}.pdb"  # Assumes file names like 'chroma.pdb', 'chromb.pdb', etc.
    file_path = os.path.join(path_data, file_name)

    if not os.path.exists(file_path):
        print(f"File not found: {file_path}, skipping...")
        return  # Skip this chromosome if the file does not exist

    n, coord = select_chrom(ord(chrom[-1]) - ord('a'), path_data)
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
        print(f"Analysis in {cell} on {chrom} at {n_atoms} atoms")
        coord_cut = coord[:n_atoms]
        dist, N_contacts[t], index = geom_distance(coord_cut, parameters['cutoff'], parameters['neighbors'])

        try:
            mat, stats = get_matrix(index, chrom)
            Parallel[t], Series[t], Cross[t] = normalize_psc(stats, N_contacts[t])
            Dim_fractal[t], r2_fractalfit[t] = fractal_dimension(mat, plot_fig=0)
            G = make_graph(index)
            clustering[t] = nx.average_clustering(G)
        except Exception as e:
            print(f'WARNING: Not enough contacts for analysis on {chrom}. Error: {e}')
            continue

    # Define the complete file path
    csv_file_path = os.path.join(path_results, f'Top_parameters_{cell}_{chrom}.csv')

    # Ensure the results directory exists
    os.makedirs(os.path.dirname(csv_file_path), exist_ok=True)

    # Create DataFrame
    topology_parameters = pd.DataFrame({
        'n atoms': np.arange(start_iteration, iterations) * parameters['resolution'],
        'Parallel': Parallel*100,
        'Series': Series*100,
        'Cross': Cros*100,
        'N contacts': N_contacts,
        'Fractal dimension': Dim_fractal,
        'r squared': r2_fractalfit,
        'Clustering': clustering
    })

    # Save DataFrame to CSV using the correct path
    topology_parameters.to_csv(csv_file_path)
    print(f"Results for {cell} and {chrom} saved successfully at {csv_file_path}.")

def main():
    # Parameters
    parameters = {
        'cutoff': 1.0,
        'neighbors': 1,
        'N chromosomes': 23,
        'resolution': 1000,
        'init': 0
    }

    # Base directory paths
    base_data_path = '/Users/jialechen/Desktop/PhD/CT/Pang_2022_GenomeBiol_3D/data/pdb/'
    base_results_path = '/Users/jialechen/Desktop/PhD/CT/Pang_2022_GenomeBiol_3D/results/cumulative analysis/'
    letters = list(string.ascii_lowercase)[:parameters['N chromosomes']]
    chr_vec = [f'chr {letter}' for letter in letters]

    # Thread pool execution
    with concurrent.futures.ProcessPoolExecutor(max_workers=os.cpu_count()) as executor:
        futures = []
        for i in range(1, 2):  # Loop over cells
            cell_id = f"Cell_ID-{i:02}"
            path_data = os.path.join(base_data_path, cell_id)
            path_results = os.path.join(base_results_path, cell_id)

            for chrom in chr_vec:
                futures.append(executor.submit(analyze_chromosome, parameters, path_data, path_results, cell_id, chrom))

        # Wait for all threads to complete
        concurrent.futures.wait(futures)
        print("All analyses complete.")

if __name__ == "__main__":
    main()

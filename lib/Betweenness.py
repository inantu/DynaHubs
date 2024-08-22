import networkx as nx
import numpy as np
import glob
import matplotlib.pyplot as plt
import os
import argparse

def process_matrix_file(file_path):
    G = nx.Graph()
    with open(file_path, 'r') as file:
        for line in file:
            segid1, node1, segid2, node2, weight = line.split()
            node1 = (segid1, int(node1))
            node2 = (segid2, int(node2))
            weight = float(weight)
            G.add_edge(node1, node2, weight=1/weight)
    return G

def process_betweenness(top_betweenness):
    result_dir = os.path.join(os.getcwd(), 'result')
    contact_matrix_dir = os.path.join(result_dir, 'contact_matrix')
    betweenness_dir = os.path.join(result_dir, 'Betweenness')

    if not os.path.exists(betweenness_dir):
        os.makedirs(betweenness_dir)

    file_pattern = os.path.join(contact_matrix_dir, 'matrix_final_*.out')
    files = glob.glob(file_pattern)

    # Aggregate betweenness centrality values across all files
    all_betweenness_centrality = []
    node_centrality_dict = {}
    all_top_betweenness_nodes = []

    total_frames = len(files)  # Total number of frames (assuming each file corresponds to a frame)

    for file_path in files:
        G = process_matrix_file(file_path)
        betweenness_centrality = nx.betweenness_centrality(G, weight='weight')
        all_betweenness_centrality.extend(betweenness_centrality.values())

        # Store betweenness centrality for each file
        node_centrality_dict[file_path] = betweenness_centrality
        
        # Determine top betweenness nodes
        quantile_0_95 = np.quantile(list(betweenness_centrality.values()), 0.95)
        top_betweenness_nodes = {node for node, centrality in betweenness_centrality.items() if centrality >= quantile_0_95}
        all_top_betweenness_nodes.extend(top_betweenness_nodes)

        # Generate an output file name based on the input file name
        output_file_path = os.path.join(betweenness_dir, os.path.basename(file_path).replace('.out', '_betweenness.txt'))
        
        # Write the betweenness centrality values to a text file
        with open(output_file_path, 'w') as f:
            for node, centrality in sorted(betweenness_centrality.items(), key=lambda x: x[1], reverse=True):
                f.write(f"{node[0]}_{node[1]}\t{centrality}\n")

        print(f"Betweenness centrality written to {output_file_path}")

    if all_betweenness_centrality:
        # Determine the 95th percentile value for betweenness centrality across all files
        quantile_0_95 = np.quantile(all_betweenness_centrality, 0.95)

        # Calculate the frequency of each node being in the top 95th percentile
        node_frequencies = {}
        for node in all_top_betweenness_nodes:
            if node in node_frequencies:
                node_frequencies[node] += 1
            else:
                node_frequencies[node] = 1

        # Convert frequencies to fraction of total frames
        node_frequencies_fraction = {node: freq / total_frames for node, freq in node_frequencies.items()}

        # Write nodes with 95th percentile betweenness centrality and their frequencies to a file
        output_file_path = os.path.join(betweenness_dir, 'top_betweenness_node_frequencies.txt')
        with open(output_file_path, 'w') as f:
            f.write("Node\tFrequency\n")
            for node, frequency in node_frequencies_fraction.items():
                f.write(f"{node[0]}_{node[1]}\t{frequency:.6f}\n")

        # Write names of the residues (nodes) to a separate text file
        residues_file_path = os.path.join(betweenness_dir, 'residue_names.txt')
        with open(residues_file_path, 'w') as f:
            for node in node_frequencies_fraction.keys():
                f.write(f"{node[0]} {node[1]}\n")

        # Print confirmation message
        print(f"Nodes and their frequencies written to {output_file_path}")

        # Plot the frequencies
        nodes = list(node_frequencies_fraction.keys())
        frequencies = list(node_frequencies_fraction.values())

        plt.figure(figsize=(10, 6))
        plt.bar([f"{node[0]}_{node[1]}" for node in nodes], frequencies, color='blue')
        plt.xlabel('Nodes')
        plt.ylabel('Frequency')
        plt.title('Frequency of Nodes in the top 5% of Betweenness Centrality')
        plt.xticks(rotation=90)
        plt.tight_layout()
        output_freq_path = os.path.join(result_dir, 'frequency_plot.png')
        plt.savefig(output_freq_path, format='png', dpi=300)
        plt.show()

        # Visualize the network and highlight top betweenness nodes
        G_combined = nx.Graph()

        # Combine all graphs into a single one
        for file_path in files:
            G = process_matrix_file(file_path)
            G_combined = nx.compose(G_combined, G)

        # Calculate betweenness centrality for the combined graph
        combined_betweenness_centrality = nx.betweenness_centrality(G_combined, weight='weight')

        # Highlight top betweenness nodes
        top_nodes = [node for node, centrality in combined_betweenness_centrality.items() if centrality >= quantile_0_95]

        # Assign a color to each segment
        segments = {node[0] for node in G_combined.nodes()}
        color_map = plt.get_cmap('tab20')
        segment_colors = {seg: color_map(i / len(segments)) for i, seg in enumerate(segments)}

        # Draw the network
        pos = nx.spring_layout(G_combined)
        plt.figure(figsize=(12, 12))

        # Draw all nodes with colors based on their segment
        node_colors = [segment_colors[node[0]] for node in G_combined.nodes()]
        nx.draw_networkx_nodes(G_combined, pos, node_size=50, node_color=node_colors)

        # Draw all edges
        nx.draw_networkx_edges(G_combined, pos, alpha=0.5)

        # Highlight top betweenness nodes in red and larger size
        nx.draw_networkx_nodes(G_combined, pos, nodelist=top_nodes, node_color='r', node_size=80)

        # Draw labels for top betweenness nodes
        labels = {node: f"{node[0]}_{node[1]}" for node in top_nodes}
        nx.draw_networkx_labels(G_combined, pos, labels, font_size=12, font_color='black')

        plt.title('Network Graph with Top Betweenness Nodes Highlighted and Segments Colored')
        # Save the network as a PNG file
        output_network_path = os.path.join(result_dir, 'network_graph.png')
        plt.savefig(output_network_path, format='png', dpi=300)

        plt.show()
    else:
        print("No betweenness centrality values to process.")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process betweenness centrality.')
    parser.add_argument('--top_betweenness', type=int, default=5, help='Top betweenness percentile (default: 5).')

    args = parser.parse_args()

    process_betweenness(args.top_betweenness)

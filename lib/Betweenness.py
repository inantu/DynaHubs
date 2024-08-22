import networkx as nx
import numpy as np
import glob
import matplotlib.pyplot as plt
import argparse
import os

def process_betweenness(top_betweenness):
    # Define the result directory based on the current working directory
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

    for file_path in files:
        G = process_matrix_file(file_path)
        betweenness_centrality = nx.betweenness_centrality(G, weight='weight')
        all_betweenness_centrality.extend(betweenness_centrality.values())

        # Store betweenness centrality for each file
        node_centrality_dict[file_path] = betweenness_centrality

        # Generate an output file name based on the input file name
        output_file_path = file_path.replace('.out', '_betweenness.txt')
        
        # Write the betweenness centrality values to a text file
        with open(output_file_path, 'w') as f:
            for node, centrality in sorted(betweenness_centrality.items(), key=lambda x: x[1], reverse=True):
                f.write(f"{node[0]}_{node[1]}\t{centrality}\n")

        print(f"Betweenness centrality written to {output_file_path}")

    # Determine the specified percentile value for betweenness centrality
    quantile = np.quantile(all_betweenness_centrality, top_betweenness / 100)

    # Collect nodes with betweenness centrality above the specified percentile
    top_betweenness_nodes = set()

    for file_path in files:
        G = process_matrix_file(file_path)
        betweenness_centrality = nx.betweenness_centrality(G, weight='weight') 
        for node, centrality in betweenness_centrality.items():
            if centrality >= quantile:
                top_betweenness_nodes.add((node, centrality))

    # Write nodes with percentile betweenness centrality and their frequencies to a file
    output_file_path = os.path.join(result_dir, 'top_betweenness_node_frequencies.txt')
    with open(output_file_path, 'w') as f:
        f.write("Node\tFrequency\n")
        for node, frequency in node_frequencies.items():
            f.write(f"{node[0]}_{node[1]}\t{frequency}\n")
    
    # Write names of the residues (nodes) to a separate text file
    residues_file_path = os.path.join(betweenness_dir, 'residue_names.txt')
    with open(residues_file_path, 'w') as f:
        for node in node_frequencies.keys():
            f.write(f"{node[0]} {node[1]}\n")
    
    # Print confirmation message
    print(f"Nodes and their frequencies written to {output_file_path}")

    # Plot the frequencies
    nodes = list(node_frequencies.keys())
    frequencies = list(node_frequencies.values())

    plt.figure(figsize=(10, 6))
    plt.bar([f"{node[0]}_{node[1]}" for node in nodes], frequencies, color='blue')
    plt.xlabel('Nodes')
    plt.ylabel('Frequency')
    plt.title('Frequency of Nodes in the top percentile of Betweenness Centrality')
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
    top_nodes = [node for node, centrality in combined_betweenness_centrality.items() if centrality >= quantile]

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

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate betweenness centrality from contact matrices.')
    parser.add_argument('--top_betweenness', type=float, default=5.0, help='Top betweenness percentile value (default: 5).')

    args = parser.parse_args()
    
    process_betweenness(args.top_betweenness)

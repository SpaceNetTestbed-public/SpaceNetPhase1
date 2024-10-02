''''

SpaceNet Graph Weight Function

AUTHOR:         Bruce Barbour
                Virginia Tech

DESCRIPTION:    This Python script weighting functions for the topology graphs.


CONTENTS:       G_weight_fun / (STARTS AT 55)
                    ...
'''

# =================================================================================== #
# ------------------------------- IMPORT PACKAGES ----------------------------------- #
# =================================================================================== #

import numpy as np
import networkx as nx

# =================================================================================== #
# ----------------------------- G WEIGHING UTILITY ---------------------------------- #
# =================================================================================== #

def G_purely_ISLs(   
                    G                       : nx.Graph,
                    num_sat                 : int,
                    source_dest_nodes       : tuple,
                    connectivity_matrix     : np.ndarray,
                    metrics                 : np.ndarray,
                ) -> nx.Graph:
    """
    Generates a weighted topology Graph that is biased towards satellite-satellite connections, with an exception to ground endpoints.

    Args:
        G (nx.Graph):                       Initial graph
        num_sat (int):                      Number of satellite nodes
        source_dest_nodes (tuple):          Source and destination endpoint nodes
        connectivity_matrix (np.ndarray):   Matrix representing the connectivity between nodes
        metrics (np.ndarray):               Matrix of point-to-point metric values to be used in graph weighing
    
    Returns:
        nx.Graph:                           Weighted, biased NetworkX graph
    """

    # Calculate edge weights
    for i in range(len(connectivity_matrix)):
        for j in range(len(connectivity_matrix[i])):
            
            # Check for a connection
            if connectivity_matrix[i][j] == 1:

                # Check if 'i' or 'j' indices are a terrestrial node
                if i >= num_sat or j >= num_sat:

                    # Check if either are an endpoint
                    if i in source_dest_nodes or j in source_dest_nodes:
                        G.add_edge(i, j, weight=int(metrics[i][j]*1e3))

                    # If not, add bias weight
                    else:
                        G.add_edge(i, j, weight=int(1e6))

                # Check if both are satellites
                elif i < num_sat and j < num_sat:
                    G.add_edge(i, j, weight=int(metrics[i][j]*1e3))
         
    # Return graph
    return G

# =================================================================================== #
# ------------------------------- MAIN FUNCTION ------------------------------------- #
# =================================================================================== #

def topology_G(   
                satellites              : list,
                source_dest_nodes       : tuple,
                connectivity_matrix     : np.ndarray,
                metrics                 : np.ndarray = None,
                criterion               : int = 2
              ) -> nx.Graph:
    """
    Initializes the NetworkX graph and assigns the weighing function.

    Args:
        satellites (list):                  List of satellite nodes
        source_dest_nodes (tuple):          Source and destination endpoint nodes
        connectivity_matrix (np.ndarray):   Matrix representing the connectivity between nodes
        metrics (np.ndarray):               Matrix of point-to-point metric values to be used in graph weighing (default=None)
        criterion (int):                    Routing criterion for corresponding edge weight assignment
                                            (0) ISL only
                                            (1) Terrestrial only
                                            (2) Integrated satellite/terrestrial (default)
    
    Returns:
        nx.Graph:                           Weighted, biased NetworkX graph
    """

    # Initialize NetworkX graph
    topology_G = nx.Graph()

    # Determine number of satellite nodes
    num_sat = len(satellites)

    # Determine number of nodes
    num_nodes = len(connectivity_matrix)

    # Check if metrics is None, then default to hops
    if not metrics:
        metrics = np.ones((num_nodes, num_nodes))

    # TERRESTRIAL ONLY =======================================================================================
    if criterion == 1:

        # Assign nodes to graph
        for n in range(num_nodes - num_sat):
            topology_G.add_node(n + num_sat)

        # Compute the edge weights
        adj_i = 0
        adj_j = 0
        for i in range(num_nodes - num_sat):
            for j in range(num_nodes - num_sat):
                adj_i = i + num_sat
                adj_j = j + num_sat
                if connectivity_matrix[adj_i][adj_j] == 1:
                    topology_G.add_edge(adj_i, adj_j, weight=int(metrics[adj_i][adj_j]*1e3))

    # SATELLITE OR ISTN ONLY =================================================================================
    elif criterion == 0 or 2:
        
        # Assign nodes to graph
        for n in range(num_nodes):
            topology_G.add_node(n)

        # Compute the edge weights
        if criterion == 0:
            topology_G = G_purely_ISLs(G=topology_G, num_sat=num_sat, source_dest_nodes=source_dest_nodes,
                                       connectivity_matrix=connectivity_matrix, metrics=metrics)
        else:
            # Mixed (ISL + Ground links)
            for i in range(num_nodes):
                for j in range(len(connectivity_matrix[i])):
                    if connectivity_matrix[i][j] == 1:
                        topology_G.add_edge(i, j, weight=int(metrics[i][j]*1e3))

    # Return weighted graph
    return topology_G
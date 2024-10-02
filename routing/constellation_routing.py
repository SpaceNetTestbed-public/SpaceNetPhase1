import time
import os

import threading

import networkx as nx
import matplotlib.pyplot as plt
import bellmanford as bf
from multiprocessing import Pool

import sys
sys.path.append("./")
from routing.routing_utils import *
from routing.G_weight_fun import *







def gs_routing_worker(data_path, gs_sat, links_updated, num_of_satellites, satellites_by_index, list_of_Intf_IPs, constellation_routes):
    """
    Performs routing for a ground station (gs) satellite.

    Args:
        data_path (str): The path to the data directory.
        gs_sat (str): The ground station satellite identifier.
        links_updated (bool): Indicates whether the links have been updated.
        num_of_satellites (int): The total number of satellites.
        satellites_by_index (list): A list of satellites indexed by their numbers.
        list_of_Intf_IPs (list): A list of interface IP addresses.
        constellation_routes (list): A list of constellation routes.

    Returns:
        None
    """
    update_gsl_routing_cmd = [] # Starts empty list to store GSL ip route commands
    gs_number = int(gs_sat[1]) % num_of_satellites
    gs_ip = get_gs_ip(list_of_Intf_IPs, "gs" + str(gs_number) + "-eth0").split("/")[0] # Get the IP address of the ground station

    gs_network_address = get_network_address(gs_ip) # Get the network address of the ground station
    #thread_list = [] # unused
    for i in range(num_of_satellites):
        route_to_sat_GW = find_route_between_src_dest(i, gs_sat[0], constellation_routes) # Find the route between the current satellite and GS satellite using list of available routes
        if route_to_sat_GW != -1: # If a route is found
            parameters = get_static_route_parameter([route_to_sat_GW], links_updated, list_of_Intf_IPs, satellites_by_index)#; # Get the static route parameters for the route using the satellites in route and the list of satellite interface IPs
            if len(parameters) > 0: # If parameters can be generated
                update_gsl_routing_cmd.append("sat" + str(i) + ",ip route add " + str(gs_network_address) + "/28 via " + str(parameters[2][:-3]) + " dev " + str(parameters[3])) # Append the GSL ip route command to the list of GSL ip route commands
            else:
                print(route_to_sat_GW) # If parameters cannot be generated, print the returned value to the console

    if len(update_gsl_routing_cmd) > 0: # If there are GSL ip route commands
        for update in update_gsl_routing_cmd:
            update_routes = update.split(",") # Split the GSL ip route command into its components (sat1539,ip route del 10.2.2.112 via 10.2.6.130/28 dev sat1539-eth4)
            file = open(data_path + "/cmd_files/" + update_routes[0] + "_routes.sh", 'a') # Open a file to write the GSL ip route commands to
            file.writelines(update_routes[1].strip() + " & \n") # Write the GSL ip route command to the file (append with & to run in background)
            file.close() # Close the file




def update_GSL_thread(sat_id, change, constellation_routes, links_updated, list_of_Intf_IPs, satellites_by_index, gs_network_address, update_gsl_routing_cmd):
    """
    UNUSED?
    Update the GSL (Ground Station Link) routing table for a given satellite.

    Args:
        sat_id (int): The ID of the satellite.
        change (tuple): A tuple containing the source and destination satellite IDs, and the change type. Format: (dest_sat_id, ???, change_type, change_type) <== VERIFY!
        constellation_routes (list): List of constellation routes.
        links_updated (bool): Flag indicating if links have been updated. <== VERIFY!
        list_of_Intf_IPs (list): List of interface IP addresses.
        satellites_by_index (dict): Dictionary mapping satellite IDs to their index.
        gs_network_address (str): The network address of the ground station.
        update_gsl_routing_cmd (list): List to store the GSL routing commands.

    Returns:
        None
    """
    route_to_sat_GW = find_route_between_src_dest(sat_id, change[0], constellation_routes) # Find the route between the current satellite and GS satellite using list of available routes
    if route_to_sat_GW != -1: # If a route is found
        parameters = get_static_route_parameter([route_to_sat_GW], links_updated, list_of_Intf_IPs, satellites_by_index) #; # Get the static route parameters for the route using the satellites in route and the list of satellite interface IPs
        if len(parameters) > 0: # If parameters can be generated
            if change[2] == 0 and change[3] == 1: # If the change[2] is type 0 AND change[3] is type 1
                update_gsl_routing_cmd.append("ip route add "+str(gs_network_address)+" via "+str(parameters[2])+" dev "+str(parameters[3])) # Append the GSL ip route 'add' command to the list of GSL ip route commands
            elif change[2] == 1 and change[3] == 0: # If the change[2] is type 1 AND change[3] is type 0
                update_gsl_routing_cmd.append("ip route del "+str(gs_network_address)+" via "+str(parameters[2])+" dev "+str(parameters[3])) # Append the GSL ip route 'del' command to the list of GSL ip route commands
    else: # If no route is found
        if sat_id != int(change[0]): # If the satellite ID is not equal to the destination satellite ID (VERIFY)
            print("Error: cannot find the route between sat", sat_id, " and sat", change[0]) # Print an error message to the console



# ================================================================================================
# FLOYD-WARSHALL ALG. IMPLEMENTATION (INITIAL ROUTING)
# ================================================================================================
def initial_routing_fw(satellites, connectivity_matrix, metric, source_dest_nodes, criterion=2):
    """
    Perform initial routing for a constellation network using Floyd-Warshall algorithm.

    Args:
        satellites (list):          List of satellite nodes.
        connectivity_matrix (list): Matrix representing the connectivity between nodes.
        metric (list):              Matrix representing the metric used for minimum-cost routing
        source_dest_nodes (tuple):  Default is None. If provided, then return output of shortest-path first route between Source/Destination nodes
        criterion (int):            Routing criterion for corresponding edge weight assignment
                                    (0) ISL only
                                    (1) Terrestrial only
                                    (2) Integrated satellite/terrestrial (default)

    Returns:
        dict: Dictionary of static routes with keys as (source, destination) tuples and values as lists of nodes in the path.
    """

    # mega_constellation_graph = nx.Graph()  # Create a new NetworkX graph for the constellation network
    # #for n in range(len(satellites) + len(ground_stations)):  # For each satellite and ground station
    # for n in range(len(connectivity_matrix)):  # For all nodes in the connectivity matrix (satellites, ground stations, and endpoints)
    #     mega_constellation_graph.add_node(n)  # Add the satellite or ground station to the graph

    # for i in range(len(connectivity_matrix)):  # For each row in the connectivity matrix
    #     for j in range(len(connectivity_matrix[i])):  # For each column in the connectivity matrix
    #         if connectivity_matrix[i][j] == 1:  # If there is a connection between the nodes
    #             mega_constellation_graph.add_edge(i, j, weight=1)#int(latency[i][j]))  # Add an edge between the nodes with a weight of 1

    mega_constellation_graph = topology_G(satellites=satellites, source_dest_nodes=source_dest_nodes, connectivity_matrix=connectivity_matrix,
                                          metrics=metric, criterion=criterion)

    # Use Floyd-Warshall algorithm to find shortest paths between all pairs of nodes
    pred, _ = nx.floyd_warshall_predecessor_and_distance(mega_constellation_graph, weight="weight")

    # Find optimal route if provided
    optimal_output = None
    if source_dest_nodes:   # Check if source exists
        src, dest = source_dest_nodes
        try:
            optimal_output = nx.reconstruct_path(src, dest, pred)
        except KeyError:
            optimal_output = None

    if criterion > 0:
        max_iter_value = len(connectivity_matrix) # Iterate over results for all nodes
    elif criterion == 0:
        max_iter_value = len(satellites) # Iterate over results and only look at satellite nodes for the routing
        
    static_routes = {}
    for i in range(max_iter_value):  # Iterate only over satellite nodes as sources
        for j in range(max_iter_value):
            if i != j and (i in pred and j in pred[i]):
                static_routes[(i, j)] = nx.reconstruct_path(i, j, pred)

    if source_dest_nodes:
        return (static_routes, optimal_output, nx.path_weight(mega_constellation_graph, optimal_output, weight="weight"))
    else:
        return static_routes



def static_routing_update_commands(static_routes, links, list_of_Intf_IPs, satellites):
    """
    Prints BASH ip route commands to stdout for every static route provided.

    Args:
        static_routes (list): List of static routes.
        links (list): List of links.
        list_of_Intf_IPs (list): List of interface IPs.
        satellites (list): List of satellites.

    Returns:
        None
    """
    for route in static_routes: # For each route in the list of static routes
        current_route = route[0]#; # Get the current route
        if len(current_route) > 2: # If the route has more than 2 nodes
            src_node, next_hop_node, dest_node, last_hop_node = current_route[0], current_route[1], current_route[len(current_route)-1], current_route[len(current_route)-2] # Set the source, next hop, destination, and last hop node variables

            src_node = "sat"+str(src_node) if src_node < len(satellites) else "gs"+str(src_node%len(satellites)) # Set the source node to the satellite or ground station number
            next_hop_node = "sat"+str(next_hop_node) if next_hop_node < len(satellites) else "gs"+str(next_hop_node%len(satellites)) # Set the next hop node to the satellite or ground station number
            dest_node = "sat"+str(dest_node) if dest_node < len(satellites) else "gs"+str(dest_node%len(satellites)) # Set the destination node to the satellite or ground station number
            last_hop_node = "sat"+str(last_hop_node) if last_hop_node < len(satellites) else "gs"+str(last_hop_node%len(satellites)) # Set the last hop node to the satellite or ground station number

            src_node_intf = ""
            dest_node_intf= ""
            next_h_node_intf = ""
            last_h_node_intf = ""

            for link in links: # For each link in the list of links
                if str(src_node)+str("-") in link and str(next_hop_node)+str("-") in link: # If the source node and next hop node are in the link
                    intfs = link.split(":") # Split the link into its component interfaces
                    if str(src_node) in intfs[0] and str(next_hop_node) in intfs[1]: # If the source node is in the first interface and the next hop node is in the second interface
                        src_node_intf = intfs[0] # Set the source node interface to the first interface
                        next_h_node_intf = intfs[1] # Set the next hop node interface to the second interface
                    elif str(src_node) in intfs[1] and str(next_hop_node) in intfs[0]: # If the source node is in the second interface and the next hop node is in the first interface
                        src_node_intf = intfs[1] # Set the source node interface to the second interface
                        next_h_node_intf = intfs[0] # Set the next hop node interface to the first interface

                if str(dest_node)+str("-") in link and str(last_hop_node)+str("-") in link: # If the destination node and last hop node are in the link
                    intfs = link.split(":") # Split the link into its component interfaces
                    if str(dest_node) in intfs[0] and str(last_hop_node) in intfs[1]: # If the destination node is in the first interface and the last hop node is in the second interface
                        dest_node_intf = intfs[0] # Set the destination node interface to the first interface
                        last_h_node_intf = intfs[1] # Set the last hop node interface to the second interface
                    elif str(dest_node) in intfs[1] and str(last_hop_node) in intfs[0]: # If the destination node is in the second interface and the last hop node is in the first interface
                        dest_node_intf = intfs[1] # Set the destination node interface to the second interface
                        last_h_node_intf = intfs[0] # Set the last hop node interface to the first interface

            if dest_node_intf != "" and next_h_node_intf != "" and src_node_intf !="": # If the destination node interface, next hop node interface, and source node interface are not empty
                cmd_on_src_node  = "ip route add "+get_network_address(get_node_intf_ip(dest_node_intf, list_of_Intf_IPs))+"/28 via "+get_node_intf_ip(next_h_node_intf, list_of_Intf_IPs)+" dev "+src_node_intf # Set the command on the source node to add the network address of the destination node via the next hop node interface
                print(cmd_on_src_node) # Print the command to the console

            if src_node_intf != "" and last_h_node_intf != "" and dest_node_intf != "": # If the source node interface, last hop node interface, and destination node interface are not empty
                cmd_on_dest_node = "ip route add "+get_network_address(get_node_intf_ip(src_node_intf, list_of_Intf_IPs))+"/28 via "+get_node_intf_ip(last_h_node_intf, list_of_Intf_IPs)+" dev "+dest_node_intf # Set the command on the destination node to add the network address of the source node via the last hop node interface
                print(cmd_on_dest_node) # Print the command to the console







def get_static_route_parameter_optimised(route, links, list_of_Intf_IPs, satellites):
    """
    Get the optimized static route parameters.

    Args:
        route (list): The route to be optimized.
        links (dict): The dictionary of links.
        list_of_Intf_IPs (dict): The dictionary of interface IPs.
        satellites (list): The list of satellites.

    Returns:
        list: The optimized static route parameters in the format src_node, dest_network_ip, next_hop_ip, src_node_intf, dest_node, src_network_ip, last_hop_ip, dest_node_intf.
    """
    parameters = [] # Initialize the parameters list
    current_route = route[0]#; # Get the current route
    src_node        = "" # Initialize the source node, next hop node, destination node, and last hop node variable strings
    next_hop_node   = ""
    dest_node       = ""
    last_hop_node   = ""
    link            = ""

    if len(current_route) > 2: # If the route has more than 2 nodes
        src_node, next_hop_node, dest_node, last_hop_node = current_route[0], current_route[1], current_route[len(current_route)-1], current_route[len(current_route)-2] # Set the source, next hop, destination, and last hop node variables
    elif len(current_route) == 2: # If the route has 2 nodes
        src_node, next_hop_node, dest_node, last_hop_node = current_route[0], current_route[1], current_route[1], current_route[0] # Set the source, next hop, destination, and last hop node variables appropriately

    if src_node != "": # If the source node is not empty
        src_node = "sat"+str(src_node) if src_node < len(satellites) else "gs"+str(src_node%len(satellites)) # Set the source node to the satellite or ground station number
        next_hop_node = "sat"+str(next_hop_node) if next_hop_node < len(satellites) else "gs"+str(next_hop_node%len(satellites)) # Set the next hop node to the satellite or ground station number
        dest_node = "sat"+str(dest_node) if dest_node < len(satellites) else "gs"+str(dest_node%len(satellites)) # Set the destination node to the satellite or ground station number
        last_hop_node = "sat"+str(last_hop_node) if last_hop_node < len(satellites) else "gs"+str(last_hop_node%len(satellites)) # Set the last hop node to the satellite or ground station number

        src_node_intf = "" # Initialize the source node interface, destination node interface, next hop node interface, and last hop node interface variable strings
        dest_node_intf= ""
        next_h_node_intf = ""
        last_h_node_intf = ""

        key1 = str(src_node)+"_"+str(next_hop_node) # Set key1 as the source node and next hop node
        key2 = str(next_hop_node)+"_"+str(src_node) # Set key2 as the next hop node and source node
        if links.get(key1) is not None: # If key1 is in the links dictionary
            link = links[key1][0] # Set the link to the value of key1
        elif links.get(key2) is not None: # If key2 is in the links dictionary
            link = links[key2][0] # Set the link to the value of key2
        if link != "": # If the link is not empty
            intfs = link.split(":") # Split the link into its component interfaces
            if str(src_node) in intfs[0] and str(next_hop_node) in intfs[1]: # If the source node is in the first interface and the next hop node is in the second interface
                src_node_intf = intfs[0] # Set the source node interface to the first interface
                next_h_node_intf = intfs[1] # Set the next hop node interface to the second interface
            elif str(src_node) in intfs[1] and str(next_hop_node) in intfs[0]: # If the source node is in the second interface and the next hop node is in the first interface
                src_node_intf = intfs[1] # Set the source node interface to the second interface
                next_h_node_intf = intfs[0] # Set the next hop node interface to the first interface


        key1 = str(last_hop_node)+"_"+str(dest_node) # Set key1 as the last hop node and destination node
        key2 = str(dest_node)+"_"+str(last_hop_node) # Set key2 as the destination node and last hop node
        if links.get(key1) is not None: # If key1 is in the links dictionary
            link = links[key1][0] # Set the link to the value of key1
        elif links.get(key2) is not None: # If key2 is in the links dictionary
            link = links[key2][0] # Set the link to the value of key2
        if link != "": # If the link is not empty
            intfs = link.split(":") # Split the link into its component interfaces
            if str(dest_node) in intfs[0] and str(last_hop_node) in intfs[1]: # If the destination node is in the first interface and the last hop node is in the second interface
                dest_node_intf = intfs[0] # Set the destination node interface to the first interface
                last_h_node_intf = intfs[1] # Set the last hop node interface to the second interface
            elif str(dest_node) in intfs[1] and str(last_hop_node) in intfs[0]: # If the destination node is in the second interface and the last hop node is in the first interface
                dest_node_intf = intfs[1] # Set the destination node interface to the second interface
                last_h_node_intf = intfs[0] # Set the last hop node interface to the first interface

        if dest_node_intf != "": # If the destination node interface is not empty
            dest_ip_address = list_of_Intf_IPs[str(dest_node_intf)][0] # Set the destination IP address to the value of the destination node interface in the list of interface IPs
            next_hop_ip = list_of_Intf_IPs[str(next_h_node_intf)][0] # Set the next hop IP address to the value of the next hop node interface in the list of interface IPs
            out_interface = src_node_intf # Set the output interface to the source node interface
            # print dest_ip_address
            dest_nw_ip = get_network_address(dest_ip_address.split("/")[0])+"/28" # Set the destination network IP to the network address of the destination IP address
        else:
            print("Error: No link between ", str(last_hop_node), " and ", str(dest_node)) # If destination node interface is empty, print an error message to the console
            exit()

        if src_node_intf != "": # If the source node interface is not empty
            src_ip_address = list_of_Intf_IPs[str(src_node_intf)][0] # Set the source IP address to the value of the source node interface in the list of interface IPs
            last_hop_ip = list_of_Intf_IPs[str(last_h_node_intf)][0] # Set the last hop IP address to the value of the last hop node interface in the list of interface IPs
            out_interface_2 = dest_node_intf # Set the output interface to the destination node interface
            # print src_ip_address
            src_nw_ip = get_network_address(src_ip_address.split("/")[0])+"/28" # Set the source network IP to the network address of the source IP address
        else:
            print("Error: No link between ", str(src_node), " and ", str(next_hop_node)) # If source node interface is empty, print an error message to the console
            exit()

        parameters.append(src_node) # Append the source node, destination network IP, next hop IP, source node interface, destination node, source network IP, last hop IP, and destination node interface to the parameters list
        parameters.append(dest_nw_ip)
        parameters.append(next_hop_ip)
        parameters.append(out_interface)

        parameters.append(dest_node)
        parameters.append(src_nw_ip)
        parameters.append(last_hop_ip)
        parameters.append(out_interface_2)

        # parameters.append(src_node)
        # parameters.append(dest_node)
        # parameters.append(next_hop_node)
        # parameters.append(last_hop_node)
        #
        # parameters.append(src_nw_ip)
        # parameters.append(dest_nw_ip)
        # parameters.append(next_hop_ip)
        # parameters.append(last_hop_ip)
        #
        # parameters.append(out_interface)
        # parameters.append(out_interface_2)

    return parameters # Return the parameters list









def get_static_route_parameter(route, links, list_of_Intf_IPs, satellites):
    """
    Get the static route parameters for a given route.

    Args:
        route (list): The route to get the parameters for.
        links (list): The list of links between nodes.
        list_of_Intf_IPs (list): The list of interface IPs for each node.
        satellites (list): The list of satellite nodes.

    Returns:
        list: The static route parameters in the format: [src_node, dest_node, next_hop_node, last_hop_node, src_network_ip, dest_network_ip, next_hop_ip, last_hop_ip, src_node_intf, dest_node_intf]

    Raises:
        ValueError: If unable to find an interface for any step in the route
    """
    parameters = [] # Initialize the parameters list
    current_route = route[0]#; # Get the current route
    if len(current_route) > 2: # If the route has more than 2 nodes
        src_node, next_hop_node, dest_node, last_hop_node = current_route[0], current_route[1], current_route[len(current_route)-1], current_route[len(current_route)-2] # Set the source, next hop, destination, and last hop node variables

        src_node = "sat"+str(src_node) if src_node < len(satellites) else "gs"+str(src_node%len(satellites)) # Set the source node to the satellite or ground station number
        next_hop_node = "sat"+str(next_hop_node) if next_hop_node < len(satellites) else "gs"+str(next_hop_node%len(satellites)) # Set the next hop node to the satellite or ground station number
        dest_node = "sat"+str(dest_node) if dest_node < len(satellites) else "gs"+str(dest_node%len(satellites)) # Set the destination node to the satellite or ground station number
        last_hop_node = "sat"+str(last_hop_node) if last_hop_node < len(satellites) else "gs"+str(last_hop_node%len(satellites)) # Set the last hop node to the satellite or ground station number

        src_node_intf = "" # Initialize the source node interface, destination node interface, next hop node interface, and last hop node interface variable strings
        dest_node_intf= "" 
        next_h_node_intf = ""
        last_h_node_intf = ""

        for link in links:  # For each link in the list of links
            if str(src_node)+str("-") in link and str(next_hop_node)+str("-") in link: # If the source node and next hop node are in the link
                intfs = link.split(":") # Split the link into its component interfaces
                if str(src_node) in intfs[0] and str(next_hop_node) in intfs[1]: # If the source node is in the first interface and the next hop node is in the second interface
                    src_node_intf = intfs[0] # Set the source node interface to the first interface
                    next_h_node_intf = intfs[1] # Set the next hop node interface to the second interface
                elif str(src_node) in intfs[1] and str(next_hop_node) in intfs[0]: # If the source node is in the second interface and the next hop node is in the first interface
                    src_node_intf = intfs[1] # Set the source node interface to the second interface
                    next_h_node_intf = intfs[0] # Set the next hop node interface to the first interface

            if str(dest_node)+str("-") in link and str(last_hop_node)+str("-") in link:  # If the destination node and last hop node are in the link
                intfs = link.split(":") # Split the link into its component interfaces
                if str(dest_node) in intfs[0] and str(last_hop_node) in intfs[1]: # If the destination node is in the first interface and the last hop node is in the second interface
                    dest_node_intf = intfs[0] # Set the destination node interface to the first interface
                    last_h_node_intf = intfs[1] # Set the last hop node interface to the second interface
                elif str(dest_node) in intfs[1] and str(last_hop_node) in intfs[0]: # If the destination node is in the second interface and the last hop node is in the first interface
                    dest_node_intf = intfs[1] # Set the destination node interface to the second interface
                    last_h_node_intf = intfs[0] # Set the last hop node interface to the first interface

        if dest_node_intf != "": # If the destination node interface is not empty
            dest_nw_ip = get_network_address(get_node_intf_ip(dest_node_intf, list_of_Intf_IPs).split("/")[0])+"/28" # Set the destination network IP to the network address of the destination node interface
            next_hop_ip = get_node_intf_ip(next_h_node_intf, list_of_Intf_IPs) # Set the next hop IP to the IP address of the next hop node interface
            out_interface = src_node_intf # Set the output interface to the source node interface
        else: # If the destination node interface is empty
            print("Error: No link between ", str(last_hop_node), " and ", str(dest_node)) # Print an error message to the console
            exit()

        if src_node_intf != "": # If the source node interface is not empty
            src_nw_ip = get_network_address(get_node_intf_ip(src_node_intf, list_of_Intf_IPs).split("/")[0])+"/28" # Set the source network IP to the network address of the source node interface
            last_hop_ip = get_node_intf_ip(last_h_node_intf, list_of_Intf_IPs) # Set the last hop IP to the IP address of the last hop node interface
            out_interface_2 = dest_node_intf # Set the output interface to the destination node interface
        else: # If the source node interface is empty
            print("Error: No link between ", str(src_node), " and ", str(next_hop_node)) # Print an error message to the console
            exit()

        parameters.append(src_node) # Append the source node, destination node, next hop node, last hop node, source network IP, destination network IP, next hop IP, last hop IP, source node interface, and destination node interface to the parameters list
        parameters.append(dest_node)
        parameters.append(next_hop_node)
        parameters.append(last_hop_node)

        parameters.append(src_nw_ip)
        parameters.append(dest_nw_ip)
        parameters.append(next_hop_ip)
        parameters.append(last_hop_ip)

        parameters.append(out_interface)
        parameters.append(out_interface_2)

    if len(current_route) == 2: # If the route has 2 nodes
        src_node, next_hop_node, dest_node, last_hop_node = current_route[0], current_route[1], current_route[1], current_route[0] # Set the source, next hop, destination, and last hop node variables appropriately

        src_node = "sat"+str(src_node) if src_node < len(satellites) else "gs"+str(src_node%len(satellites)) # Set the source node to the satellite or ground station number
        next_hop_node = "sat"+str(next_hop_node) if next_hop_node < len(satellites) else "gs"+str(next_hop_node%len(satellites)) # Set the next hop node to the satellite or ground station number
        dest_node = "sat"+str(dest_node) if dest_node < len(satellites) else "gs"+str(dest_node%len(satellites)) # Set the destination node to the satellite or ground station number
        last_hop_node = "sat"+str(last_hop_node) if last_hop_node < len(satellites) else "gs"+str(last_hop_node%len(satellites)) # Set the last hop node to the satellite or ground station number

        src_node_intf = "" # Initialize the source node interface, destination node interface, next hop node interface, and last hop node interface variable strings
        dest_node_intf= ""
        next_h_node_intf = ""
        last_h_node_intf = ""

        for link in links: # For each link in the list of links
            if str(src_node)+str("-") in link and str(next_hop_node)+str("-") in link: # If the source node and next hop node are in the link
                intfs = link.split(":") # Split the link into its component interfaces
                if str(src_node) in intfs[0] and str(next_hop_node) in intfs[1]: # If the source node is in the first interface and the next hop node is in the second interface
                    src_node_intf = intfs[0] # Set the source node interface to the first interface
                    next_h_node_intf = intfs[1] # Set the next hop node interface to the second interface
                elif str(src_node) in intfs[1] and str(next_hop_node) in intfs[0]: # If the source node is in the second interface and the next hop node is in the first interface
                    src_node_intf = intfs[1] # Set the source node interface to the second interface
                    next_h_node_intf = intfs[0] # Set the next hop node interface to the first interface

            if str(dest_node)+str("-") in link and str(last_hop_node)+str("-") in link: # If the destination node and last hop node are in the link
                intfs = link.split(":") # Split the link into its component interfaces
                if str(dest_node) in intfs[0] and str(last_hop_node) in intfs[1]: # If the destination node is in the first interface and the last hop node is in the second interface
                    dest_node_intf = intfs[0] # Set the destination node interface to the first interface
                    last_h_node_intf = intfs[1] # Set the last hop node interface to the second interface
                elif str(dest_node) in intfs[1] and str(last_hop_node) in intfs[0]: # If the destination node is in the second interface and the last hop node is in the first interface
                    dest_node_intf = intfs[1] # Set the destination node interface to the second interface
                    last_h_node_intf = intfs[0] # Set the last hop node interface to the first interface

        if dest_node_intf != "": # If the destination node interface is not empty
            dest_nw_ip = get_network_address(get_node_intf_ip(dest_node_intf, list_of_Intf_IPs).split("/")[0])+"/28" # Set the destination network IP to the network address of the destination node interface
            next_hop_ip = get_node_intf_ip(next_h_node_intf, list_of_Intf_IPs) # Set the next hop IP to the IP address of the next hop node interface
            out_interface = src_node_intf # Set the output interface to the source node interface
        else: # If the destination node interface is empty
            print("Error: No link between ", str(last_hop_node), " and ", str(dest_node)) # Print an error message to the console
            exit()

        if src_node_intf != "": # If the source node interface is not empty
            src_nw_ip = get_network_address(get_node_intf_ip(src_node_intf, list_of_Intf_IPs).split("/")[0])+"/28" # Set the source network IP to the network address of the source node interface
            last_hop_ip = get_node_intf_ip(last_h_node_intf, list_of_Intf_IPs) # Set the last hop IP to the IP address of the last hop node interface
            out_interface_2 = dest_node_intf # Set the output interface to the destination node interface
        else: # If the source node interface is empty
            print("Error: No link between ", str(src_node), " and ", str(next_hop_node)) # Print an error message to the console
            exit()

        parameters.append(src_node) # Append the source node, destination node, next hop node, last hop node, source network IP, destination network IP, next hop IP, last hop IP, source node interface, and destination node interface to the parameters list
        parameters.append(dest_node)
        parameters.append(next_hop_node)
        parameters.append(last_hop_node)

        parameters.append(src_nw_ip)
        parameters.append(dest_nw_ip)
        parameters.append(next_hop_ip)
        parameters.append(last_hop_ip)

        parameters.append(out_interface)
        parameters.append(out_interface_2)

    return parameters # Return the parameters list








def find_route_between_src_dest(src_sat, dest_sat, constellation_routes):
    """
    Finds the route between the source satellite and the destination satellite in the given constellation routes.

    Args:
        src_sat (int): The source satellite ID.
        dest_sat (int): The destination satellite ID.
        constellation_routes (dict): A dictionary containing the routes for each satellite in the constellation.

    Returns:
        list or int: The route between the source and destination satellites as a list of satellite IDs, or -1 if no route is found.
    """
    if src_sat < dest_sat: # If the source satellite number is less than the destination satellite number
        for route in constellation_routes[src_sat]: # For each route in the list of routes for the source satellite
            # print route
            current_route = route[0][:] # Makes shallow copy of route[0]
            if "sat"+str(src_sat) == "sat"+str(current_route[0]) and "sat"+str(dest_sat) == "sat"+str(current_route[len(current_route)-1]): # If the source satellite is the first node in the route and the destination satellite is the last node in the route
                return current_route # Return the current route

    # This is added because in constellation_routes and in initial_routes, we only store the route from x to y but not the route for y to x.
    # We do that to minimize the calcuations. Therefore, sometimes if src_sat > dest_sat, we need to search in the opposite direction
    if src_sat > dest_sat: # If the source satellite number is greater than the destination satellite number
        for route in constellation_routes[dest_sat]: # For each route in the list of routes for the destination satellite
            current_route = route[0][:] # Makes shallow copy of route[0]
            if "sat"+str(src_sat) == "sat"+str(current_route[len(current_route)-1]) and "sat"+str(dest_sat) == "sat"+str(current_route[0]): # If the source satellite is the last node in the route and the destination satellite is the first node in the route
                retr_route = route[0][:] # Makes shallow copy of route[0]
                retr_route.reverse() # Reverse the route
                return retr_route # Return the reversed route

    return -1 # Return -1 if no route is found









def get_gs_ip(list_of_Intf_IPs, gs):
    """
    Get the IP address associated with a given ground station.

    Args:
        list_of_Intf_IPs (list): A list of dictionaries containing interface and IP information.
        gs (str): The name of the ground station.

    Returns:
        str: The IP address associated with the ground station, or -1 if not found.
    """
    for pair in list_of_Intf_IPs: # For each pair in the list of interface IPs
        if gs in pair["Interface"]: # If the ground station is in the interface
            return pair["IP"] # Return the IP address associated with the ground station

    return -1 # Return -1 if the ground station is not found










def gs_routing(data_path, gs_statellite_pair, links_updated, num_of_satellites, satellites_by_index, list_of_Intf_IPs, constellation_routes, main_configurations, border_gateway):
    """
    Generates BASH commands for ground station routing and writes them to disk.

    Args:
        data_path (str): The path to the data.
        gs_statellite_pair (list): List of ground station-satellite pairs.
        links_updated (list): List of updated links.
        num_of_satellites (int): Number of satellites.
        satellites_by_index (dict): Dictionary of satellites by index.
        list_of_Intf_IPs (dict): Dictionary of interface IPs.
        constellation_routes (list): List of constellation routes.
        main_configurations (dict): Main configurations.
        border_gateway (str): Border gateway.

    Returns:
        None
    """
    update_gsl_routing_cmd = [] # Initialize the update GSL routing command list
    for gs_sat in gs_statellite_pair: # For each ground station-satellite pair
        gs_number = int(gs_sat[1])%num_of_satellites # Get the ground station number

        if "gs"+str(gs_number) != border_gateway: # If the ground station is not the border gateway (What is the 'border_gateway'? Is this a ground station that connects to another network's ground station?)
            key = str("gs"+str(gs_number)+"-eth0") # Set the key to the ground station number and eth0
        else: # If the ground station is the border gateway
            key = str("gs"+str(gs_number)+"-eth1") # Set the key to the ground station number and eth1

        # print key
        if list_of_Intf_IPs.get(key) is None: # If the key is not in the list of interface IPs
            print("error -- no ip for this ground station", gs_number) # Print an error message to the console
            return

        gs_ip = list_of_Intf_IPs[key][0].split("/")[0] # Set the ground station IP to the value of the key in the list of interface IPs
        gs_network_address = get_network_address(gs_ip) # Set the ground station network address to the network address of the ground station IP
        thread_list = [] # Initialize the thread list
        for i in range(num_of_satellites): # For each satellite
            route_to_sat_GW = find_route_between_src_dest(i, gs_sat[0], constellation_routes) # Find the route to the satellite gateway
            # if key == "gs1-eth1":
            #     print route_to_sat_GW
            if route_to_sat_GW != -1: # If the route to the satellite gateway is not -1
                parameters = get_static_route_parameter_optimised([route_to_sat_GW], links_updated, list_of_Intf_IPs, satellites_by_index)#; # Get the static route parameters
                # if key == "gs1-eth1":
                #     print parameters
                if len(parameters) > 0: # If the parameters list is not empty
                    update_gsl_routing_cmd.append("sat"+str(i)+",ip route add "+str(gs_network_address)+"/28 via "+str(parameters[2][:-3])+" dev "+str(parameters[3])) # Append the GSL ip routing add command to the update GSL routing command list
                    if main_configurations["constellation"]["routing"]["interDomain_routing"] == 1 and "gs"+str(gs_number) == main_configurations["constellation"]["routing"]["border_gateway"]: # If inter-domain routing is enabled and the ground station is the border gateway
                        update_gsl_routing_cmd.append("sat"+str(i)+",ip route add "+str(main_configurations["constellation"]["routing"]["other_constellation_ip_range"])+"/20 via "+str(parameters[2][:-3])+" dev "+str(parameters[3])) # Append the GSL ip routing add command to the update GSL routing command list
                else: # If the parameters list is empty
                    print("-----> ", route_to_sat_GW) # Print the route to the satellite gateway to the console

    if len(update_gsl_routing_cmd) > 0: # If the update GSL routing command list is not empty
        for update in update_gsl_routing_cmd: # For each update in the update GSL routing command list
            update_routes = update.split(",") # Split the update into its component parts (#sat1539,ip route del 10.2.2.112 via 10.2.6.130/28 dev sat1539-eth4)
            file = open(data_path+"/cmd_files/"+update_routes[0]+"_routes.sh", 'a') # Open the file for appending
            file.writelines(update_routes[1].strip()+" & \n") # Write the update to the file and append & to run the command in the background
            file.close() # Close the file










def gs_routing_parallel(data_path, gs_statellite_pair, links_updated, num_of_satellites, satellites_by_index, list_of_Intf_IPs, constellation_routes, num_of_threads):
    """
    Perform parallel ground station routing.

    Args:
        data_path (str): The path to the data.
        gs_statellite_pair (list): List of ground station-satellite pairs.
        links_updated (bool): Flag indicating if links have been updated.
        num_of_satellites (int): Number of satellites.
        satellites_by_index (dict): Dictionary mapping satellite index to satellite object.
        list_of_Intf_IPs (list): List of interface IP addresses.
        constellation_routes (dict): Dictionary of constellation routes.
        num_of_threads (int): Number of threads to use for parallel execution.
    """
    gs_routing_args = [] # Initialize the ground station routing arguments list
    thread_list = [] # Initialize the thread list

    num_thread = num_of_threads#; # Set the number of threads to the number of threads to use for parallel execution

    for gs_sat in gs_statellite_pair: # For each ground station-satellite pair
        thread = threading.Thread(target=gs_routing_worker, args=(data_path, gs_sat, links_updated, num_of_satellites, satellites_by_index, list_of_Intf_IPs, constellation_routes)) # Create a new thread for the ground station routing worker
        thread_list.append(thread) # Append the thread to the thread list

    for thread in thread_list: # For each thread in the thread list
        thread.start() # Start the thread
    for thread in thread_list: # For each thread in the thread list
        thread.join() # Join the thread











def lightweight_routing(data_path, route_changes, links_updated, num_of_satellites, satellites_by_index, list_of_Intf_IPs, constellation_routes, t_time, border_gateway):
    """
    Generates BASH commands for node routing in the given time and writes them to disk.

    Args:
        data_path (str): The path to the data directory.
        route_changes (list): List of route changes.
        links_updated (list): List of updated links.
        num_of_satellites (int): Number of satellites.
        satellites_by_index (dict): Dictionary mapping satellite index to satellite name.
        list_of_Intf_IPs (dict): Dictionary mapping interface name to IP addresses.
        constellation_routes (dict): Dictionary of current topology routes.
        t_time (datetime): The current time.
        border_gateway (str): The border gateway name.

    Returns:
        None
    """
    update_gsl_routing_cmd = [] # Initialize the update GSL routing command list
    isl_ch = 0 # Initialize the number of ISL changes
    gsl_ch = 0 # Initialize the number of GSL changes
    allchanges_log = open(data_path+"/allchanges_log_"+str(t_time.utc_strftime())+"_.txt", "a") # Open the all changes log file for appending
    for change in route_changes: # For each change in the route changes
        start = round(time.time()*1000) # Set the start time to the current time
        # print "the updated route --> ", change
        allchanges_log.write(str(change[0])+","+str(change[1])+","+str(change[2])+","+str(change[3])+"\n") # Write the change to the all changes log file

        # #####
        # If there is a change the GSL links
        # #####
        if change[0] < num_of_satellites and change[1] >= num_of_satellites: # If the first element of the change is a satellite number and the second element of the change is a ground station number
            # 1. Get the network address of the changing ground station
            gs_number = int(change[1])%num_of_satellites # Get the ground station number
            if "gs"+str(gs_number) != border_gateway: # If the ground station is not the border gateway
                key = str("gs"+str(gs_number)+"-eth0") # Set the key to the ground station number and eth0
            else: # If the ground station is the border gateway
                key = str("gs"+str(gs_number)+"-eth1") # Set the key to the ground station number and eth1

            if list_of_Intf_IPs.get(key) is None: # If the key is not in the list of interface IPs
                print("error -- no ip for this ground station", gs_number) # Print an error message to the console
                return

            gs_ip = list_of_Intf_IPs[key][0].split("/")[0] # Set the ground station IP to the value of the key in the list of interface IPs
            gs_network_address = get_network_address(gs_ip) # Set the ground station network address to the network address of the ground station IP

            ## Original comments ##
            # gs_ip = get_gs_ip(list_of_Intf_IPs, "gs"+str(gs_number)+"-eth0").split("/")[0]
            # gs_network_address = get_network_address(gs_ip)

            # 2. Find the route to the current gateway of that ground station. The current gateway is at change[0].
            #    The ground station is at change[1]. constellation_routes has all the current topology routes.
            # 3. If there is a route from satellite (i) to the gateway of that ground station:
            #    3.1. Check the details of that route, specifically,
            #           (1) The IP address of the next hop to reach that gateway --> parameters[2]
            #           (2) The output interface to reach that gateway           --> parameters[3]
            #    3.2. Based on the type of the route update: Is it route addition or deletion  (change[3])
            #         We generate the IP route command.
            ## End Original comments ##
            thread_list = [] # Initialize the thread list
            for i in range(num_of_satellites): # For each satellite

                route_to_sat_GW = find_route_between_src_dest(i, change[0], constellation_routes) # Find the route to the satellite gateway

                if route_to_sat_GW != -1: # If the route to the satellite gateway is not -1

                    parameters = get_static_route_parameter_optimised([route_to_sat_GW], links_updated, list_of_Intf_IPs, satellites_by_index)#; # Get the static route parameters

                    if len(parameters) > 0: # If the parameters list is not empty
                        if change[2] == 0 and change[3] == 1: # If the change[2] is type 0 and the route modification is 1
                            update_gsl_routing_cmd.append("sat"+str(i)+",ip route add "+str(gs_network_address)+"/28 via "+str(parameters[2])+" dev "+str(parameters[3])) # Append the GSL ip routing add command to the update GSL routing command list
                        elif change[2] == 1 and change[3] == 0: # If the change[2] is type 1 and the route modification is 0
                            update_gsl_routing_cmd.append("sat"+str(i)+",ip route del "+str(gs_network_address)+"/28 via "+str(parameters[2])+" dev "+str(parameters[3])) # Append the GSL ip routing delete command to the update GSL routing command list
                else: # If the route to the satellite gateway is -1
                    if i!= int(change[0]): # If i is not equal to change[0] (current gateway?)
                        print("Error: cannot find the route between sat", i, " and sat", change[0]) # Print an error message to the console
            gsl_ch += 1 # Increment the number of GSL changes

        # #####
        # If there is a change the ISL links
        # #####
        elif change[0] < num_of_satellites and change[1] < num_of_satellites: # If the first two elements of change are satellite numbers
            isl_ch += 1 # Increment the number of ISL changes

        end = round(time.time()*1000) # Set the end time to the current time
        # print "one iteration of change takes --- ", end-start, "ms "

    if len(update_gsl_routing_cmd) > 0: # If the update GSL routing command list is not empty
        start = round(time.time()*1000) # Set the start time to the current time
        if os.path.isdir(data_path+"/routes_updates_"+str(t_time.utc_strftime())) == False: # If the routes updates directory does not exist
            os.mkdir(data_path+"/routes_updates_"+str(t_time.utc_strftime())) # Create the routes updates directory

        for f in os.listdir(data_path+"/routes_updates_"+str(t_time.utc_strftime())): # For each file in the routes updates directory
            os.remove(os.path.join(data_path+"/routes_updates_"+str(t_time.utc_strftime()), f)) # Remove the file

        updates_log = open(data_path+"/routes_updates_"+str(t_time.utc_strftime())+"/routing_updates_"+str(t_time.utc_strftime())+"_.txt", "w") # Open the updates log file for writing
        for update in update_gsl_routing_cmd: # For each update in the update GSL routing command list
            updates_log.write(str(update) + "\n") # Write the update to the updates log file
            update_routes = update.split(",") # Split the update into its component parts (#sat1539,ip route del 10.2.2.112 via 10.2.6.130/28 dev sat1539-eth4)
            file = open(data_path+"/routes_updates_"+str(t_time.utc_strftime())+"/"+update_routes[0]+"_routes.sh", 'a') # Open the file for appending
            file.writelines(update_routes[1].strip()+"\n") # Write the update to the file
            file.close() # Close the file

        cmd_path = data_path+"/routes_updates_"+str(t_time.utc_strftime()) # Set the command path to the routes updates directory
        cmdd =  'chmod +x "'+cmd_path+'"'+'/*.sh' # Set the command to change the permissions of the files in the command path
        os.system(cmdd) # Run the command
        end = round(time.time()*1000) # Set the end time to the current time
        # print "writing the routes update to file takes --- ", end-start, "ms "

        updates_log.close() # Close the updates log file
        allchanges_log.close() # Close the all changes log file

    # print gsl_ch, isl_ch





def check_changes_in_topology(last, new):
    """
    Check for changes in the topology between two states.

    Args:
        last (list): The last known topology state.
        new (list): The new topology state.

    Returns:
        list: A list of tuples representing the changes in the topology.
              Each tuple contains the row index, column index, previous value,
              and new value of the changed element.
    """
    changes = [] # Initialize the changes list
    for i in range(len(new)): # For each element in the new topology
        for j in range(len(new[i])): # For each sub-element in the selected element
            if new[i][j] != last[i][j]: # If the new sub-element is not equal to the last sub-element
                changes.append((i, j, last[i][j], new[i][j])) # Append the row index, column index, previous value, and new value to the changes list
                # print i,j, last[i][j], new[i][j]

    return changes # Return the changes list

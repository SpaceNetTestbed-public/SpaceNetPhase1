''''

SpaceNet Utils/

AUTHOR:         Mohamed M. Kassem, Ph.D.
                University of Surrey

EDITOR:         Bruce Barbour
                Virginia Tech

DESCRIPTION:    This Python script supplies the primary utility functions for the SpaceNet simulator.


CONTENTS:       FILE SYSTEM/ (STARTS AT 80)
                    get_recent_TLEs_using_timestamp(path, timestamp, constellation)
                    get_recent_TLEs_using_datetime(path, datetime, constellation)
                    read_IProute_files_thread(routes, initial_routes)
                    save_topology(connectivity_matrix, links_charateristics, main_configurations, timestamp)
                    save_routes(routes, main_configurations, timestamp)
                    save_optimal_path(optimal_path, main_configurations, timestap)

                PARSING/ (STARTS AT 275)
                    parse_config_file_yml(filepath, filename)                     
                    parse_resiliency_experiment_parameters(main_configurations)
                    parse_config_file(filepath, filename)
                    parse_connectivity_matrix_n_charateristics(time_utc, conn_mat_size, topology_path)
                    parse_topology_routes(topology_route_path, num_of_satellites, time_utc)
                    parse_resiliency_experiment_parameters(main_configurations)
                
                SATELLITE/ (STARTS AT 555)
                    arrange_satellites(path, orbital_data, satellites_by_name, main_configurations, simulation_start_time, satellites_by_index, tle_timestamp)
                    reload_tles(path_of_recent_TLE, main_configurations)
                    get_gs_sat_pairs(connectivity_matrix, num_of_satellites)
                    get_sats_by_index(filename)
                    get_sats_by_name(filename)

                CONVERSION/ (STARTS AT 812)
                    convert_time_utc_to_unix(time_utc)
                    convert_time_utc_to_ymdhms(time_utc)

                NETWORK/ (STARTS AT 870)
                    check_changes_in_link_charateristics(old_links_charateristics, new_links_charateristics)
                    merge_link_link_charateristics(latency_changes, capacity_changes)
                    apply_topology_updates_to_mininet(path, net, topology_changes, num_of_satellites, time_utc_inc)
                    apply_link_updates_to_mininet(net, latency_changes, capacity_changes, num_of_satellites, time_utc_inc)

                PERFORMANCE TESTING/ (STARTS AT 1103)
                    iperf_app(path, net, main_configurations, list_of_Intf_IPs)
                    ping_app(path, net, main_configurations, list_of_Intf_IPs)
                    run_application(path, net, main_configurations, list_of_Intf_IPs)

                RESILIENCY APP/ (STARTS AT 1222)  
                    check_time_to_deploy_RE(resiliency_satellite_timestamp, year, month, day, hour, minute, seconds)

'''


# =================================================================================== #
# ------------------------------- IMPORT PACKAGES ----------------------------------- #
# =================================================================================== #

import threading
import os
from datetime import *
import time
import re
from mininet.net import Mininet
from mininet.node import Node, OVSKernelSwitch, Controller, RemoteController
from mininet.cli import CLI
from mininet.link import TCLink
from mininet.link import *
from mininet.topo import Topo
from mininet.log import setLogLevel, info
from mininet.node import OVSController
import yaml
import sys
sys.path.append("../")
from mobility.read_live_tles import *
from utils.file_utils import *

# =================================================================================== #
# ---------------------------------- FILE SYSTEM ------------------------------------ #
# =================================================================================== #
def get_recent_TLEs_using_timestamp(
                                        path                : str, 
                                        timestamp           : float, 
                                        constellation       : str
                                   ) -> str:
    """
    Find the most recent TLE path that corresponds to the provided unix time.

    Args:
        path (str):             Directory of TLE files
        timestamp (float):      Unix time
        constellation (str):    Type of constellation, typically defined as 'starlink'

    Returns:
        str:                    File path of TLE
    """

    # Initialize variables
    recent_file = ""
    timestamp_diff = 999999999999
    directory = path+constellation+'_tles'

    # Search for TLE in directory
    for filename in os.listdir(directory):
        f = os.path.join(directory, filename)
        if os.path.isfile(f):
            file_timestamp = int(filename.split("_")[1])
            if abs(int(timestamp-file_timestamp)) < timestamp_diff and int(timestamp-file_timestamp) <= 86400:
                timestamp_diff = abs(int(timestamp-file_timestamp))
                recent_file = f

    # Return file path
    #print(timestamp_diff, file_timesamp, timestamp)
    return recent_file


def get_recent_TLEs_using_datetime(
                                    path                : str,  
                                    datetime            : str, 
                                    constellation       : str
                                  ) -> str:
    """
    Find the most recent TLE path that corresponds to the provided UTC datetime.

    Args:
        path (str):             Directory of TLE files
        datetime (str):         UTC datetime
        constellation (str):    Type of constellation, typically defined as 'starlink'

    Returns:
        str:                    File path of TLE
    """

    # Create timescale object
    ts = load.timescale()

    # Sort the time components
    y,mon,d,h,min,s = datetime.split(",")[0], datetime.split(",")[1], datetime.split(",")[2], datetime.split(",")[3], datetime.split(",")[4], datetime.split(",")[5]

    # Convert from UTC to unix
    time_utc = ts.utc(int(y), int(mon), int(d), int(h), int(min), float(s))
    time_timestamp = convert_time_utc_to_unix(time_utc)

    # Return file path using converted unix
    return get_recent_TLEs_using_timestamp(path, time_timestamp, constellation)


def read_IProute_files_thread(
                                routes              : list, 
                                initial_routes      : list
                             ) -> list:
    """
    Appends new routes to the current list of routes.

    Args:
        routes (list):          New list of routes to be appended
        initial_routes (list):  Current list of routes

    Returns:
        list:                   Appended route list
    """

    # Iterate over new list and append each new route to current list
    for route in routes:
        route_new = []
        route_list = re.split(", | |\n", route)
        route = [int(r) for r in route_list if r.strip()]
        route_new.append(route)
        initial_routes.append(route_new)


def save_topology(
                    connectivity_matrix         : list, 
                    links_charateristics        : dict, 
                    operator_name               : str, 
                    timestamp                   : int,
                    connectivity_matrix_path    : str
                 ):
    """
    Saves the link characteristics (latency and bandwidth) for each sat/gs pair in the topology.

    Args:
        connectivity_matrix (list):     Two-dimensional matrix list of sat/gs connections, where each row index (i) corresponds 
                                        a single sat/gs in the sorted list and each column index (j) are each of the sat/gs in the sorted list. 
                                        Any elements of the indexpair (i, j) that are 1 is a connected pair. Discludes loopback.
        link_characteristics (dict):    Dictionary of 2 two-dimensional matrices for computed latency and bandwidth, corresponding
                                        to connected pairs of the indices (i, j), i.e. where element == 1
        operator_name (str):            Constellation/operator name
        timestamp (int):                Unix time
        connectivity_matrix_path (str): Path to output the connectivity matrix files                       

    Returns:
        Saves the topology as a .txt file.
    """

    # Initialize list
    existing_links = []

    # Generate a new file
    file_path = connectivity_matrix_path+operator_name+"/"
    check_create_path(file_path)
    file_name = "topology_"+timestamp+".txt"
    #f = open(connectivity_matrix_path+operator_name+"/topology_"+timestamp+".txt", "a")
    f = open(file_path + file_name, "a")

    # Iterate over the connectivity matrix list
    for i in range(len(connectivity_matrix)):
        for j in range(len(connectivity_matrix[i])):
            if connectivity_matrix[i][j] == 1:
                if i!=j and (i, j) not in existing_links:
                   write_this = str(i)+","+str(j)+","+str(round(links_charateristics["latency_matrix"][i][j],2))+","+str(round(links_charateristics["throughput_matrix"][i][j],2))+"\n"
                   #write_this = str(i)+","+str(j)+","+str(round(links_charateristics["latency_matrix"][i][j],2))+","+str(round(links_charateristics["throughput_matrix"][i][j],2))+","+str(round(links_charateristics["distance_matrix"][i][j],2))+"\n"
                   f.write(write_this)
                   existing_links.append((i, j))
    
    # Close file to minimize memory leaks
    f.close()



def extract_connectivity(  
                            topology_path               : str, 
                            conn_mat_size               : int  
                        ):
    # Extracting connectivity matrix from topology files already generated
    connectivity_matrix = [[0 for _ in range(conn_mat_size)] for r in range(conn_mat_size)]
    latency_matrix = [[0.0 for _ in range(conn_mat_size)] for _ in range(conn_mat_size)]
    throughput_matrix = [[0.0 for _ in range(conn_mat_size)] for _ in range(conn_mat_size)]
    distance_matrix = [[0.0 for _ in range(conn_mat_size)] for _ in range(conn_mat_size)]
    with open(topology_path, 'r') as links:
        for link_data in links:
            data  = link_data.split(",")
            d = data[3].split("\n")
            data[3] = d[0]
            connectivity_matrix[int(data[0])][int(data[1])] = 1
            latency_matrix[int(data[0])][int(data[1])] = float(data[2])
            throughput_matrix[int(data[0])][int(data[1])] = data[3]
            try:
                distance_matrix[int(data[0])][int(data[1])] = data[4]
            except:
                pass

    link_characteristics = {"latency_matrix": latency_matrix, "throughput_matrix": throughput_matrix, "distance_matrix": distance_matrix}

    return connectivity_matrix, link_characteristics
            


def save_routes(
                    routes                  : list, 
                    operator_name           : str, 
                    timestamp               : int,
                    routing_file_path       : str
               ):
    """
    Saves all possible optimal routes from the Bellman-Ford (BF) algorithm for each satellite and ground 
    station node in the topology.

    Args:
        routes (dict):                  Dictionary of all possible routes determined by BF algorithm starting 
                                        from any satellite/ground to any sequence of connected 
                                        satellite/ground nodes
        operator_name (str):            Constellation/operator name
        timestamp (int):                Unix time
        routing_file_path (str):        Path to output the routing files

    Returns:
        Saves the routes as a .txt file.
    """

    # Generate a new file
    file_path = routing_file_path+operator_name+"/"
    check_create_path(file_path)
    file_name = "routes_"+timestamp+".txt"
    #routes_log = open(routing_file_path+operator_name+"/routes_"+timestamp+".txt", "a")
    routes_log = open(file_path + file_name, "a")

    # Iterate over the routes list
    for _, route in routes.items():
        current_route = ', '.join(map(str, route))
        routes_log.write(current_route + "\n")
    
    # Close file to minimize memory leaks
    routes_log.close()


def save_optimal_path(
                        optimal_path            : list, 
                        timestamp               : int,
                        operator_name           : str, 
                        optimal_file_path       : str
                     ):
    """
    Saves a single optimal route determined by BF algorithm between a source and destination node.

    Args:
        optimal_path (list):            List of the optimal route between source/destination nodes
        timestamp (list):               Unix time as a list
        operator_name (str):            Constellation/operator name
        optimal_file_path (str):        Path to output optimal path files

    Returns:
        Saves the route as a .txt file.
    """

    # Generate a new file
    file_path = optimal_file_path+operator_name+"/"
    check_create_path(file_path)
    file_name = "best_path_"+("_".join(timestamp))+".txt"
    #optimal_log = open(optimal_file_path+operator_name+"/best_path_"+("_".join(timestamp[:3]))+".txt", "a")
    optimal_log = open(file_path + file_name, "a")
    
    # Iterate over the optimal path list
    if optimal_path:
        optimal_log.write("(" + ("_".join(timestamp)) + "): " + str(optimal_path)[1:-1] + "\n")
    else:
        optimal_log.write("(" + ("_".join(timestamp)) + "): " + "Unreachable\n")

    # Close file to minimize memory leaks
    optimal_log.close()

def save_optimal_weights(
                        optimal_weights            : list, 
                        timestamp               : int,
                        operator_name           : str, 
                        optimal_w_path          : str
                     ):
    
    # Generate a new file
    file_path = optimal_w_path+operator_name+"/"
    check_create_path(file_path)
    file_name = "optimal_weights_"+("_".join(timestamp[0:3]))+".txt"
    #optimal_log = open(optimal_file_path+operator_name+"/best_path_"+("_".join(timestamp[:3]))+".txt", "a")
    optimal_log = open(file_path + file_name, "a")
    
    # Iterate over the optimal path list
    if optimal_weights:
        optimal_log.write("(" + ("_".join(timestamp)) + "): " + str(optimal_weights) + "\n")
    else:
        optimal_log.write("(" + ("_".join(timestamp)) + "): " + "Unreachable\n")

    # Close file to minimize memory leaks
    optimal_log.close()

"""
def update_node_index(t2t_dict, node_index_file_path, timestamp, operator_name):
    # Open existing file for appending
    file_path = node_index_file_path+operator_name+"/"
    file_name = "nodeindex_"+timestamp+".txt"
    # Open existing file and find the highest used index value
    nodeindex_log_r = open(file_path + file_name, "r")
    max_index = 0
    for line in nodeindex_log_r:
        index = int(line.split(":")[0])
        if index > max_index:
            max_index = index
    nodeindex_log_r.close()
    next_index = max_index + 1
    endpoint_gid_list = []
    nodeindex_log_a = open(node_index_file_path+operator_name+"/nodeindex_"+timestamp+".txt", "a")
    for key in t2t_dict.keys():
        if 'type' in t2t_dict[key] and t2t_dict[key]['type'] == 'endpoint' and 'gid' in t2t_dict[key]:
            endpoint_gid_list.append(t2t_dict[key]['gid'])
            
    endpoint_gid_list = list(set(endpoint_gid_list)) # Ensure no duplicates
    endpoint_gid_list.sort() # Sort the list so that the alias is consistent
    for endpoint_gid in endpoint_gid_list:
        alias_prefix = "IE-"
        nodeindex_log_a.write(str(next_index)+":"+alias_prefix+str(endpoint_gid)+"\n")
        next_index += 1
    nodeindex_log_a.close() # Close file to minimize memory leaks
"""
def save_node_index_and_terrestrial_info(
                                            satellites_by_index     : dict, 
                                            ground_stations         : list, 
                                            node_index_file_path    : str,
                                            terrestrial_file_path   : str,
                                            timestamp               : int,
                                            operator_name           : str,
                                            t2t_dict                : dict = None
                                        ):
    """
    Saves the matching node indices and their corresponding aliases.

    Args:
        satellites_by_index (dict):     Satellites_by_index (dict): satellites sorted by index
                                        satellite/ground nodes
        ground_stations (list):         List of supplied ground stations
        node_index_file_path (str):     Path to output the node index matching file
        terrestrial_file_path (str):    Path to output the terrestrial node info file
        timestamp (int):                Unix time
        operator_name (str):            Constellation/operator name
        t2t_dict (dict):                Dictionary of terrestrial-to-terrestrial connectio

    Returns:
        Saves the node indices and terrestrial information as .txt files.
    """

    # Generate a new file for node indices
    file_path = node_index_file_path+operator_name+"/"
    check_create_path(file_path)
    file_name = "nodeindex_"+timestamp+".txt"
    nodeindex_log_write = open(file_path + file_name, "w") # Open file in write mode, overwriting any existing content
    nodeindex_log_write.close()

    # Generate a new file for terrestrial information
    file_path2 = terrestrial_file_path+"/"
    check_create_path(file_path2)
    file_name2 = "terrestrial_"+timestamp+".txt"
    terrestrial_log_write = open(file_path2 + file_name2, "w")
    terrestrial_log_write.close()

    # Append to files
    nodeindex_log = open(node_index_file_path+operator_name+"/nodeindex_"+timestamp+".txt", "a")
    terrestrial_log = open(terrestrial_file_path+"terrestrial_"+timestamp+".txt", "a")

    # Iterate over the satellites_by_index
    for sat_indx, sat_alias in satellites_by_index.items():
        nodeindex_log.write(str(sat_indx)+":"+str(sat_alias)+"\n")

    # Iterate over ground station list
    for gs in ground_stations:
        if gs['type'] == 9: # 9 indicates a gateway ground station
            alias_prefix = "GW-"
        else:
            alias_prefix = "CT-"  # all else are designated as customer terminals (ct)
        nodeindex_log.write(str(1+sat_indx+gs['gid'])+":"+alias_prefix+str(gs['gid'])+"\n")
        terrestrial_log.write(
                                str(1+sat_indx+gs['gid']) +
                                ":"+alias_prefix+str(gs['gid']) +
                                ":"+str(gs['name']) +
                                ":"+str(gs['latitude_degrees_str']) +
                                ":"+str(gs['longitude_degrees_str']) +
                                "\n"
                             )

    # If t2t_dict is included, iterate over the dictionary and append endpoint aliases to the file
    if t2t_dict is not None:
        alias_prefix = "IE-" # IE indicates an internet endpoint
        for id in t2t_dict.keys():
            if 'type' in t2t_dict[id] and t2t_dict[id]['type'] == 'endpoint':
                gs_coord = t2t_dict[id]['coordinates']
                gs_name = str(t2t_dict[id]['friendly_name']) if 'friendly_name' in t2t_dict.keys() else str(t2t_dict[id]['name'])
                nodeindex_log.write(str(1+sat_indx+t2t_dict[id]['gid'])+":"+alias_prefix+str(t2t_dict[id]['gid'])+"\n")
                terrestrial_log.write(
                                        str(1+sat_indx+t2t_dict[id]['gid']) +
                                        ":"+alias_prefix+str(t2t_dict[id]['gid']) +
                                        ":"+gs_name +
                                        ":"+str(gs_coord[0]) +
                                        ":"+str(gs_coord[1]) +
                                        "\n"
                                     )
    # Close file to minimize memory leaks
    nodeindex_log.close()
    terrestrial_log.close()

def save_cpu_time(
                    cpu_runtime     : float, 
                    timestamp       : int,
                    operator_name   : str, 
                    cpu_time_path   : str
                 ):
    """
    Saves the matching node indices and their corresponding aliases.

    Args:
        cpu_runtime (float):            CPU clock runtime, in seconds
        timestamp (list):               Unix time as a list
        operator_name (str):            Constellation/operator name
        cpu_time_path (str):            Path to output CPU clock runtime file  

    Returns:
        Saves the CPU clock runtime as a .txt file.
    """

    # Generate a new file
    cpu_clock_log = open(cpu_time_path+operator_name+"/cpu_clockruntime_"+("_".join(timestamp[:3]))+".txt", "w")
    
    # Iterate over the optimal path list
    if cpu_runtime != None:
        cpu_clock_log.write(str(cpu_runtime)+"\n")

    # Close file to minimize memory leaks
    cpu_clock_log.close()


# =================================================================================== #
# ------------------------------------ PARSING -------------------------------------- #
# =================================================================================== #
def parse_config_file_yml(
                            filepath        : str, 
                            configpath      : str,
                            main_config     : str,
                            satconfigpath   : str
                         ) -> dict:
    """
    Parse a YAML configuration file and load its contents into a dictionary.

    Args:
        filepath (str): Path to the directory containing the configuration file
        filename (str): Name of the configuration file

    Returns:
        dict:           Dictionary containing the parsed configuration parameters
    """


    # Obtain the constellation YAML file from main configuration YAML
    with open(filepath+"/"+filename, "r") as main_yml:
        main_cfg = yaml.safe_load(main_yml)
        const_filename = main_cfg.get('ConstellationName')

    # Open the YAML file in read mode
    with open(filepath+"/../config_files/sat_config_files/"+const_filename+".yaml", "r") as ymlfile:

        # Use yaml.safe_load to parse the YAML content into a dictionary
        cfg = yaml.safe_load(ymlfile)

    # Return the parsed configurations
    return cfg


def parse_config_file(
                        filepath    : str, 
                        filename    : str
                     ) -> dict:
    """
    Parse a configuration file and extract key-value pairs.

    Args:
        filepath (str): Path to the directory containing the configuration file
        filename (str): Name of the configuration file

    Returns:
        dict:           Dictionary containing parsed configuration parameters
    """

    # Default configuration parameters
    configurations = {
                        "topology_path": "", 
                        "topology_routes_path": "", 
                        "simulation_results_n_data": "", 
                        "simulation_start_time": "", 
                        "simulation_time(second)": 0, 
                        "mode": 0, 
                        "simulation_step(second)": 0, 
                        "Fresh_run": False, 
                        "ground_stations": "./", 
                        "inclination": 0, 
                        "constellation": "starlink", 
                        "tle_file": "", 
                        "number_of_orbits": 0, 
                        "number_of_sat_per_orbit": 0, 
                        "altitude": 0, 
                        "elevation_angle": 0, 
                        "experiment": 2, 
                        "constellation_ip_range": "", 
                        "other_constellation_ip_range": "", 
                        "association_criteria": "BASED_ON_DISTANCE_ONLY_MININET" ,
                        "exit_gw": "",
                        "interDomain_routing": 0, 
                        "False_run_archieve_path_foldername": "",
                        "Debug": 1
                     }
    
    # Open and read the configuration file
    configFile = open(filepath+"/"+filename, 'r')
    configs = configFile.readlines()

    # Parse each line in the configuration file
    for config in configs:

        # Split string to array
        config_parameters = config.split("=")

        # Check if the value is a digit (int)
        if config_parameters[1].strip().isdigit():
            configurations[str(config_parameters[0])]=int(config_parameters[1].strip())
        
        # Check if the value is either "False" or "True" (bool)
        elif config_parameters[1].strip() == "False" or config_parameters[1].strip() == "True":
            if config_parameters[1].strip() == "False":
                configurations[str(config_parameters[0])] = False
            elif config_parameters[1].strip() == "True":
                configurations[str(config_parameters[0])] = True
        else:
            configurations[str(config_parameters[0])]=config_parameters[1].strip()

    # Return parsed configurations
    return configurations


def parse_connectivity_matrix_n_charateristics(
                                                time_utc        : object, 
                                                conn_mat_size   : int, 
                                                topology_path   : str
                                              ) -> dict:
    """
    Parse and retrieve connectivity matrix, links latency, and links capacity information at a specific time

    Args:
        time_utc (object):          Specified time expressed in UTC using Skyfield time object
        conn_mat_size (int):        Size of the connectivity matrix
        topology_path (str):        Path to the directory containing topology files

    Returns:
        dict:                       Dictionary containing parsed connectivity matrix, links latency, and links capacity
    """

    # Convert the provided UTC timestamp to individual time components
    time_components         = convert_time_utc_to_ymdhms(time_utc)
    topology_file_found     = 0

    # Initialize matrices for connectivity, links latency, and links capacity
    connectivity_matrix     = [[0 for c in range(conn_mat_size)] for r in range(conn_mat_size)]
    links_latency           = [[0 for c in range(conn_mat_size)] for r in range(conn_mat_size)]
    links_capacity          = [[0 for c in range(conn_mat_size)] for r in range(conn_mat_size)]

    # Iterate through files in the specified directory
    for filename in os.listdir(topology_path):

        # Combine the directory paths
        f = os.path.join(topology_path, filename)

        # Check if the current item is a file
        if os.path.isfile(f):

            # Check if the components match the target time
            if int(time_components[0]) == int(filename.split("_")[1]) and int(time_components[1]) == int(filename.split("_")[2]) and int(time_components[2]) == int(filename.split("_")[3]) and int(time_components[3]) == int(filename.split("_")[4]) and int(time_components[4]) == int(filename.split("_")[5]) and int(time_components[5]) == int(float(filename.split("_")[6][:-4])):
                
                # Initialize variable
                topology_file_found = 1

                # Open and read the topology file
                topology_file = open(f, 'r')
                lines = topology_file.readlines()

                # Process each line in the topology file
                for line in lines:

                    # Split string to array
                    link_config                                                                     = line.split(",")

                    # Update the connectivity matrix
                    connectivity_matrix[int(float(link_config[0]))][int(float(link_config[1]))]     = 1
                    connectivity_matrix[int(float(link_config[1]))][int(float(link_config[0]))]     = 1

                    # Update links latency matrices
                    links_latency[int(float(link_config[0]))][int(float(link_config[1]))]           = round(float(link_config[2]),0)
                    links_latency[int(float(link_config[1]))][int(float(link_config[0]))]           = round(float(link_config[2]),0)

                    # Update links capacity matrices
                    links_capacity[int(float(link_config[0]))][int(float(link_config[1]))]          = round(float(link_config[3]),0)
                    links_capacity[int(float(link_config[1]))][int(float(link_config[0]))]          = round(float(link_config[3]),0)

                # Break after processing the first matching file 
                break

    # Check if a topology file was found
    if topology_file_found == 0:
        print(time_components)
        print("[Error] No Topology file available ... check the simulation step resolution")
        return -1

    # Return parsed connectivity data
    return {
            "connectivity_matrix": connectivity_matrix,
            "links_latency":       links_latency,
            "links_capacity":      links_capacity
            }


def parse_topology_routes(
                            topology_route_path : str, 
                            num_of_satellites   : int, 
                            time_utc            : object
                         ) -> dict:
    """
    Parse and retrieve topology routes information at a specified time.

    Args:
        topology_route_path (str):  Path to the directory containing topology route files
        num_of_satellites (int):    Number of satellite nodes in network
        time_utc (object):          Specified time expressed in UTC using Skyfield time object

    Returns:
        dict:                       Dictionary containing parsed topology route information.
    """

    # Convert the provided UTC timestamp to individual time
    time_components         = convert_time_utc_to_ymdhms(time_utc)
    routes_file_found       = 0

    # Iterate through files in the specified directory
    for filename in os.listdir(topology_route_path):

        # Combine the directory paths
        f = os.path.join(topology_route_path, filename)

        # Check if the current item is a file
        if os.path.isfile(f):

            # Check if the components match the target time
            if int(time_components[0]) == int(filename.split("_")[1]) and int(time_components[1]) == int(filename.split("_")[2]) and int(time_components[2]) == int(filename.split("_")[3]) and int(time_components[3]) == int(filename.split("_")[4]) and int(time_components[4]) == int(filename.split("_")[5]) and int(time_components[5]) == int(float(filename.split("_")[6][:-4])):
                
                # Initialize data structures
                constellation_routes    = {m: [] for m in range(num_of_satellites)}
                initial_routes          = []
                thread_list             = []
                routes_file_found       = 1
                Rfilename               = f
                route_file              = open(Rfilename, 'r')
                routes                  = route_file.readlines()
                num_thread              = 1000
                sublist_len             = len(routes)//num_thread

                # Divide the routes into sublists for parallel processing
                for i in range(0, len(routes), sublist_len):
                    subroutes = routes[i:i+sublist_len]
                    thread = threading.Thread(target=read_IProute_files_thread, args=(subroutes, initial_routes))
                    thread_list.append(thread)

                # Start and join threads for parallel route parsing
                for thread in thread_list:
                    thread.start()
                for thread in thread_list:
                    thread.join()

                # Copy and organize the initial routes into the constellation_routes dictionary
                copy_initial_routes     = initial_routes[:]
                for route in copy_initial_routes:
                    if len(route[0]) != 0:
                        current_route = route[0]
                        i = current_route[0]
                        constellation_routes[i].append(route)

                # Break after processing the first matching file
                break

    # Check if a routes file was found
    if routes_file_found == 0:
        print("[Error] No Route file available ... check the simulation step resolution")
        return -1

    # Return parsed routes data
    return {
                "All_PreConfigured_routes": initial_routes,
                "Routes_per_satellites": constellation_routes
            }


def parse_resiliency_experiment_parameters(
                                            main_configurations : dict
                                          ) -> dict:
    """
    Logs the start and stop timestamps of resiliency event of each affected satellite from the simulation configuration file as a dictionary.

    Args:
        main_configurations (dict): Simulation definitions from the YAML configuration file

    Returns:
        dict:                       Resiliency event timestamps for each affected satellite
    """

    # Initialize empty dictionary
    resiliency_satellite_timestamp = {}

    # Iterate over affected satellites
    for sats, s_timestamps, e_timestamps in zip(main_configurations["resiliency"]["affected_satellites"], main_configurations["resiliency"]["s_timestamps"], main_configurations["resiliency"]["e_timestamps"]):
        resiliency_satellite_timestamp[sats] = (s_timestamps,e_timestamps)

    # Return the dictionary of resiliency timestamps per satellite
    return resiliency_satellite_timestamp


# =================================================================================== #
# ----------------------------------- SATELLITE ------------------------------------- #
# =================================================================================== #
def arrange_satellites(
                        orbital_data            : dict, 
                        satellites_by_name      : dict, 
                        sat_config              : dict,
                        operator_name           : str,
                        satellites_by_index     : {},
                        timestamp               : object,
                        tle_timestamp           : str,
                        sat_orbit_file_path     : str
                      ) -> dict:
    """
    Arranges the satellites by index by transforming a TLE file into sorted orbits.

    Args:
        orbital_data (dict):                Extracted orbital parameters for each satellite in TLE
        satellites_by_name (dict):          Satellite information arranged by name
        sat_config (dict):                  Constellation configuration file
        operator_name (name):               Constellation/operator name
        satellites_by_index (empty dic):    Satellite information arranged by index, given as an empty dictionary
        timestamp (object):                 Skyfield object datetime
        tle_timestamp (str):                Attached unix timestamp of the TLE file
        sat_orbit_file_path (str):          Path to output the satellite orbit files

    Returns:
        dict:                               Appended satellites_by_index, with sorted satellites in their orbits
    """

    # Generate a new file
    file_path = sat_orbit_file_path+operator_name+"/"
    check_create_path(file_path)
    file_name = "sorted_satellites_within_orbit_"+tle_timestamp+".txt"
    #f = open(sat_orbit_file_path+operator_name+"/sorted_satellites_within_orbit_"+tle_timestamp+".txt", "a")
    f0 = open(file_path + file_name, "w")
    f0.close()
    f = open(file_path + file_name, "a")
    
    # Debugging purposes
    if sat_config["Debug"] == 1:
        print("..... Phase-1: Constellation Orbits:")

    # Initialize satellite names according to constellation naming conversion
    satellites_sorted_in_orbits = []

    # Iterate over the number of orbits in first shell
    for i in range(sat_config["shell1"]["orbits"]):
        sorted = []
        satellites_in_orbit = []
        cn = 0
        for data in orbital_data:
            if i == int(orbital_data[str(data)][0]):
                satellites_in_orbit.append(satellites_by_name[str(data.split(" ")[0])])
                cn +=1

        # Sort the satellites in orbit and append them to list
        sorted = sort_satellites_in_orbit(satellites_in_orbit, timestamp)
        satellites_sorted_in_orbits.append(sorted)

        # Debugging purposes
        if sat_config["Debug"] == 1:
            print(".......... Orbit no.    "+str(i)+"    ->  "+str(cn)+" satellites")

        # Debugging purposes
        if sat_config["Debug"] == 1:
            for s in sorted:
                write_this = str(i)+" "+str(s.name)+" "+str(orbital_data[str(s.name)])+"\n"
                f.write(write_this)
    
    # Close debugging file to minimize memory leak
    f.close()
   
    # Generate new file for sorted satellites in their orbits
    file = open(sat_orbit_file_path+"orbits_satellites.txt", 'w')
    sat_index = -1
    orbit_id = 0

    # Append new information to provided empty dictionary
    for orbit in satellites_sorted_in_orbits:
        orbit_id += 1
        for i in range(len(orbit)):
            sat_index += 1
            satellites_by_index[sat_index] = orbit[i].name.split(" ")[0]
            file.write(str(orbit_id)+"\t"+str(sat_index)+"\n")

    # Close file to minimize memory leak
    file.close()

    # Return dictionary
    return {"sorted satellite in orbits": satellites_sorted_in_orbits,
            "satellites by index": satellites_by_index
            }


def reload_tles(
                    path_of_recent_TLE      : str, 
                    main_configurations     : dict
               ) -> dict:
    """
    Retrieves the orbit and satellite information from a new TLE file.

    Args:
        path_of_recent_TLE (str):   Directory of new TLE file
        main_configurations (dict): Simulation definitions from the YAML configuration file

    Returns:
        dict:                       Large dictionary of satellite information, their sorted and arranged orbits, and number of satellites
                                    from the new TLE file
    """

    # Get the timestamp of new TLE
    tle_timestamp = path_of_recent_TLE.split("_")[2]

    # Load TLE data using Skyfield
    satellites = load.tle_file(path_of_recent_TLE)

    # Arrange satellites by their naming convention
    satellites_by_name = {sat.name.split(" ")[0]: sat for sat in satellites}

    # Initialize empty dictionary
    satellites_by_index = {}

    # Get the orbital plane classifications
    orbital_data  = get_orbital_planes_classifications(path_of_recent_TLE, main_configurations["constellation"]["operator"], main_configurations["constellation"]["shell1"]["orbits"], main_configurations["constellation"]["shell1"]["sat_per_orbit"], main_configurations["constellation"]["shell1"]["inclination"])
    
    # Arrange satellites by index and within their sorted orbits
    arranged_sats = arrange_satellites(orbital_data, satellites_by_name, main_configurations, time_utc ,satellites_by_index, tle_timestamp)
    
    # Split the information based on arranged indicies and sortment
    satellites_by_index = arranged_sats["satellites by index"]
    satellites_sorted_in_orbits = arranged_sats["sorted satellite in orbits"]

    # Number of sat and gs nodes
    num_of_satellites = len(orbital_data)
    num_of_ground_stations = len(ground_stations)

    # Debugging purposes
    if main_configurations["Debug"] == 1:
        print("................................. Re Loading the new TLE files ...........................")
        print(".......... total number of satellites = ", num_of_satellites)
        print(".......... total number of ground_stations = ", num_of_ground_stations)

    # Return new dictionary of satellite information
    return {"orbital_data": orbital_data,
            "satellites_by_name": satellites_by_name,
            "satellites_by_index": satellites_by_index,
            "satellites_sorted_in_orbits": satellites_sorted_in_orbits,
            "num_of_satellites": num_of_satellites
    }


def get_gs_sat_pairs(
                        connectivity_matrix : list, 
                        num_of_satellites   : int
                    ) -> list:
    """
    This function generates pairs of ground stations and satellites that are connected.

    Args:
        connectivity_matrix (list): A 2D matrix representing the connectivity between satellites and ground stations
                                    If connectivity_matrix[i][j] == 1, it means satellite i and ground station j are connected
        num_of_satellites (int):    Rotal number of satellites in constellation

    Returns:
        list:  List of tuples where each tuple represents a pair of connected satellite and ground station
               The first element of the tuple is the satellite index and the second element is the ground station index
    """

    # Initialize an empty list to store the pairs
    pairs = []

    # Iterate over the rows of the connectivity matrix
    for i in range(len(connectivity_matrix)):

        # Iterate over the columns of the connectivity matrix
        for j in range(len(connectivity_matrix[i])):

            # If the element at position (i, j) is 1 and i is less than the number of satellites and j is greater or equal to the number of satellites
            if connectivity_matrix[i][j] == 1 and i < num_of_satellites and j >= num_of_satellites:

                # Append the pair (i, j) to the list of pairs
                pairs.append((i, j))

    # Return the list of pairs
    return pairs


def get_sats_by_index(filename: str) -> list:
    """
    This function reads a file containing satellite names and returns a list of satellites.

    Args:
        filename (str):     The name of the file containing the satellite names
                            Each line in the file represents a satellite

    Returns:
        satellites (list):  A list of satellite names; the index of each satellite in the list corresponds to 
                            its index in the file
    """

    # Open the file in read mode
    satsFile = open(filename, 'r')

    # Read all lines from the file
    lines = satsFile.readlines()

    # Initialize an empty list to store the satellite names
    satellites = []

    # Iterate over the lines in the file
    for i in range(len(lines)):

        # Strip the newline character from the end of the line and append the resulting string to the list of satellites
        satellites.append(lines[i].strip())

    # Return the list of satellites
    return satellites


def get_sats_by_name(filename: str) -> list:
    """
    This function reads a file containing satellite names and returns a list of satellites.

    Args:
        filename (str):     The name of the file containing the satellite names
                            Each line in the file represents a satellite

    Returns:
        satellites (list):  A list of satellite names; the order of satellites in the list corresponds 
                            to their order in the file
    """

    # Open the file in read mode
    satsFile = open(filename, 'r')

    # Read all lines from the file
    lines = satsFile.readlines()

    # Initialize an empty list to store the satellite names
    satellites = []

    # Iterate over the lines in the file
    for i in range(len(lines)):

        # Strip the newline character from the end of the line and append the resulting string to the list of satellites
        satellites.append(lines[i].strip())

    # Return the list of satellites
    return satellites


# =================================================================================== #
# ----------------------------------- CONVERSION ------------------------------------ #
# =================================================================================== #

def convert_time_utc_to_unix(
                                time_utc    : object
                            ) -> float:
    """
    Converts a datetime expressed as Skyfield time object to unix as a float.

    Args:
        time_utc (object):  Specified time expressed in UTC using Skyfield time object

    Returns:
        float:              Unix time
    """

    # Convert time object to a string
    time_utc_string = time_utc.utc_strftime()

    # Convert string to 'datetime' object
    time_datetime = datetime.strptime(time_utc_string, "%Y-%m-%d %H:%M:%S %Z")

    # Convert the 'datetime' object to a float in local time
    time_timestamp = time.mktime(time_datetime.timetuple())

    # Convert the float back to UTC from local time
    datetime.now().isoformat()
    datetime.utcnow().isoformat()
    time.altzone
    time_timestamp = time_timestamp - time.altzone

    # Return conversion to unix time
    return time_timestamp


def convert_time_utc_to_ymdhms(
                                time_utc    : object
                              ) -> tuple:
    """
    Converts a datetime expressed as Skyfield time object to epoch time components as a tuple.

    Args:
        time_utc (object):  Specified time expressed in UTC using Skyfield time object

    Returns:
        tuple:              Time components of an epoch split into (Y, M, D, H, M, S)
    """

    # Convert time object to a string
    time_utc_string = time_utc.utc_strftime()

    # Split time information
    date, t_time, _ = time_utc_string.split(" ")

    # Get the Y/M/D information from date
    year, month, day = date.split("-")

    # Get the hour:minute:second from time
    hour, minute, second = t_time.split(":")

    # Return as a tuple
    return (year, month, day, hour, minute, second)


# =================================================================================== #
# ------------------------------------ NETWORK -------------------------------------- #
# =================================================================================== #
def check_changes_in_link_charateristics(
                                            old_links_characteristics   : list, 
                                            new_links_characteristics   : list
                                        ) -> list:
    """
    Compares whether there's any contrast between link characteristics at times (t-1) and (t).

    Args:
        old_links_characteristics (list):   Link characteristics at (t-1)
        new_links_characteristics (list):   Link characteristics at (t)

    Returns:
        list:                               Log of link characteristics at times (t-1) and (t) and their corresponding indices
    """

    # Initialize updated list
    changes = []

    # Iterate over and check for contrast between two link sets
    for i in range(len(new_links_characteristics)):
        for j in range(len(new_links_characteristics[i])):
            if new_links_characteristics[i][j] != old_links_characteristics[i][j]:
                changes.append((i, j, old_links_characteristics[i][j], new_links_characteristics[i][j]))

    # Return updated list
    return changes


def merge_link_link_charateristics(
                                    latency_changes     : list, 
                                    capacity_changes    : list
                                  ) -> tuple:
    """
    *TO DO: Not sure what this is doing 

    Args:
        latency_changes (list):     Log of latency values at times (t-1) and (t) and their corresponding indices   
        capacity_changes (list):    Log of capacity values at times (t-1) and (t) and their corresponding indices

    Returns:
        tuple:                      *TO DO: Not sure what this actually outputs               
    """

    # Initialize merged changes
    merged_changes = []
    merged_changes2 = []

    # Merge changes from latency with corresponding changes in capacity
    for lchange in latency_changes:
        found = 0
        for cchange in capacity_changes:

            # Check if the changes in latency and capacity happen at same 2 instances
            if (lchange[0] == cchange[0]) and (lchange[1] == cchange[1]):
                merged_changes.append((lchange[0], lchange[1], lchange[3], cchange[3]))
                found = 1
        if found == 0:
            merged_changes.append((lchange[0], lchange[1], lchange[3]))

    # Merge changes from capacity with corresponding merged changes
    for change in capacity_changes:
        found = 0
        for mchange in merged_changes:

            # Check if the changes in latency and capacity happen at same 2 instances
            if (change[0] == mchange[0]) and (change[1] == mchange[1]):
                found = 1
        if found == 0:
            merged_changes2.append((change[0], change[1], change[3]))

    # Return merged changes in both cases
    return (merged_changes, merged_changes2)


def apply_topology_updates_to_mininet(
                                        path                : str, 
                                        net                 : object, 
                                        topology_changes    : list,  
                                        num_of_satellites   : int, 
                                        time_utc_inc        : object
                                     ) -> object:
    """
    Apply changes in the topology to Mininet.

    Args:
        path (str):                 Path to the bash script files
        net (obj):                  Mininet network object
        topology_changes (list):    List of topology changes, each represented as a tuple (node1, node2, add_link, remove_link)
        num_of_satellites (int):    Number of satellite nodes in network
        time_utc_int (obj):         Object representing a timestamp increment

    Returns:
        net:                        Updated Mininet network object
    """

    # Iterate through each topology change
    for change in topology_changes:

        # Determine the node names based on the change
        node1 = "sat"+str(change[0]) if int(change[0]) < num_of_satellites else "gs"+str(int(change[0])%num_of_satellites)
        node2 = "sat"+str(change[1]) if int(change[1]) < num_of_satellites else "gs"+str(int(change[1])%num_of_satellites)

        # Get Mininet node objects based on node names
        net_node1 = net.getNodeByName(node1)
        net_node2 = net.getNodeByName(node2)

        # If the change indicates removal of a link
        if change[2] == 1 and change[3] == 0:

            # Check if a link exists before deleting it
            if net.linksBetween(net_node1, net_node2):
                net.delLinkBetween(net_node1, net_node2)
                print("[Info] the link between ", str(node1), "and", str(node2), "is deleted ...")
            
            else: # Print an error message if the link does not exist
                print("[Error] the link does not exist ... check apply_updates_to_mininet function", str(node1), "--", str(node2))

    # Iterate through each topology change again
    for change in topology_changes:

        # Determine the node names based on the change
        node1 = "sat"+str(change[0]) if int(change[0]) < num_of_satellites else "gs"+str(int(change[0])%num_of_satellites)
        node2 = "sat"+str(change[1]) if int(change[1]) < num_of_satellites else "gs"+str(int(change[1])%num_of_satellites)

        # Get Mininet node objects based on node names
        net_node1 = net.getNodeByName(node1)
        net_node2 = net.getNodeByName(node2)

        # If the change indicates addition of a link
        if change[2] == 0 and change[3] == 1:

            # Add a link between the nodes
            net.addLink(net_node1, net_node2, cls=TCLink)
            # print("[Info] the link between ", str(node1), "and", str(node2), "is added ...")

    # Iterate through satellite nodes
    for i in range(0, num_of_satellites):

         # Get Mininet node object for the satellite node
        sat_node = net.getNodeByName("sat"+str(i))

        # Generate a command to execute a shell script for route updates
        command = "../"+path+"/routes_updates_"+str(time_utc_inc.utc_strftime())+"/sat"+str(i)+"_routes.sh"
        
        # pkill_command = "pkill -f "+"*_routes.sh"
        # sat_node.cmd(pkill_command)

        # Execute the command on the satellite node
        sat_node.cmd(command)

     # Return the updated Mininet network object
    return net


def apply_link_updates_to_mininet(
                                    net                 : object, 
                                    latency_changes     : list,
                                    capacity_changes    : list,  
                                    num_of_satellites   : int, 
                                    time_utc_inc        : object
                                 ) -> object:
    """
    Apply changes in the link characteristics to Mininet.

    Args:
        net (obj):                  Mininet network object
        latency_changes (list):     List of latency changes, each represented as a tuple (node1, node2, latency, capacity)
        capacity_changes (list):    List of capacity changes, each represented as a tuple (node1, node2, capacity)
        num_of_satellites (int):    Number of satellite nodes in network
        time_utc_int (obj):         Object representing a timestamp increment

    Returns:
        net:                        Updated Mininet network object
    """

    # Iterate through each latency change
    for change in latency_changes:

        # Determine the node names based on the change
        node1 = "sat"+str(change[0]) if int(change[0]) < num_of_satellites else "gs"+str(int(change[0])%num_of_satellites)
        node2 = "sat"+str(change[1]) if int(change[1]) < num_of_satellites else "gs"+str(int(change[1])%num_of_satellites)

        # Get Mininet node objects based on node names
        net_node1 = net.getNodeByName(node1)
        net_node2 = net.getNodeByName(node2)

        # Check if a link exists before updating latency
        if net.linksBetween(net_node1, net_node2):

            # Get the link object
            link = net.linksBetween(net_node1, net_node2)

            # Configure the latency of both interfaces in the link
            link[0].intf1.config(delay=str(change[2])+'ms')
            link[0].intf2.config(delay=str(change[2])+'ms')
            # print("[Info] the latency of link between ", str(node1), "and", str(node2), "has changed to", str(change[2]), "ms")

            # Check if there's a latency change
            if len(change) == 4:
                link[0].intf1.config(bw=float(change[3]))
                link[0].intf2.config(bw=float(change[3]))
                # print("[Info] the capacity of link between ", str(node1), "and", str(node2), "has changed to", str(change[3]), "Mbps")

    # Iterate through each capacity change
    for change in capacity_changes:

        # Determine the node names based on the change
        node1 = "sat"+str(change[0]) if int(change[0]) < num_of_satellites else "gs"+str(int(change[0])%num_of_satellites)
        node2 = "sat"+str(change[1]) if int(change[1]) < num_of_satellites else "gs"+str(int(change[1])%num_of_satellites)

        # Get Mininet node objects based on node names
        net_node1 = net.getNodeByName(node1)
        net_node2 = net.getNodeByName(node2)

        # Check if a link exists before updating capacity
        if net.linksBetween(net_node1, net_node2):

            # Get the link object
            link = net.linksBetween(net_node1, net_node2)

            # Configure the capacity of both interfaces in the link
            link[0].intf1.config(bw=float(change[2]))
            link[0].intf2.config(bw=float(change[2]))

            #print("[Info] the capacity of link between ", str(node1), "and", str(node2), "has changed to", str(change[2]), "Mbps")

    # Return the updated Mininet network object
    return net


# =================================================================================== #
# ------------------------------- PERFORMANCE TESTING ------------------------------- #
# =================================================================================== #
def iperf_app(
                path                    : str, 
                net                     : object,
                main_configurations     : dict,
                list_of_Intf_IPs        : dict
             ) -> object:
    """
    Generates and configures IPerf commands for throughput testing using Mininet.

    Args:
        path (str):                 Path to the bash script files (unused)
        net (object):               Mininet network object
        main_configurations (dict): Simulation definitions from the YAML configuration file
        list_of_Intf_IPs (dict):    Interface IPs corresponding to each node

    Returns:
        net:                        Updated Mininet network object
    """

    # Extract source and destination nodes from main configurations
    node1 = main_configurations["application"]["source"]
    node2 = main_configurations["application"]["destination"]

    # Get Mininet node objects based on node names
    net_node1 = net.getNodeByName(node1)
    net_node2 = net.getNodeByName(node2)

    # Extract source IP address from the list of interface IPs
    source_ip = list_of_Intf_IPs[str(node2)+"-eth0"][0].split("/")[0]

    # Generate the IPerf command for the source node
    iperf_command = "iperf -c "+str(source_ip)+" -p 5201 -i1 -t"+str(main_configurations["application"]["duration"])+ " > "+str(main_configurations["application"]["result_out"])+ " &"
    
    # print iperf_command, node1, node2

    # Start IPerf server on the destination node
    net_node2.cmd("iperf -s -p 5201 > iperf_src_out &")

    # Execute the IPerf command on the source node
    net_node1.cmd(iperf_command)

    # Return the updated Mininet network object
    return net


def ping_app(
                path                    : str, 
                net                     : object,
                main_configurations     : dict,
                list_of_Intf_IPs        : dict
            ) -> object:
    """
    Generates and configures Ping commands for latency testing using Mininet.

    Args:
        path (str):                 Path to the bash script files (unused)
        net (object):               Mininet network object
        main_configurations (dict): Simulation definitions from the YAML configuration file
        list_of_Intf_IPs (dict):    Interface IPs corresponding to each node

    Returns:
        net:                        Updated Mininet network object
    """

    # Extract source and destination nodes from main configurations
    node1 = main_configurations["application"]["source"]
    node2 = main_configurations["application"]["destination"]

    # Get Mininet node objects based on node names
    net_node1 = net.getNodeByName(node1)
    # net_node2 = net.getNodeByName(node2)

    # Extract source IP address from the list of interface IPs
    source_ip = list_of_Intf_IPs[str(node2)+"-eth0"][0].split("/")[0]

    # Generate the Ping command for the source node
    ping_command = "ping "+str(source_ip)+" > "+str(main_configurations["application"]["result_out"])+" &"

    # Execute the IPerf command on the source node
    net_node1.cmd(ping_command)

    # Return the updated Mininet network object
    return net


def run_application(
                        path                    : str, 
                        net                     : object,
                        main_configurations     : dict,
                        list_of_Intf_IPs        : dict
                   ) -> object:
    """
    Selects the network performance test defined by simulation configuration file and returns the updated Mininet object.

    Args:
        path (str):                 Path to the bash script files (unused)
        net (object):               Mininet network object
        main_configurations (dict): Simulation definitions from the YAML configuration file
        list_of_Intf_IPs (dict):    Interface IPs corresponding to each node

    Returns:
        net:                        Updated Mininet network object
    """

    # IPerf:    Throughput testing
    if main_configurations["application"]["type"] == "iperf":
        net = iperf_app(path, net, main_configurations, list_of_Intf_IPs)

    # Ping:     Latency testing
    if main_configurations["application"]["type"] == "ping":
        net = ping_app(path, net, main_configurations, list_of_Intf_IPs)

    # Return the updated Mininet network object
    return net


# =================================================================================== #
# --------------------------------- RESILIENCY APP ---------------------------------- #
# =================================================================================== #
def check_time_to_deploy_RE(
                                resiliency_satellite_timestamp  : dict, 
                                year                            : int, 
                                month                           : int, 
                                day                             : int, 
                                hour                            : int, 
                                minute                          : int, 
                                seconds                         : int 
                           ) -> list:
    """
    Check satellites that are scheduled for resiliency disruption event in the specified time.

    Args:
        resiliency_satellite_timestamp (dict):  Resiliency event timestamps for each affected satellite
        year (int):                             Year time component
        month (int):                            Month time component
        day (int):                              Day time component
        hour (int):                             Hour time component
        minute (int):                           Minute time component
        seconds (int):                          Seconds time component

    Returns:
        list:                                   List of satellite IDs scheduled for resiliency disruption event at the specified time
    """

    # Initialize output as empty list
    ret_sats = []

    # Iterate through each satellite in the disruption schedule
    for sats in resiliency_satellite_timestamp:

        # Split the timestamp into individual components
        RET_start_y,RET_start_m,RET_start_d,RET_start_h,RET_start_min,RET_start_s = resiliency_satellite_timestamp[sats][0].split(",")[0], resiliency_satellite_timestamp[sats][0].split(",")[1], resiliency_satellite_timestamp[sats][0].split(",")[2], resiliency_satellite_timestamp[sats][0].split(",")[3], resiliency_satellite_timestamp[sats][0].split(",")[4], resiliency_satellite_timestamp[sats][0].split(",")[5]
        
        # Check if the components match the target time
        if (
            int(RET_start_y) == int(year) and 
            int(RET_start_m) == int(month) and 
            int(RET_start_d) == int(day) and 
            int(RET_start_h) == int(hour) and 
            int(RET_start_min) == int(minute) and 
            int(RET_start_s) == int(seconds)
           ):
            
            # If a match is found, append the satellite ID to the result list
            ret_sats.append(sats[3:])

    # Return the list of satellite IDs scheduled for disruption event
    return ret_sats
# =================================================================================== #
# ------------------------------- IMPORT PACKAGES ----------------------------------- #
# =================================================================================== #

from tqdm import tqdm
from utils import *
import numpy as np
import re
from mobility.read_live_tles import *
from mobility.mobility_utils import *
from mobility.read_gs import *
from routing.routing_utils import *
from routing.constellation_routing import *
from utils.utils import *
from library import spacenet_yaml_config

# =================================================================================== #
# ---------------------------------- INPUT VARS ------------------------------------- #
# =================================================================================== #

find_optimal_routes         = True

# =================================================================================== #
# ---------------------------------- PARSE VARS ------------------------------------- #
# =================================================================================== #

config_file_path            = "config_files/"
config_file_name            = "main_mn_config.yaml"
sat_config_sub_path         = "sat_config_files/"
tle_file_path               = "utils/"
data_filepath               = "/mnt/c/Users/BluBoy/Desktop/Professional/Git/Repositories/dynamic-topology-generator/"
output_filepath             = data_filepath+"output/"
connectivity_matrix_path    = output_filepath+"connectivity/"
routing_file_path           = output_filepath+"routing/"
arranged_sat_file_path      = output_filepath+"general/"
sat_orbit_file_path         = output_filepath+"satellites_orbits/"
optimal_file_path           = output_filepath+"optimal_routes/"

# =================================================================================== #
# -------------------------------- MAIN FUNCTION ------------------------------------ #
# =================================================================================== #

def main():

    # Parse the main configurations from the YAML file
    main_config, sat_config = spacenet_yaml_config.load_sim_and_constellation_config_file(config_file_path, config_file_name, sat_config_sub_path)
    operator_name = re.match(r'[a-zA-Z]+', main_config["ConstellationName"]).group(0)

    # Load the timescale and initialize variables
    ts = load.timescale()
    inc = 0
    indx = 1
    time_resolution_in_seconds = sat_config["EpochIntervalDuration"]
    simulation_length = sat_config["EpochIntervalDuration"] * sat_config["EpochIntervalCount"]

    # Split the start time from the configurations into individual components
    epoch_start = (
    sat_config["EpochStartYear"],
    sat_config["EpochStartMonth"],
    sat_config["EpochStartDay"],
    sat_config["EpochStartHour"],
    sat_config["EpochStartMinute"],
    sat_config["EpochStartSecond"]
    )

    # Convert the start time to UTC and Unix timestamp
    time_utc = ts.utc(*map(int, epoch_start))
    time_timestamp = convert_time_utc_to_unix(time_utc)
    print(epoch_start)

    # Get the path of the most recent TLE file based on the timestamp
    path_of_recent_TLE  = get_recent_TLEs_using_timestamp(tle_file_path, time_timestamp, operator_name)
    tle_timestamp       = path_of_recent_TLE.split("_")[2]
    print(operator_name+tle_timestamp)

    # Load the satellites from the TLE file
    satellites = load.tle_file(path_of_recent_TLE)

    # Create dictionaries of satellites by name and index
    satellites_by_name = {sat.name.split(" ")[0]: sat for sat in satellites}
    satellites_by_index = {}

    # Read the ground stations from the file specified in the configurations
    ground_stations = read_gs(sat_config["GroundStationFile"])

    # Get the orbital data and arrange the satellites in the orbits
    orbital_data  = get_orbital_planes_classifications(path_of_recent_TLE, operator_name, sat_config["shell1"]["orbits"], sat_config["shell1"]["sat_per_orbit"], sat_config["shell1"]["inclination"], sat_config["shell1"]["altitude"])
    arranged_sats = arrange_satellites(orbital_data, satellites_by_name, sat_config, operator_name, satellites_by_index, time_utc, tle_timestamp, arranged_sat_file_path, sat_orbit_file_path)
    satellites_by_index = arranged_sats["satellites by index"]
    satellites_sorted_in_orbits = arranged_sats["sorted satellite in orbits"]

    # Get the total number of satellites and ground stations
    num_of_satellites = len(orbital_data)
    num_of_ground_stations = len(ground_stations)

    # Print debug information if enabled in the configurations
    if sat_config["Debug"] == 1:
        print(".......... Total number of satellites = ", num_of_satellites)
        print(".......... Total number of ground_stations = ", num_of_ground_stations, "\n")

    time_hist = np.arange(0.0, simulation_length, time_resolution_in_seconds)

    # Loop over the time history, update the topology and save it in a file
    for inc in tqdm(time_hist, total=len(time_hist), desc=r'.......... Creating topology files'):
        
        # Update the time
        indx += 1

        # Get the source and destination nodes
        source_node         = num_of_satellites + int(''.join(filter(str.isdigit, sat_config["Source"])))
        destination_node    = num_of_satellites + int(''.join(filter(str.isdigit, sat_config["Destination"])))
        optimal_path_nodes  = [source_node, destination_node]

        # Convert the updated time to UTC and Unix timestamp
        time_utc_inc = ts.utc(*map(int, epoch_start[:-1]), epoch_start[-1]+inc)
        y, mon, d, h, min, s = convert_time_utc_to_ymdhms(time_utc_inc)

        # Update the size of the connectivity matrix
        conn_mat_size = num_of_satellites + num_of_ground_stations

        # Initialize the connectivity matrix
        connectivity_matrix = [[0 for _ in range(conn_mat_size)] for r in range(conn_mat_size)]

        # Add ISLs to the connectivity matrix
        connectivity_matrix = mininet_add_ISLs(connectivity_matrix, satellites_sorted_in_orbits, satellites_by_name, satellites_by_index, "SAME_ORBIT_AND_GRID_ACROSS_ORBITS", time_utc_inc)

        # Add GSLs to the connectivity matrix
        connectivity_matrix = mininet_add_GSLs_parallel(connectivity_matrix, satellites_by_name, satellites_by_index, ground_stations, 2, sat_config["AssociationCritGSL"], time_utc_inc, sat_config, operator_name)

        # Calculate the link characteristics for GSLs and ISLs
        links_characteristics = calculate_link_characteristics_for_gsls_isls(connectivity_matrix, satellites_by_index, satellites_by_name, ground_stations, time_utc_inc)

        # Save the topology
        save_topology(connectivity_matrix, links_characteristics, operator_name, str(y)+"_"+str(mon)+"_"+str(d)+"_"+str(h)+"_"+str(min)+"_"+str(float(s)), connectivity_matrix_path)

        # Pre-compute the routing tables
        if find_optimal_routes:
            all_possible_routes, optimal_route = initial_routing_fw(satellites_by_index, ground_stations, connectivity_matrix, links_characteristics["latency_matrix"], optimal_path_nodes)
        else:
            all_possible_routes = initial_routing_fw(satellites_by_index, ground_stations, connectivity_matrix, links_characteristics["latency_matrix"], None)

        # Save the routes
        save_routes(all_possible_routes, operator_name, str(y)+"_"+str(mon)+"_"+str(d)+"_"+str(h)+"_"+str(min)+"_"+str(float(s)), routing_file_path)

        # Save the optimal routes between provided src/dest
        save_optimal_path(optimal_route, [str(y), str(mon), str(d), str(h), str(min), str(float(s))], operator_name, optimal_file_path)



if __name__ == '__main__':
    main()
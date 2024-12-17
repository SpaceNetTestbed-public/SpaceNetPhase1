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
from utils.gateway_utils import *
from concurrent.futures import ProcessPoolExecutor as ProcessExecutor
from resource_monitor.top_logger import TOP_LOGGER
import multiprocessing
import signal
import atexit

# =================================================================================== #
# ---------------------------------- INPUT VARS ------------------------------------- #
# =================================================================================== #

find_optimal_routes         = True
use_multiprocessing         = True
global_arranged_sats        = None
global_satellites_by_name   = None
plot_ground_stations        = False
run_resource_logger         = False

# =================================================================================== #
# ---------------------------------- PARSE VARS ------------------------------------- #
# =================================================================================== #

config_file_path            = "config_files/"
config_file_name            = "main_mn_config.yaml"
sat_config_sub_path         = "sat_config_files/"

def topology_generation(inc, sat_config, 
                        ts, epoch_start, 
                        num_of_satellites, 
                        num_of_ground_stations, 
                        ground_stations,
                        optimal_path_nodes, 
                        operator_name, 
                        main_config, 
                        t2t_dict, 
                        connectivity_matrix_path, 
                        routing_file_path,  
                        optimal_file_path,
                        optimal_weight_path,
                        cpu_time_path, M=None, e=None, n_seq=0, dotd=None):
        
        # Start CPU timer
        t0_it = time.perf_counter_ns()

        # Update the time
        num_of_ground_stations = len(ground_stations)
        arranged_sats = global_arranged_sats
        satellites_by_index = arranged_sats["satellites by index"]
        satellites_sorted_in_orbits = arranged_sats["sorted satellite in orbits"]
        satellites_by_name = global_satellites_by_name
        # sat1 = satellites_by_name[satellites_by_index[22]]
        # sat2 = satellites_by_name[satellites_by_index[23]]

        # Convert the updated time to UTC and Unix timestamp
        time_utc_inc = ts.utc(*map(int, epoch_start[:-1]), epoch_start[-1]+inc)
        y, mon, d, h, min, s = convert_time_utc_to_ymdhms(time_utc_inc)

        # Update the size of the connectivity matrix
        conn_mat_size = num_of_satellites + num_of_ground_stations

        # Initialize the connectivity matrix
        connectivity_matrix = [[0 for _ in range(conn_mat_size)] for r in range(conn_mat_size)]
        
        # Add ISLs to the connectivity matrix
        connectivity_matrix = mininet_add_ISLs(connectivity_matrix, satellites_sorted_in_orbits, satellites_by_name, satellites_by_index, sat_config["TopologyISL"], time_utc_inc, M, e, n_seq, dotd=dotd)

        # Add GSLs to the connectivity matrix
        connectivity_matrix = mininet_add_GSLs_parallel(connectivity_matrix, satellites_by_name, satellites_by_index, ground_stations, 2, sat_config["AssociationCritGSL"], time_utc_inc, sat_config, operator_name)

        # Calculate the link characteristics for GSLs and ISLs
        links_characteristics = calculate_link_characteristics_for_gsls_isls(connectivity_matrix, satellites_by_index, satellites_by_name, ground_stations, time_utc_inc)

        # Add t2t links to the connectivity matrix, if enabled
        if "TopoCrit" in main_config and int(main_config["TopoCrit"]) > 0:
            connectivity_matrix, links_characteristics, t2t_dict = add_t2t_links_to_connectivity_matrix(connectivity_matrix, links_characteristics, satellites_by_index, ground_stations, t2t_dict)
            #if inc == time_hist_initial: # Have first timestep update the node index file with Internet Endnodes (they have not yet been added)
                #update_node_index(t2t_dict, node_index_file_path, tle_timestamp, operator_name) # Update the node index file with Internet Endnodes (they have not yet been added)
        
        """
        topfile_path = connectivity_matrix_path+operator_name+"/topology_"+str(y)+"_"+str(mon)+"_"+str(d)+"_"+str(h)+"_"+str(min)+"_"+str(float(s))+".txt"
        if t2t_dict == None:
            connectivity_matrix, links_characteristics = extract_connectivity(topfile_path, conn_mat_size)
        else:
            connectivity_matrix, links_characteristics = extract_connectivity(topfile_path, conn_mat_size + len(t2t_dict))
        #print(topfile_path)
        """
        # Assign the metrics for routing
        metric_type = None # default (hops)
        if "RouteWeight" in main_config and str(main_config["RouteWeight"]):
            metric_type = str(main_config["RouteWeight"])
            if metric_type == "latency":
                metrics = links_characteristics["latency_matrix"]
            elif metric_type == "capacity":
                metrics = links_characteristics["throughput_matrix"]
            elif metric_type == "distance":
                metrics = links_characteristics["distance_matrix"]
            elif metric_type == "hops":
                metrics = None
        
        # Pre-compute the routing tables
        if find_optimal_routes:
            all_possible_routes, optimal_route, net_optimal_weight = initial_routing_fw(satellites_by_index, connectivity_matrix, metrics, optimal_path_nodes, criterion)
        else:
            all_possible_routes = initial_routing_fw(satellites_by_index, connectivity_matrix, metrics, None, criterion)
        print(net_optimal_weight)
        # Stop CPU timer
        dt_it = (time.perf_counter_ns() - t0_it) * 1e-9 # Convert to seconds
        
        # Save the topology
        if os.path.exists(connectivity_matrix_path+operator_name+"/topology_"+str(y)+"_"+str(mon)+"_"+str(d)+"_"+str(h)+"_"+str(min)+"_"+str(float(s))+".txt"): # Check if file already exists, if so then rewrite
            os.remove(connectivity_matrix_path+operator_name+"/topology_"+str(y)+"_"+str(mon)+"_"+str(d)+"_"+str(h)+"_"+str(min)+"_"+str(float(s))+".txt")
        save_topology(connectivity_matrix, links_characteristics, operator_name, str(y)+"_"+str(mon)+"_"+str(d)+"_"+str(h)+"_"+str(min)+"_"+str(float(s)), connectivity_matrix_path)
        
        # Save the routes
        if os.path.exists(routing_file_path+operator_name+"/routes_"+str(y)+"_"+str(mon)+"_"+str(d)+"_"+str(h)+"_"+str(min)+"_"+str(float(s))+".txt"): # Check if file already exists, if so then rewrite
            os.remove(routing_file_path+operator_name+"/routes_"+str(y)+"_"+str(mon)+"_"+str(d)+"_"+str(h)+"_"+str(min)+"_"+str(float(s))+".txt")
        save_routes(all_possible_routes, operator_name, str(y)+"_"+str(mon)+"_"+str(d)+"_"+str(h)+"_"+str(min)+"_"+str(float(s)), routing_file_path)
        
        # Save the optimal routes between provided src/dest
        save_optimal_path(optimal_route, [str(y), str(mon), str(d), str(h), str(min), str(float(s))], operator_name, optimal_file_path)

        save_optimal_weights(net_optimal_weight, [str(y), str(mon), str(d), str(h), str(min), str(float(s))], operator_name, optimal_weight_path)
        
        # Save CPU clock runtime
        #save_cpu_time(dt_it, [str(y), str(mon), str(d), str(h), str(min), str(float(s))], operator_name, cpu_time_path)

# =================================================================================== #
# -------------------------------- MAIN FUNCTION ------------------------------------ #
# =================================================================================== #

def main():

    # Set global variable
    global criterion
    criterion = 2 # default

    # Parse the main configurations from the YAML file
    main_config, sat_config = spacenet_yaml_config.load_sim_and_constellation_config_file(config_file_path, config_file_name, sat_config_sub_path)
    operator_name = re.match(r'[a-zA-Z]+', main_config["ConstellationName"]).group(0)

    # Path configuration
    output_filepath             = sat_config["OutputFilePath"]
    use_weather_data            = bool(main_config["UseWeatherData"]) if "UseWeatherData" in main_config else True
    run_resource_logger         = bool(main_config["MonitorResource"]) if "MonitorResource" in main_config else False
    if output_filepath[-1] == "/":
        output_filepath = output_filepath[:-1] # Remove the last slash if it exists
    gs_file_path                = sat_config["GroundStationFile"]
    tle_file_path               = sat_config["TLEFilePath"]
    connectivity_matrix_path    = output_filepath+"/connectivity_matrix/"
    routing_file_path           = output_filepath+"/routing/"
    sat_orbit_file_path         = output_filepath+"/satellites_orbits/"
    node_index_file_path        = output_filepath+"/node_indices/"
    terrestrial_file_path       = output_filepath+"/terrestrial_info/"
    optimal_file_path           = output_filepath+"/optimal_routes/"
    optimal_weight_path           = output_filepath+"/optimal_weights/"
    cpu_time_path               = output_filepath+"/cpu_time/"
    resource_path               = output_filepath+"/resource/"

    # Start the subprocess for resource logging
    if run_resource_logger:
        print("\n.......... Initiating resource logger")
        resource_log_process = multiprocessing.Process(target=TOP_LOGGER, args=(1, resource_path, '1584_10_10'))
        resource_log_process.start()
        atexit.register(lambda: os.kill(resource_log_process.pid, signal.SIGTERM))
        time.sleep(10)

    # Get the source and destination nodes
    source_node         = int(main_config["SourceDeviceName"]) #num_of_satellites + int(''.join(filter(str.isdigit, sat_config["Source"])))
    destination_node    = int(main_config["DestDeviceName"]) #num_of_satellites + int(''.join(filter(str.isdigit, sat_config["Destination"])))
    optimal_path_nodes  = [source_node, destination_node]

    # Load the timescale and initialize variables
    ts = load.timescale()
    inc = 0
    time_resolution_in_seconds = sat_config["EpochIntervalDuration"]
    simulation_length = time_resolution_in_seconds * sat_config["EpochIntervalCount"]

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

    # Get the path of the most recent TLE file based on the timestamp
    path_of_recent_TLE  = get_recent_TLEs_using_timestamp(tle_file_path, time_timestamp, operator_name)
    tle_timestamp       = path_of_recent_TLE.split("_")[2]
    print(tle_timestamp)
    print("\n\n..... Phase-0: Configuration Set-up:")
    print(".......... Operator Name: \t\t", operator_name)
    print(".......... Start Epoch: \t\t", datetime.fromtimestamp(int(tle_timestamp)).strftime('%B %d, %Y %H:%M:%S UTC'))
    print(".......... End Epoch: \t\t\t", datetime.fromtimestamp(int(tle_timestamp)+int(simulation_length)).strftime('%B %d, %Y %H:%M:%S UTC'))
    print(".......... Simulation Step-Size: \t", sat_config["EpochIntervalDuration"], "s")
    print(".......... Simulation Interval Count: \t", sat_config["EpochIntervalCount"])
    print(".......... Simulation Length: \t\t", simulation_length, "s") 
    print(".......... TLE File: \t\t\t", path_of_recent_TLE, "\n")

    # Load the satellites from the TLE file
    satellites = load.tle_file(path_of_recent_TLE)

    # Create dictionaries of satellites by name and index
    satellites_by_name = {sat.name.split(" ")[0]: sat for sat in satellites}
    satellites_by_index = {}

    # Read the ground stations from the file specified in the configurations
    ground_stations = read_gs(gs_file_path)
    num_assigned_gs = len(ground_stations)

    # If using t2t links, generage t2t dictionary, then add Gateways to ground stations
    t2t_dict = None
    if "TopoCrit" in main_config and int(main_config["TopoCrit"]) > 0:
        print(".......... Using T2T links. Collecting settings")
        t2t_settings = get_t2t_settings(main_config, output_filepath)
        print(".......... T2T settings collected. Loading T2T dictionary")
        t2t_dict = load_t2t_dict(t2t_settings)
        num_gateways = 0
        num_endpoints = 0
        for key in t2t_dict:
            if 'type' in t2t_dict[key]:
                if t2t_dict[key]['type'] == 'gateway':
                    num_gateways += 1
                if t2t_dict[key]['type'] == 'endpoint':
                    num_endpoints += 1
        # Add gateways to ground station list
        print(f".......... T2T dictionary loaded: Adding {num_gateways} Gateways to ground stations; {num_endpoints} Endpoints loaded.\n")
        ground_stations, t2t_dict = add_gateway_gs(ground_stations, t2t_dict) # Add gateways to ground stations (t2t_dict is updated with gid values for gateways and endpoints)
        criterion = int(main_config["TopoCrit"])
        
    # Get the orbital data and arrange the satellites in the orbits
    orbital_data  = get_orbital_planes_classifications(path_of_recent_TLE, operator_name, sat_config["shell1"]["orbits"], sat_config["shell1"]["sat_per_orbit"], sat_config["shell1"]["inclination"], sat_config["shell1"]["altitude"])
    arranged_sats = arrange_satellites(orbital_data, satellites_by_name, sat_config, operator_name, satellites_by_index, time_utc, tle_timestamp, sat_orbit_file_path)
    satellites_by_index = arranged_sats["satellites by index"]
    satellites_sorted_in_orbits = arranged_sats["sorted satellite in orbits"]

    # Get the total number of satellites and ground stations
    num_of_satellites = len(orbital_data)
    num_of_ground_stations = len(ground_stations)
    num_of_terrestrials = num_of_ground_stations
    if t2t_dict:
        num_of_terrestrials = len(t2t_dict) + num_assigned_gs

    # Print debug information if enabled in the configurations
    if sat_config["Debug"] == 1:
        print(".......... Total number of satellites = ", num_of_satellites)
        print(".......... Total number of ground stations = ", num_of_ground_stations)
        print(".......... Total number of terrestrial nodes = ", num_of_terrestrials)
        print(".......... Total number of network nodes = ", num_of_satellites+num_of_terrestrials)
        print(".......... Phase-1 complete.\n")

    # Instantiate simulation time history
    time_hist = np.arange(0.0, simulation_length, time_resolution_in_seconds)

    # Init motif, dotd
    if sat_config["TopologyISL"] == "MOTIF":
        M, e = motif_find_m_se_e(satellites_sorted_in_orbits, satellites_by_name, satellites_by_index, time_utc)
    else:
        M = None
        e = None
    n_seq = 0
    dotd = DoTD_History(num_of_satellites, len(time_hist))

    if plot_ground_stations:
        plot_nodes_links_ground_stations(None, ground_stations, None, operator_name, None, None, None, t2t_dict)

    # Start topology generation
    print("..... Phase-2: Building topology and connectivity matrices:")

    # Check if there's any files that exist
    y, mon, d, h, min, s = convert_time_utc_to_ymdhms(ts.utc(*map(int, epoch_start)))
    if  os.path.exists(connectivity_matrix_path+operator_name+"/topology_"+str(y)+"_"+str(mon)+"_"+str(d)+"_"+str(h)+"_"+str(min)+"_"+str(float(s))+".txt") \
        or os.path.exists(routing_file_path+operator_name+"/routes_"+str(y)+"_"+str(mon)+"_"+str(d)+"_"+str(h)+"_"+str(min)+"_"+str(float(s))+".txt") \
        or os.path.exists(optimal_file_path+operator_name+"/best_path_"+("_".join([str(y), str(mon), str(d)]))+".txt"):
            user_response = input(f"\033[91m.......... Files for this simulation already exists. Do you want to overwrite them?\033[0m (y/[n]): ") or 'n'
            if user_response.lower() == 'n':
                return
            else: print("\033[94m", end="")

    # Start CPU clock timer
    cpu_clock_tot_t0 = time.perf_counter_ns()

    # Save satellite and ground station indices
    save_node_index_and_terrestrial_info(satellites_by_index, ground_stations, node_index_file_path, terrestrial_file_path, tle_timestamp, operator_name, t2t_dict)

    # Remove any existing files
    if os.path.exists(optimal_file_path+operator_name+"/best_path_"+("_".join([str(y), str(mon), str(d)]))+".txt"):
        os.remove(optimal_file_path+operator_name+"/best_path_"+("_".join([str(y), str(mon), str(d)]))+".txt")
    if os.path.exists(cpu_time_path+operator_name+"/cpu_clockruntime_"+("_".join([str(y), str(mon), str(d)]))+".txt"):
        os.remove(cpu_time_path+operator_name+"/cpu_clockruntime_"+("_".join([str(y), str(mon), str(d)]))+".txt")
    
    # Start main simulation process
    global global_arranged_sats, global_satellites_by_name # Make these global for multiprocessing, but being used regardless
    global_arranged_sats = arranged_sats
    global_satellites_by_name = satellites_by_name
    if use_multiprocessing:
        print(f".......... Using multi-process execution for topology generation.\n.......... Operation started {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        # Before starting concurrent execution, get weather conditions for all ground stations to avoid excessive/unnecesary API calls
        print(".......... Preemptively getting weather data for all ground stations")
        from link.link_utils import get_weather_info
        recv_cnt = 0
        recving_weather_data = True
        for ground_station in tqdm(ground_stations, desc=".......... Getting weather data"):
            if recving_weather_data and use_weather_data:
                gs_lat = float(ground_station["latitude_degrees_str"])
                gs_lon = float(ground_station["longitude_degrees_str"])
                weather_data = get_weather_info(gs_lat, gs_lon)
                if weather_data != "":
                    ground_station["weather_data"] = weather_data
                    recv_cnt += 1
                    # Wait 1 second to avoid API rate limit
                    if recv_cnt % 25 == 0:
                        time.sleep(1)
                else:
                    recving_weather_data = False # Stop trying to get weather data
                    ground_station["weather_data"] = ""
            else:
                ground_station["weather_data"] = ""
        if not use_weather_data:
            print(".......... User decided not to use weather data for simulation")
        else:
            print(f".......... Weather data received for {recv_cnt} ground stations")
        
        # Execute simulation process
        with ProcessExecutor(max_workers=20) as executor:
            results = list(tqdm(executor.map(topology_generation,
                                             time_hist,
                                             [sat_config]*len(time_hist),
                                             [ts]*len(time_hist),
                                             [epoch_start]*len(time_hist),
                                             [num_of_satellites]*len(time_hist),
                                             [num_of_ground_stations]*len(time_hist),
                                             [ground_stations]*len(time_hist),
                                             [optimal_path_nodes]*len(time_hist),
                                             [operator_name]*len(time_hist),
                                             [main_config]*len(time_hist),
                                             [t2t_dict]*len(time_hist),                                                                                    
                                             [connectivity_matrix_path]*len(time_hist),
                                             [routing_file_path]*len(time_hist),
                                             [optimal_file_path]*len(time_hist),
                                             [optimal_weight_path]*len(time_hist),
                                             [cpu_time_path]*len(time_hist)),
                                total=len(time_hist), desc=r'.......... Computing network'))
    else:
        
        # Loop over the time history, update the topology and save it in a file
        for inc in tqdm(time_hist, total=len(time_hist), desc=r'.......... Computing network'):
            
            topology_generation(inc, sat_config, ts, epoch_start, num_of_satellites, num_of_ground_stations, ground_stations, optimal_path_nodes, operator_name, main_config, t2t_dict, connectivity_matrix_path, routing_file_path, optimal_file_path, optimal_weight_path, cpu_time_path, M, e, n_seq, dotd)
            
            
            """
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

            # Add t2t links to the connectivity matrix, if enabled
            if "Use_t2t" in main_config and bool(main_config["Use_t2t"]) == True:
                connectivity_matrix, links_characteristics, t2t_dict = add_t2t_links_to_connectivity_matrix(connectivity_matrix, links_characteristics, satellites_by_index, ground_stations, t2t_dict)

            # Save the topology
            if os.path.exists(connectivity_matrix_path+operator_name+"/topology_"+str(y)+"_"+str(mon)+"_"+str(d)+"_"+str(h)+"_"+str(min)+"_"+str(float(s))+".txt"): # Check if file already exists, if so then rewrite
                os.remove(connectivity_matrix_path+operator_name+"/topology_"+str(y)+"_"+str(mon)+"_"+str(d)+"_"+str(h)+"_"+str(min)+"_"+str(float(s))+".txt")
            save_topology(connectivity_matrix, links_characteristics, operator_name, str(y)+"_"+str(mon)+"_"+str(d)+"_"+str(h)+"_"+str(min)+"_"+str(float(s)), connectivity_matrix_path)

            # Pre-compute the routing tables
            # Add flag to include ground_stations in route calculations if using t2t links 
            if "Use_t2t" in main_config and bool(main_config["Use_t2t"]) == True:
                global route_to_gs
                route_to_gs = True
                
            if find_optimal_routes:
                all_possible_routes, optimal_route = initial_routing_fw(satellites_by_index, ground_stations, connectivity_matrix, links_characteristics["latency_matrix"], links_characteristics["distance_matrix"], optimal_path_nodes, route_to_gs)
            else:
                all_possible_routes = initial_routing_fw(satellites_by_index, ground_stations, connectivity_matrix, links_characteristics["distance_matrix"], None, route_to_gs)

            # Save the routes
            if os.path.exists(routing_file_path+operator_name+"/routes_"+str(y)+"_"+str(mon)+"_"+str(d)+"_"+str(h)+"_"+str(min)+"_"+str(float(s))+".txt"): # Check if file already exists, if so then rewrite
                os.remove(routing_file_path+operator_name+"/routes_"+str(y)+"_"+str(mon)+"_"+str(d)+"_"+str(h)+"_"+str(min)+"_"+str(float(s))+".txt")
            save_routes(all_possible_routes, operator_name, str(y)+"_"+str(mon)+"_"+str(d)+"_"+str(h)+"_"+str(min)+"_"+str(float(s)), routing_file_path)

            # Save the optimal routes between provided src/dest
            if inc == time_hist[0] and os.path.exists(optimal_file_path+operator_name+"/best_path_"+("_".join([str(y), str(mon), str(d)]))+".txt"): # Check if file already exists, if so then rewrite
                os.remove(optimal_file_path+operator_name+"/best_path_"+("_".join([str(y), str(mon), str(d)]))+".txt")
            save_optimal_path(optimal_route, [str(y), str(mon), str(d), str(h), str(min), str(float(s))], operator_name, optimal_file_path)
            """
    
    # Stop CPU clock timer and save total time
    executor.shutdown()
    cpu_clock_tot_dt = (time.perf_counter_ns() - cpu_clock_tot_t0) * 1e-9
    save_cpu_time("TOTSIM:"+str(cpu_clock_tot_dt), [str(y), str(mon), str(d), str(h), str(min), str(float(s))], operator_name, cpu_time_path)

    # Update progress
    print("\033[0m.......... Phase-2 complete. See the results under: "+output_filepath+"\n\n")

    # Kill resource logger process
    if run_resource_logger:
        os.kill(resource_log_process.pid, signal.SIGTERM)


if __name__ == '__main__':
    main()
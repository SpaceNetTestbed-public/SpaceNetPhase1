from skyfield.api import wgs84, load
import math
import threading

import sys
sys.path.append("../")
from link.link_utils import *


def calc_max_gsl_length(
                        sat_config,
                        operator_name
                        ):
    """
    Calculates the maximum Ground Station-to-Satellite Link (GSL) length

    Args:
        main_configurations (dict): simulation definitions from the YAML configuration file

    Returns:
        max_gsl_length_m (float): maximum gs-sat link length (in meters)
    """
    
    # Initialize return variable
    max_gsl_length_m = -1 
    
    # Check for starlink operator
    # if operator_name == "starlink":
        # Set a specific value for max GSL length
        # max_gsl_length_m = 2089686.4181956202 # same number used in Hypatia code, further reasoning behind this exact value is unknown
        # (additionally, the above value does not match the value one would get using the algorithm in the else case, but applied to a starlink case)
        
        # return max_gsl_length_m
    
    # Max GSL length for non-starlink operators
    # else:

    # Calculate satellite cone radius based on altitude and elevation angle
    satellite_cone_radius = (sat_config["shell1"]["altitude"])/math.tan(math.radians(sat_config["shell1"]["elevation_angle"]))
     
    # Calculate max GSL length using cone radius and satellite altitude, convert to meters
    max_gsl_length_m =  (math.sqrt(math.pow(satellite_cone_radius, 2) + math.pow(sat_config["shell1"]["altitude"], 2)))*1000
    
    return max_gsl_length_m

# removed calc_distance_gs_sat_worker, as it was only used in mininet_add_GSLs (which has also been removed)

def get_orbit_from_sat(
                       sat_id,
                       satellites_by_index,
                       satellites_by_name,
                       satellites_sorted_in_orbits):
    """
    Provides Orbit ID from satellite ID (Function for simplifying things!)
    """

    sat_name = satellites_by_index[sat_id]

    for i, orbit in enumerate(satellites_sorted_in_orbits):
        for sat in orbit:
            this_sat_name = list(satellites_by_name.keys())[list(satellites_by_name.values()).index(sat)]
            if this_sat_name == sat_name:
                return i


def calc_distance_gs_sat_thread(
                                ground_stations, 
                                satellites_by_name, 
                                satellites_by_index, 
                                time_t, 
                                max_gsl_length_m, 
                                ground_station_satellites_in_range
                                ):
    """
    Determines which ground stations are in range of each satellite.

    Args:
        ground_stations (dict): list of ground stations
        satellites_by_name (dict): satellites sorted by name
        satellites_by_index (dict): satellites sorted by index
        time_t (datetime): timestamp corresponding to current satellite locations
        max_gsl_length_m (float): maximum gs-sat link length (in meters)
        ground_station_satellites_in_range (dict): list containing gs identifiers, sat indices, and distances in between

    Returns:
        ground_station_satellites_in_range (dict): list containing gs identifiers, sat indices, and distances in between, including newly appended data

    """

    # Iterate over each ground station
    for gs in ground_stations:
        # Iterate over the range of satellite indices
        for sid in range(len(satellites_by_index)):
            # Calculate the distance between the current ground station and satellite
            distance_m = distance_between_ground_station_satellite(gs, satellites_by_name[str(satellites_by_index[sid])], time_t)
            
            # Check if the calculated distance is within the maximum GSL length
            if distance_m <= max_gsl_length_m:
                # If in range, append a tuple to the result list
                ground_station_satellites_in_range.append((distance_m, sid, gs["gid"]))

    # Return the list of valid ground station-satellite pairs
    return ground_station_satellites_in_range

# removed calc_distance_gs_sat_worker_alan, as it was only used in mininet_add_GSLs (which has also been removed)

def distance_between_ground_station_satellite(
                                              ground_station, 
                                              satellite, 
                                              t
                                              ):
    """
    Calculates the distance between a ground station and a satellite at a specific time

    Args:
        ground_station (object): ground station
        satellite (object): satellite
        t (datetime): time corresponding to current satellite position

    Returns:
        distance (float): distance between a ground station and a satellite (in meters)
    """
    
    # Convert ground station coordinates to a WGS84 latlon object
    bluffton = wgs84.latlon(float(ground_station["latitude_degrees_str"]), float(ground_station["longitude_degrees_str"]), ground_station["elevation_m_float"])
    
    # Calculate the difference vector between the satellite and ground station
    difference = satellite - bluffton

    # Transform the difference vector to topocentric coordinates
    topocentric = difference.at(t)

    # Get the altitude, azimuth, and distance from the topocentric coordinates
    alt, az, distance = topocentric.altaz()

    # Return the distance between the ground station and satellite in meters
    return distance.m

# removed distance_between_ground_station_satellite_alan, as it was only used in calc_distance_gs_sat_worker_alan

def calc_distance_sat_sat_thread(
                                current_sat_list, 
                                satellites_by_name, 
                                satellites_by_index, 
                                satellite_sorted_in_orbits,
                                time_t, 
                                max_isl_search_length, 
                                current_sat_satellites_in_range
                                ):
    """
    Determines which nearby satellites are in range of each satellite.

    Args:
        current_sat_list (dict): list of current orbit satellite (SKYFIELD type)
        satellites_by_name (dict): satellites sorted by name
        satellites_by_index (dict): satellites sorted by index
        time_t (datetime): timestamp corresponding to current satellite locations
        max_isl_search_length (float): maximum sat-sat link length (in meters)
        current_sat_satellites_in_range (dict): list containing current_sat identifiers, sat indices, and distances in between

    Returns:
        current_sat_satellites_in_range (dict): list containing newly appended current_sat identifiers, sat indices, and distances in between

    """
    current_sat_by_name = [sat.name.split(" ")[0] for sat in current_sat_list]   # STORES SATELLITE STRINGS IN THE ORDER
    current_sat_ids = [list(satellites_by_index.keys())[list(satellites_by_index.values()).index(s)] for s in current_sat_by_name]  # STORES SATELLITE INDICES IN THE ORDER
    #print("Interior thread --> " + str(current_sat_ids))

    for itr, current_sat in enumerate(current_sat_list):
        curr_id = list(satellites_by_index.keys())[list(satellites_by_index.values()).index(current_sat_by_name[itr])]
        # Iterate over the range of satellite indices
        for num, sid in enumerate(satellites_by_index.keys()):
            
            if sid is not curr_id:   # dont store itself
                # Calculate the distance between the current ground station and satellite
                distance_m = distance_between_two_satellites(current_sat, satellites_by_name[str(satellites_by_index[sid])], time_t)
                # Check if the calculated distance is within the maximum GSL length
                if distance_m <= max_isl_search_length:
                    # If in range, append a tuple to the result list
                    current_sat_satellites_in_range.append((distance_m, sid, current_sat_by_name[itr]))

     # Return the list of valid ground station-satellite pairs
    return current_sat_satellites_in_range

def distance_between_two_satellites(
                                    satellite1, 
                                    satellite2, 
                                    t
                                    ):
    """
    Calculates the distance between two satellites

    Args:
        satellite1 (object): first relevant satellite
        satellite2 (object): second relevant satellite
        t (datetime): time corresponding to current satellite positions

    Returns:
        distance (float): distance between the two relevant satellites (in meters)
    """
    
    # Get the position of the first satellite at the given time
    position1 = satellite1.at(t)

    # Get the position of the second satellite at the given time
    position2 = satellite2.at(t)
    
    # Calculate the vector difference between the positions of the two satellites
    difference = position2 - position1

    # Calculate the distance between the two satellite positions and convert to meters
    distance = difference.distance().m

    return distance


def find_adjacent_orbit_sat( 
                            origin_sat, 
                            adj_plane, 
                            satellites_sorted_in_orbits,  
                            t
                            ):
    """
    Finds the satellite in the adjacent plane that is closest to an original satellite

    Args:
        origin_sat (object): the original satellite whose nearest neighboring satellite needs to be found
        adj_plane (int): the adjacent plane in which to search for the nearest satellite
        satellites_sorted_in_orbits (dict): list of satellites sorted by their orbit plane
        t (datetime): time corresponding to current satellite positions

    Returns:
        nearest_sat_in_adj_plane (object): satellite in the adjacent plane nearest to the original satellite
    """
    
    # Get the list of satellites in the specified adjacent plane
    adj_plane_sats = satellites_sorted_in_orbits[adj_plane]

    # Initialize variables to store the nearest satellite and set a minimum distance
    nearest_sat_in_adj_plane = -1
    min_distance = 1000000000000000

    # Iterate through satellites in the adjacent plane
    for i in range(len(adj_plane_sats)):
        
        # Calculate the distance between the original satellite and the current satellite in the adjacent plane
        distance = distance_between_two_satellites(origin_sat, adj_plane_sats[i], t)

        # Check if the calculated distance is smaller than both the current minimum distance and a threshold value
        if distance < min_distance and distance < 5016000:
            min_distance = distance # update the minimum distance
            nearest_sat_in_adj_plane = adj_plane_sats[i] # set the current adj. plane sat as the nearest to the original sat

    # if origin_sat.name == "STARLINK-1215":
    #     print(min_distance*1e-3, nearest_sat_in_adj_plane.name)

    # Return the name of the nearest satellite in the adjacent plane
    return nearest_sat_in_adj_plane.name.split(" ")[0] if nearest_sat_in_adj_plane != -1 else None


def get_current_isl_to_sats(
                            connectivity_matrix_row
                            ):
    
    satidx_list = []
    for idx, link in connectivity_matrix_row:
        if link==1:
            satidx_list.append(idx) 
    return satidx_list


def mininet_add_ISLs(
                        connectivity_matrix, 
                        satellites_sorted_in_orbits, 
                        satellites_by_name, 
                        satellites_by_index, 
                        isl_config, 
                        t
                    ):
    """
    Adds Inter-Satellite Links (ISLs) to the connectivity matrix

    Args:
        connectivity_matrix (list): 2D matrix representing the network connectivity between satellites, as well as ground stations
        satellites_sorted_in_orbits (dict): list of satellites sorted by their orbit plane
        satellites_by_name (dict): satellites sorted by name
        satellites_by_index (dict): satellites sorted by index
        isl_config (str): desired ISL configuration type
        t (datetime): time corresponding to current satellite positions

    Returns:
        connectivity_matrix (list): updated connectivity matrix, now including ISLs
    """

    # Get the number of orbits
    n_orbits = len(satellites_sorted_in_orbits)

    # Initialize the total number of satellites
    total_sat_now = 0

    # Check the ISL configuration (only one for the time being)
    if isl_config == "SAME_ORBIT_AND_GRID_ACROSS_ORBITS":
        
        # Iterate through each orbit
        for i in range(n_orbits):
           
            # Get the number of satellites in the current orbit
            n_sats_per_orbit = len(satellites_sorted_in_orbits[i])
            
            # Iterate through each satellite in the current orbit
            for j in range(n_sats_per_orbit):
                
                # Determine the index of the current satellite
                sat = total_sat_now + j
                current_sat_name = satellites_by_index[sat]
                current_sat = satellites_by_name[current_sat_name]

                # Determine the index of next satellite in same orbit
                sat_same_orbit = total_sat_now + ((j + 1) % n_sats_per_orbit)
                current_sat_same_orbit_name = satellites_by_index[sat_same_orbit]
                current_sat_same_orbit = satellites_by_name[current_sat_same_orbit_name]

                # Intra-orbit connection (Connection to all same orbit sats within threshold)
                if distance_between_two_satellites(current_sat, current_sat_same_orbit, t) < 5016000:
                    connectivity_matrix[sat][sat_same_orbit] = 1
                    connectivity_matrix[sat_same_orbit][sat] = 1
                
                # Inter-orbit connections
                # For the satellite in the next orbit
                sat_adjacent_orbit_1 = find_adjacent_orbit_sat(current_sat, (i + 1)%n_orbits, satellites_sorted_in_orbits, t)
                sat_adjacent_orbit_1_index = list(satellites_by_index.keys())[list(satellites_by_index.values()).index(sat_adjacent_orbit_1)]

                # For the satellite in the previous orbit
                sat_adjacent_orbit_2 = find_adjacent_orbit_sat(current_sat, (i - 1)%n_orbits, satellites_sorted_in_orbits, t)
                sat_adjacent_orbit_2_index = list(satellites_by_index.keys())[list(satellites_by_index.values()).index(sat_adjacent_orbit_2)]

                # Establishing connections
                connectivity_matrix[sat][sat_adjacent_orbit_1_index] = 1
                connectivity_matrix[sat_adjacent_orbit_1_index][sat] = 1
                connectivity_matrix[sat][sat_adjacent_orbit_2_index] = 1
                connectivity_matrix[sat_adjacent_orbit_2_index][sat] = 1

            # Update the current total number of satellites
            total_sat_now += n_sats_per_orbit

    # Simple distance based forward sat connections
    elif isl_config == "DISTANCE_BASED_SAME_AND_ACROSS_ORBITS":

        number_of_threads = 4
        max_isl_conn = 80
        
        # Setting maximum ISL length
        max_isl_search_length = int(5016000/2)

        numsats_per_orb = [len(orbs) for orbs in satellites_sorted_in_orbits]
        
        """
        OPTIMIZING CODESPACE STARTS
        """
        ######################## FINDING NEARBY SATS BASED ON THREADING ALLOCATION OF ORBIT #############################
        # Calculate number of pools and orbits per thread pool (for parallel execution)
        number_of_pools = n_orbits/number_of_threads
        num_of_orbits_per_pool = n_orbits/number_of_pools

        # Initialize list to store results for each pool
        current_sat_satellites_in_range = [[] for _ in range(int(number_of_pools+1))]

        # Create thread list
        thread_list = []
        count = 0

        # Divide same-orbit satellites into pools and create threads
        for pool in range(int(number_of_pools)):
            orb_index = int(num_of_orbits_per_pool)*pool
            same_orbit_sat_list = satellites_sorted_in_orbits[orb_index:orb_index+int(num_of_orbits_per_pool)]  # List of Skyfield type
            subsat_list = [sts for orb in same_orbit_sat_list for sts in orb]
            #print(subsat_list)
            total_sat_name_list = [sat.name.split(" ")[0] for sat in subsat_list]  # List of strs
            thread = threading.Thread(target=calc_distance_sat_sat_thread, args=(subsat_list, satellites_by_name, satellites_by_index, satellites_sorted_in_orbits, t, max_isl_search_length, current_sat_satellites_in_range[count]))
            thread_list.append(thread)
            count += 1

        # Start and join threads for parallel execution
        for thread in thread_list:
            thread.start()
        for thread in thread_list:
            thread.join()
        
        current_sat_satellites_in_range_flatten = [sats for sat_list in current_sat_satellites_in_range for sats in sat_list]
        for i in range(n_orbits):
            # Get sats and the number of satellites in the current orbit
            same_orbit_sat_list = satellites_sorted_in_orbits[i]    # List of Skyfield type
            total_sat_name_list = [sat.name.split(" ")[0] for sat in same_orbit_sat_list]  # List of strs
            adj_orb = [(i-1)%len(satellites_sorted_in_orbits), (i+1)%len(satellites_sorted_in_orbits)]
            adj_sat_ids = []
            adj_sat_ids_flatten = []
            for orb in adj_orb:
                adj_sats = satellites_sorted_in_orbits[orb]
                adj_sats_id = [list(satellites_by_index.keys())[list(satellites_by_index.values()).index(adj_sats[m].name.split(" ")[0])] for m in range(len(adj_sats))]
                adj_sat_ids.append(adj_sats_id)
                adj_sat_ids_flatten.extend(adj_sats_id)
            
            same_orbit_sats_id = [list(satellites_by_index.keys())[list(satellites_by_index.values()).index(s)] for s in total_sat_name_list]
            # print(i, same_orbit_sats_id)
            current_sat_satellites_in_range_flatten = [sats for sat_list in current_sat_satellites_in_range for sats in sat_list]
            for curr_sid_name in total_sat_name_list:  # Iterating over sats in the current ith orbit
                temp_isl_list = []
                sat = list(satellites_by_index.keys())[list(satellites_by_index.values()).index(curr_sid_name)]  # sat id based on name
                for sat_data in current_sat_satellites_in_range_flatten:
                    if sat_data:
                        if sat_data[2] == curr_sid_name:
                            temp_isl_list.append([sat_data[0],sat_data[1],sat])  # Taking all the isl connections generated from the thread for the current sat in the orbit (distance, neighbour_sat_id, current_sat_id)
                temp_isl_list = sorted(temp_isl_list, key=lambda x:x[0])  # sort from smallest to largest distance ISL
                #temp_isl_list = temp_isl_list[:max_isl_conn]
                same_counter = 0
                adj_flag1 = 0
                adj_flag2 = 0
                non_adj_counter = 0
                store_orbits = []   # Just stores the relevant orbits for each ISL connection per sat
                #print("Number of connections for " + curr_sid_name + "(" + str(sat) + ")" + " : " + str(len(temp_isl_list)))
                for sat_neighbour_index in temp_isl_list:
                    if sat_neighbour_index[1] in same_orbit_sats_id: # Takes two sats from the same orbit (j orbit)
                        if same_counter<2:
                            # Establishing connections
                            connectivity_matrix[sat][sat_neighbour_index[1]] = 1
                            connectivity_matrix[sat_neighbour_index[1]][sat] = 1
                            same_counter += 1
                            if same_counter==1:
                                store_orbits.append(i)
                            
                    elif sat_neighbour_index[1] in adj_sat_ids_flatten: # Takes two sats from the next near neighbours 

                        if sat_neighbour_index[1] in adj_sat_ids[0] and not adj_flag1:  # Takes one sat from j-1 orbit
                            # Establishing connections
                            connectivity_matrix[sat][sat_neighbour_index[1]] = 1
                            connectivity_matrix[sat_neighbour_index[1]][sat] = 1
                            adj_flag1 = 1
                            store_orbits.append(adj_orb[0])
                        elif sat_neighbour_index[1] in adj_sat_ids[1] and not adj_flag2:  # Takes one sat from j+1 orbit
                            # Establishing connections
                            connectivity_matrix[sat][sat_neighbour_index[1]] = 1
                            connectivity_matrix[sat_neighbour_index[1]][sat] = 1
                            adj_flag2 = 1
                            store_orbits.append(adj_orb[1])

                    else:  # If sats are neither in same orbit or the two closest orbit (k, l orbits)
                        neighbour_sat_orbit_id = get_orbit_from_sat(sat_neighbour_index[1], satellites_by_index, satellites_by_name, satellites_sorted_in_orbits)
                        if non_adj_counter<2 and neighbour_sat_orbit_id not in store_orbits:
                            # Establishing connections
                            connectivity_matrix[sat][sat_neighbour_index[1]] = 1
                            connectivity_matrix[sat_neighbour_index[1]][sat] = 1
                            non_adj_counter += 1
                            store_orbits.append(neighbour_sat_orbit_id)

                    if same_counter==2 and adj_flag1 and adj_flag2 and non_adj_counter==1:
                        # sat_count += 1
                        # print(sat_count, " Reached!")
                        break
        """
        OPTIMIZING CODESPACE ENDS
        """
        # Iterate through each orbit
        # sat_count = 0
        """
        for i in range(n_orbits):
            ######################## FINDING NEARBY SATS BASED ON THREADING INDIVIDUALLY FOR EACH ORBIT #############################
            # Get sats and the number of satellites in the current orbit
            same_orbit_sat_list = satellites_sorted_in_orbits[i]    # List of Skyfield type
            total_sat_name_list = [sat.name.split(" ")[0] for sat in same_orbit_sat_list]  # List of strs
            n_sats_per_orbit = len(same_orbit_sat_list)
                
            # Calculate number of pools and satellites per thread pool (for parallel execution)
            ####### (VERY IMPORTANT: FOR REAL TLES THIS SECTION WOULD CREATE ALOT OF PROBLEMS, SINCE EACH ORBIT HAS DIFFERENT NUMER OF SATS)
            number_of_pools = n_sats_per_orbit/number_of_threads
            num_of_sat_per_pool = n_sats_per_orbit/number_of_pools

            # Initialize list to store results for each pool
            current_sat_satellites_in_range = [[] for _ in range(int(number_of_pools+1))]

            # Create thread list
            thread_list = []
            count = 0

            # Divide same-orbit satellites into pools and create threads
            for pool in range(int(number_of_pools)):
                index = int(num_of_sat_per_pool)*pool
                if pool == int(number_of_pools)-1:      # Check for the last pool and assign all the remaining sats that were not able to be assigned because of num_of_sat_per_pool not a perfect integer
                    subsat_list = same_orbit_sat_list[index:]
                else:
                    subsat_list = same_orbit_sat_list[index:index+int(num_of_sat_per_pool)]
                thread = threading.Thread(target=calc_distance_sat_sat_thread, args=(subsat_list, satellites_by_name, satellites_by_index, satellites_sorted_in_orbits, t, max_isl_search_length, current_sat_satellites_in_range[count]))
                thread_list.append(thread)
                count += 1

            # Start and join threads for parallel execution
            for thread in thread_list:
                thread.start()
            for thread in thread_list:
                thread.join()

            # if i==32:
            #     print(current_sat_satellites_in_range[-1])
            ######################### THREADING DONE ########################
            # s = 0
            # for curr_sid_name in total_sat_name_list:
            #     c = 0
            #     for pool in current_sat_satellites_in_range:
            #         for sat_data in pool:
            #             if sat_data is not []:
            #                 d, sid, curr_sid = sat_data
            #                 if curr_sid == curr_sid_name:
            #                     c += 1
            #     print(i, s, curr_sid_name, c)
            #     s += 1

            adj_orb = [(i-1)%len(satellites_sorted_in_orbits), (i+1)%len(satellites_sorted_in_orbits)]
            adj_sat_ids = []
            adj_sat_ids_flatten = []
            for orb in adj_orb:
                adj_sats = satellites_sorted_in_orbits[orb]
                adj_sats_id = [list(satellites_by_index.keys())[list(satellites_by_index.values()).index(adj_sats[m].name.split(" ")[0])] for m in range(len(adj_sats))]
                adj_sat_ids.append(adj_sats_id)
                adj_sat_ids_flatten.extend(adj_sats_id)
            
            same_orbit_sats_id = [list(satellites_by_index.keys())[list(satellites_by_index.values()).index(s)] for s in total_sat_name_list]
            # print(i, same_orbit_sats_id)
            current_sat_satellites_in_range_flatten = [sats for sat_list in current_sat_satellites_in_range for sats in sat_list]
            for curr_sid_name in total_sat_name_list:  # Iterating over sats in the current ith orbit
                temp_isl_list = []
                sat = list(satellites_by_index.keys())[list(satellites_by_index.values()).index(curr_sid_name)]
                for sat_data in current_sat_satellites_in_range_flatten:
                    if sat_data:
                        if sat_data[2] == curr_sid_name:
                            temp_isl_list.append([sat_data[0],sat_data[1],sat])  # Taking all the isl connections generated from the thread for the current sat in the orbit (distance, neighbour_sat_id, current_sat_id)
                temp_isl_list = sorted(temp_isl_list, key=lambda x:x[0])  # sort from smallest to largest distance ISL
                #temp_isl_list = temp_isl_list[:max_isl_conn]
                same_counter = 0
                adj_flag1 = 0
                adj_flag2 = 0
                non_adj_counter = 0
                store_orbits = []   # Just stores the relevant orbits for each ISL connection per sat
                #print("Number of connections for " + curr_sid_name + "(" + str(sat) + ")" + " : " + str(len(temp_isl_list)))
                for sat_neighbour_index in temp_isl_list:
                    if sat_neighbour_index[1] in same_orbit_sats_id: # Takes two sats from the same orbit (j orbit)
                        if same_counter<2:
                            # Establishing connections
                            connectivity_matrix[sat][sat_neighbour_index[1]] = 1
                            connectivity_matrix[sat_neighbour_index[1]][sat] = 1
                            same_counter += 1
                            if same_counter==1:
                                store_orbits.append(i)
                            
                    elif sat_neighbour_index[1] in adj_sat_ids_flatten: # Takes two sats from the next near neighbours 

                        if sat_neighbour_index[1] in adj_sat_ids[0] and not adj_flag1:  # Takes one sat from j-1 orbit
                            # Establishing connections
                            connectivity_matrix[sat][sat_neighbour_index[1]] = 1
                            connectivity_matrix[sat_neighbour_index[1]][sat] = 1
                            adj_flag1 = 1
                            store_orbits.append(adj_orb[0])
                        elif sat_neighbour_index[1] in adj_sat_ids[1] and not adj_flag2:  # Takes one sat from j+1 orbit
                            # Establishing connections
                            connectivity_matrix[sat][sat_neighbour_index[1]] = 1
                            connectivity_matrix[sat_neighbour_index[1]][sat] = 1
                            adj_flag2 = 1
                            store_orbits.append(adj_orb[1])

                    else:  # If sats are neither in same orbit or the two closest orbit (k, l orbits)
                        neighbour_sat_orbit_id = get_orbit_from_sat(sat_neighbour_index[1], satellites_by_index, satellites_by_name, satellites_sorted_in_orbits)
                        if non_adj_counter<2 and neighbour_sat_orbit_id not in store_orbits:
                            # Establishing connections
                            connectivity_matrix[sat][sat_neighbour_index[1]] = 1
                            connectivity_matrix[sat_neighbour_index[1]][sat] = 1
                            non_adj_counter += 1
                            store_orbits.append(neighbour_sat_orbit_id)

                    if same_counter==2 and adj_flag1 and adj_flag2 and non_adj_counter==1:
                        # sat_count += 1
                        # print(sat_count, " Reached!")
                        break
        """
                #print("same counter --> " + str(same_counter) + "  flag 1 --> " + str(adj_flag1) + "  flag 2 --> " + str(adj_flag2))
            # print("Orbit " + str(i) + " : DONE")
        print("ISLs added to connectivity matrix")

    # Return the updated connectivity matrix
    return connectivity_matrix

def retrieve_GS_by_type(ground_stations, gs_type):
    """
    Retrieves a list of ground stations based on their type
    Gateway: 0 - connects satellite network to the internet
    Customer Terminal: 1 - connects end users to the satellite network
    Endpoint: 2 - Internet locations connected to gateways

    Args:
        ground_stations (dict): list of ground stations
        gs_type (int): type of ground station (0: gateway, 1: customer terminal, 2: endpoint)

    Returns:
        gs_list (list): list of ground stations of the specified type
    """
    
    # Initialize the list of ground stations
    gs_list = []

    # Iterate through the list of ground stations
    for gs in ground_stations:
        # Check if the ground station type matches the specified type
        if gs["type"] == gs_type:
            # Append the ground station to the list
            gs_list.append(gs)

    # Return the list of ground stations of the specified type
    return gs_list

def mininet_add_GSLs_parallel(
                              connectivity_matrix, 
                              satellites_by_name, 
                              satellites_by_index, 
                              ground_stations, 
                              number_of_threads, 
                              association_criteria, 
                              t, 
                              sat_config,
                              operator_name
                              ):
    """
    Adds Ground Station-Satellite Links (GSLs) to the connectivity matrix based on the desired association criteria

    Args: 
        connectivity_matrix (list): 2D matrix representing the network connectivity between satellites, as well as ground stations
        satellites_by_name (dict): satellites sorted by name
        satellites_by_index (dict): satellites sorted by index
        ground_stations (dict): list of ground stations
        number_of_threads (int): total number of threads (further clarification needed?)
        association_criteria (str): (??)
        t (datetime): time corresponding to current satellite positions
        sat_config (dict): constellation configuration from YAML
        operator_name (str): constellation/operator name
        
    Returns:
        connectivity_matrix (list): updated connectivity matrix, now including GSLs

    """
    
    # Retrieve maximum GSL length from config
    max_gsl_length_m = calc_max_gsl_length(sat_config, operator_name)

    # Check if max GSL length is valid
    if max_gsl_length_m == -1:
        if sat_config["Debug"] == 1:
            print ("[Mininet_add_GSLs] --- check the max GSL length variable ")
            return
        
    # Calculate number of pools and ground stations per thread pool (for parallel execution)
    number_of_pools = len(ground_stations)/number_of_threads
    num_of_gs_per_pool = len(ground_stations)/number_of_pools

    # Initialize list to store results for each pool
    ground_station_satellites_in_range = [[] for _ in range(int(number_of_pools+1))]

    # Create thread list
    thread_list = []
    count = 0

    # Divide ground stations into pools and create threads
    for i in range(0, len(ground_stations), int(num_of_gs_per_pool)):
        subgs_list = ground_stations[i:i+int(num_of_gs_per_pool)]
        thread = threading.Thread(target=calc_distance_gs_sat_thread, args=(subgs_list, satellites_by_name, satellites_by_index, t, max_gsl_length_m, ground_station_satellites_in_range[count]))
        thread_list.append(thread)
        count += 1

    # Start and join threads for parallel execution
    for thread in thread_list:
        thread.start()
    for thread in thread_list:
        thread.join()

    # Prepare temporary list for association criteria processing
    ground_station_satellites_in_range_temporary = []
    for list in ground_station_satellites_in_range:
        for ls in list:
            ground_station_satellites_in_range_temporary.append([[ls]])
    
    # Chooses a function to reconfigure the connectivity matrix to match the requested association criteria
    if association_criteria == "BASED_ON_DISTANCE_ONLY_MININET":
        connectivity_matrix = M_gs_sat_association_criteria_BasedOnDistance(connectivity_matrix, ground_station_satellites_in_range_temporary, ground_stations, len(satellites_by_index))
        return connectivity_matrix

    if association_criteria == "BASED_ON_DISTANCE_ONLY_MININET_ALAN":
        connectivity_matrix = M_gs_sat_no_association_criteria(connectivity_matrix, ground_station_satellites_in_range_temporary, len(satellites_by_index), satellites_by_index)
        return connectivity_matrix

    if association_criteria == "BASED_ON_LONGEST_ASSOCIATION_TIME":
        connectivity_matrix = M_gs_sat_association_criteria_MaxAssociationTime(connectivity_matrix, ground_station_satellites_in_range_temporary, ground_stations, len(satellites_by_index), satellites_by_index, satellites_by_name, max_gsl_length_m, t)
        return connectivity_matrix
    
    return -1

def M_gs_sat_no_association_criteria(
                                    connectivity_matrix, 
                                    all_gs_satellites_in_range, 
                                    num_of_satellites, 
                                    satellites_by_index
                                    ):
    """
    (??)

    Args:
        connectivity_matrix (list): 2D matrix representing the network connectivity between satellites, as well as ground stations
        all_gs_satellites_in_range (dict): list of tuples containing every ground station and the satellites in range of each of those ground stations
        num_of_satellites (int): total number of satellites
        satellites_by_index (dict): satellites sorted by index

    Returns:
        connectivitiy_matrix(list): updated connectivity matrix containing ??

    """

    ground_station_satellites_in_range = []

    for inrange_sat in all_gs_satellites_in_range:
        if len(inrange_sat[0]) != 0:
            ground_station_satellites_in_range.append(inrange_sat[0][0])

    for (az, distance_m, alt, sid, gr_id) in ground_station_satellites_in_range:
        connectivity_matrix[sid][num_of_satellites+0] = 1
        connectivity_matrix[num_of_satellites+0][sid] = 1

        print("best distance ",0, sid, satellites_by_index[sid], distance_m, az, alt)

    return connectivity_matrix

def last_visible_satellite(
                            ground_station, 
                            all_gs_satellites_in_range, 
                            satellites_by_index, 
                            satellites_by_name, 
                            max_gsl_length_m, 
                            t
                            ):
    """
    Determines the last visible satellite within a list of satellites in range of a given ground station

    Args:
        ground_station (object): relevant ground station
        all_gs_satellites_in_range (list): list of tuples containing every ground station and the satellites in range of each of those ground stations
        satellites_by_index (dict): satellites sorted by index
        satellites_by_name (dict): satellites sorted by name
        max_gsl_length_m (float): Maximum gs-sat link length (in meters)
        t (datetime): time corresponding to current satellite positions

    Returns:
        last_visible_satellite (tuple): the last visible satellite defined once by name and once by index

    """

    # Time step for each iteration
    step = 10       #in seconds

    # Extract current time and date information
    dt, leap_second = t.utc_datetime_and_leap_second()
    newscs = ((str(dt).split(" ")[1]).split(":")[2]).split("+")[0]
    date, timeN, zone = t.utc_strftime().split(" ")
    year, month, day = date.split("-")
    hour, minute, second = timeN.split(":")
    loggedTime = str(year)+","+str(month)+","+str(day)+","+str(hour)+","+str(minute)+","+str(newscs)

    # Initialize time variable
    ts = load.timescale()
    loop_t = ts.utc(int(year), int(month), int(day), int(hour), int(minute), float(newscs))

    # Identify visible satellites for the given ground station
    visible_sats = []
    for entry in all_gs_satellites_in_range:
        for val in entry:
            if len(val) > 0:
                if int(ground_station["gid"]) == int(val[0][2]):
                    visible_sats.append(val[0][1])

    # Count the number of visible satellites
    number_of_visible_sats = len(visible_sats)

    # Loop to find the last visible satellite
    cnt = 0
    while number_of_visible_sats > 1:
        # Update time for each iteration
        loop_t = ts.utc(int(year), int(month), int(day), int(hour), int(minute), float(newscs)+cnt)
        number_of_visible_sats = 0
        new_visible_sats = []
        for sat in visible_sats:
            satellite_name = satellites_by_index[sat]
            # Calculate distance between ground station and satellite
            distance = distance_between_ground_station_satellite(ground_station, satellites_by_name[satellite_name], loop_t)
            if distance <= max_gsl_length_m:
                number_of_visible_sats += 1
                new_visible_sats.append(sat)

        # Update the list of visible satellites and increment time
        visible_sats = new_visible_sats[:]
        cnt += step

    # Check if a single visible satellite is found
    if len(visible_sats) == 1:
        ground_station["next_update"] = loop_t.tt
        last_visible_satellite = (satellites_by_index[visible_sats[0]], visible_sats[0])
        return last_visible_satellite

    # Return -1 if no visible satellites are found
    return -1

def M_gs_sat_association_criteria_MaxAssociationTime(
                                                     connectivity_matrix, 
                                                     ground_station_satellites_in_range_temporary, 
                                                     ground_stations, 
                                                     num_of_satellites, 
                                                     satellites_by_index, 
                                                     satellites_by_name, 
                                                     max_gsl_length_m, 
                                                     t
                                                     ):
    """
    Configures the connectivity matrix to work with the association criteria that considers the maximum association time between nodes

    Args:
        connectivity_matrix (list): 2D matrix representing the network connectivity between satellites, as well as ground stations
        ground_station_satellites_in_range_temporary (??): ??
        ground_stations (dict): list of ground stations
        num_of_satellites (int): total number of satellites
        satellites_by_index (dict): satellites sorted by index
        satellites_by_name (dict): satellites sorted by name
        max_gsl_length_m (float): Maximum gs-sat link length (in meters)
        t (datetime): time corresponding to current satellite positions

    Returns:
        connectivitiy_matrix(list): updated connectivity matrix containing ??

    """
    
    for gs in ground_stations:
        if t.tt > gs["next_update"] or gs["next_update"] == "":
            chosen_satellite = last_visible_satellite(gs, ground_station_satellites_in_range_temporary, satellites_by_index, satellites_by_name, max_gsl_length_m, t)
            if chosen_satellite != -1:
                print("....... Current time = ", t.tt," GS#", gs["gid"], " is associated with SAT#", chosen_satellite[0]," which is named as ", chosen_satellite[1], ". The next uupdate time will be ", gs["next_update"])
                connectivity_matrix[num_of_satellites+gs["gid"]][chosen_satellite[1]] = 1
                connectivity_matrix[chosen_satellite[1]][num_of_satellites+gs["gid"]] = 1
                gs["sat_re_LAC"] = chosen_satellite[1]
            else:
                print(gs["gid"], -1)
        else:
            print("....... No updates = ", t.tt)
            connectivity_matrix[num_of_satellites+gs["gid"]][gs["sat_re_LAC"]] = 1
            connectivity_matrix[gs["sat_re_LAC"]][num_of_satellites+gs["gid"]] = 1
            continue

    return connectivity_matrix


def M_gs_sat_association_criteria_BasedOnDistance(
                                                    connectivity_matrix, 
                                                    all_gs_satellites_in_range, 
                                                    ground_stations, 
                                                    num_of_satellites, 
                                                    ):
    """
    Configures the connectivity matrix to work with the association criteria that considers the minimum distance between nodes

    Args:
        connectivity_matrix (list): 2D matrix representing the network connectivity between satellites, as well as ground stations
        all_gs_satellites_in_range (dict): list of tuples containing every ground station and the satellites in range of each of those ground stations
        ground_stations (dict): list of ground stations
        num_of_satellites (int): total number of satellites

    Returns:
        connectivitiy_matrix(list): updated connectivity matrix containing ??

    """
    
    gsl_snr = [0 for i in range(len(ground_stations))]
    gsl_latency = [0 for i in range(len(ground_stations))]
    ground_station_satellites_in_range = []

    for inrange_sat in all_gs_satellites_in_range:
        if len(inrange_sat[0]) != 0:
            ground_station_satellites_in_range.append(inrange_sat[0][0])

    # USE CASE 1 -- REMOVE for general run
    # chosen_sid_forAlan = -1
    # ######################################
    for gid in range(len(ground_stations)):
        chosen_sid = -1
        chosen_sid_list = []
        best_distance_m = 1000000000000000
        sat_wth_distance = {}
        for (distance_m, sid, gr_id) in ground_station_satellites_in_range:
            # print t.utc_strftime(), az, distance_m, alt, sid, gr_id
            if gid == gr_id:

                #$ Add top 4 closes sat-links for each GS
                sat_wth_distance[sid] = distance_m

                # if gid != 1: # USE CASE 1 -- REMOVE for general run
                # if distance_m < best_distance_m:
                #     chosen_sid = sid
                #     best_distance_m = distance_m
                # USE CASE 1 -- REMOVE for general run
                # if gid == 1:
                #     if sid == chosen_sid_forAlan:
                #         chosen_sid = sid
                #         best_distance_m = distance_m
                ######################################
        sat_wth_distance = dict(sorted(sat_wth_distance.items(),key=lambda x:x[1]))
        try:
            chosen_sid_list = list(sat_wth_distance.keys())[:4]   # Taking top 4 shortest sat-links to GS
        except:
            chosen_sid_list = list(sat_wth_distance.keys())       # Taking all sat-links to GS


        if chosen_sid_list:
            for sid_id in chosen_sid_list:
                connectivity_matrix[sid_id][num_of_satellites+gid] = 1
                connectivity_matrix[num_of_satellites+gid][sid_id] = 1

                gsl_snr[gid] = calc_gsl_snr_given_distance(best_distance_m)
                gsl_latency[gid] = best_distance_m/299792458            #speed of light

        # if chosen_sid != -1:
        #     # USE CASE 1 -- REMOVE for general run
        #     #if gid == 0:
        #     #    chosen_sid_forAlan = chosen_sid
        #     ######################################

        #     connectivity_matrix[chosen_sid][num_of_satellites+gid] = 1
        #     connectivity_matrix[num_of_satellites+gid][chosen_sid] = 1
        #     # print chosen_sid, gid, best_distance_m
        #     # if gid == 1:
        #     #     chosen_sid = chosen_sid_forAlan
        #     #     connectivity_matrix[chosen_sid][num_of_satellites+gid] = 1
        #     #     connectivity_matrix[num_of_satellites+gid][chosen_sid] = 1

        #     #print "best distance ",gid, chosen_sid, best_distance_m
        #     gsl_snr[gid] = calc_gsl_snr_given_distance(best_distance_m)
        #     gsl_latency[gid] = best_distance_m/299792458            #speed of light
        #     # print "best distance ",gid, chosen_sid, best_distance_m, gsl_latency[gid]

    return connectivity_matrix

# removed M_gs_sat_association_criteria_BasedOnDistance_alan, as it was not being used in any file or function

def calculate_link_characteristics_for_gsls_isls(
                                                connectivity_matrix, 
                                                satellites_by_index, 
                                                satellites_by_name, 
                                                ground_stations, 
                                                t
                                                ):
    """
    Calculates latency and throughput matrices for the network defined by the given connectivity matrix

    Args:
        connectivity_matrix (list): 2D matrix representing the network connectivity between satellites, as well as ground stations
        satellites_by_index (dict): satellites sorted by index
        satellites_by_name (dict): satellites sorted by name
        ground_stations (dict): list of ground stations
        t (datetime): time corresponding to current satellite positions
        
    Returns:
        latency_matrix (??): ??
        throughput_matrix (??): ??

    """
    
    # Initialize matrices for latency and throughput
    matrix_size = len(satellites_by_index)+len(ground_stations)
    latency_matrix = [[0.0 for _ in range(matrix_size)] for _ in range(matrix_size)]
    throughput_matrix = [[0.0 for _ in range(matrix_size)] for _ in range(matrix_size)]
    distance_matrix = [[0 for _ in range(matrix_size)] for _ in range(matrix_size)]
    
    # Define constants
    channel_bandwidth_downlink = 240
    channel_bandwidth_uplink = 60
    number_of_users_per_cell = 5.0
    density = 1.0/float(number_of_users_per_cell)

    # Loop through the connectivity matrix to calculate latency and throughput
    for i in range(len(connectivity_matrix)):
        for j in range(len(connectivity_matrix[i])):
            # ISL between two satellites
            if connectivity_matrix[i][j] == 1 and i < len(satellites_by_index) and j < len(satellites_by_index):
                distance_meters             = distance_between_two_satellites(satellites_by_name[str(satellites_by_index[i])], satellites_by_name[str(satellites_by_index[j])], t)
                distance_matrix[i][j]       = int(distance_meters)
                latency_matrix[i][j]        = ((distance_meters)/299792458.0)*1e3                                          #speed of light  (Units in ms)
                throughput_matrix[i][j]     = 500            #20Gbps

            # GSL between ground station and satellite
            if connectivity_matrix[i][j] == 1 and i >= len(satellites_by_index) and j < len(satellites_by_index):
                distance_meters             = distance_between_ground_station_satellite(ground_stations[i-len(satellites_by_index)], satellites_by_name[str(satellites_by_index[j])], t)
                distance_matrix[i][j]       = int(distance_meters)
                latency_matrix[i][j]        = ((distance_meters)/299792458.0)*1e3            #speed of light   (Units in ms)
                snr                         = calc_gsl_snr(satellites_by_name[str(satellites_by_index[j])], ground_stations[i-len(satellites_by_index)], t, distance_meters, "downlink")
                channel_width               = channel_bandwidth_downlink
                throughput_matrix[i][j]     = density*channel_width*(math.log(1+snr)/math.log(2))
                if throughput_matrix[i][j] > 500:
                    throughput_matrix[i][j] = 500

                # Additional check for specific conditions (further clarification?) [!!! As of now this part doesnt have significant effect !!!]
                if i-len(satellites_by_index) == 1:
                    snr                         = calc_gsl_snr(satellites_by_name[str(satellites_by_index[j])], ground_stations[i-len(satellites_by_index)], t, distance_meters, "downlink")
                    channel_width               = channel_bandwidth_downlink
                    throughput_matrix[i][j]     = density*channel_width*(math.log(1+snr)/math.log(2))
                    if throughput_matrix[i][j] > 500:
                        throughput_matrix[i][j] = 500
            
            # GSL between satellite and ground station
            if connectivity_matrix[i][j] == 1 and i < len(satellites_by_index) and j >= len(satellites_by_index):
                distance_meters             = distance_between_ground_station_satellite(ground_stations[j-len(satellites_by_index)], satellites_by_name[str(satellites_by_index[i])], t)
                distance_matrix[i][j]       = int(distance_meters)
                latency_matrix[i][j]        = ((distance_meters)/299792458.0)*1e3           #speed of light
                snr                         = calc_gsl_snr(satellites_by_name[str(satellites_by_index[i])], ground_stations[j-len(satellites_by_index)], t, distance_meters, "downlink")
                throughput_matrix[i][j]     = density*channel_bandwidth_downlink*(math.log(1+snr)/math.log(2))
                if throughput_matrix[i][j] > 500:
                    throughput_matrix[i][j] = 500

    # Return latency and throughput matrices
    return {
                "latency_matrix": latency_matrix,
                "throughput_matrix": throughput_matrix,
                "distance_matrix": distance_matrix
            }

###################################################
###################################################

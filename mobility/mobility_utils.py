from skyfield.api import wgs84, load
import math
import threading
import numpy as np
import sys
sys.path.append("../")
from link.link_utils import *
from DoTD import DoTD_History

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


def find_n_all_connected_sats(origin_sat, satellites_sorted_in_orbits, t):

    n_resp = 0
    s_e = []

    # Get the list of satellites in the specified adjacent plane
    for adj_plane in range(len(satellites_sorted_in_orbits)):

        adj_plane_sats = satellites_sorted_in_orbits[adj_plane]
        
        # Iterate through satellites in the adjacent plane
        for i in range(len(adj_plane_sats)):
            if origin_sat == adj_plane_sats[i]:
                continue
            # Calculate the distance between the original satellite and the current satellite in the adjacent plane
            distance = distance_between_two_satellites(origin_sat, adj_plane_sats[i], t)


            # Check if the calculated distance is smaller than both the current minimum distance and a threshold value
            if distance < 5016000:
                n_resp += 1
                s_e.append((adj_plane, i, adj_plane_sats[i]))

    return n_resp, s_e

def motif_find_m_se_e(satellites_sorted_in_orbits, satellites_by_name, satellites_by_index, t):
    
    min_connections_per_satellite = sys.maxsize
    S_e = []
    M = []
    e = None
    total_sat_now = 0
    for i in range(len(satellites_sorted_in_orbits)):
        
        n_sats_per_orbit = len(satellites_sorted_in_orbits[i])
        
        for j in range(n_sats_per_orbit):
            
            # Get information about the current satellite
            sat = total_sat_now + j
            current_sat = satellites_by_index[sat]
            current_sat = satellites_by_name[str(current_sat)]

            n_adj, s_e = find_n_all_connected_sats(current_sat, satellites_sorted_in_orbits, t)

            if n_adj < min_connections_per_satellite:
                min_connections_per_satellite = n_adj
                S_e = s_e
                e = (i, j, current_sat)

        total_sat_now += n_sats_per_orbit

    M = list( itertools.combinations(S_e, 2) )

    return M, e

def mininet_add_ISLs(
                        connectivity_matrix, 
                        satellites_sorted_in_orbits, 
                        satellites_by_name, 
                        satellites_by_index, 
                        isl_config, 
                        t,
                        M=None,
                        e=None,
                        n_seq=0,
                        dotd:DoTD_History = None
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

    ###
    elif isl_config == "DOTD":
        M = len(satellites_by_index)
        dotd.step()
        link_characterstics = calculate_link_characteristics_for_gsls_isls(satellites_by_index, satellites_by_name, t)
        latency_matrix    = np.array(link_characterstics["latency_matrix"])
        throughput_matrix = np.array(link_characterstics["throughput_matrix"]) / 1000
        distance_matrix   = latency_matrix / 1000 * 299792458.0
        is_connectd_matrix = (distance_matrix <= 10000*1000)
        latency_matrix = latency_matrix * is_connectd_matrix/1000
        S_max_t = dotd.S_max[dotd.t - 1]
        L_max_t = dotd.L_max[dotd.t - 1]
        S_max_t = max( S_max_t, np.max(throughput_matrix.ravel()) )
        L_max_t = max( L_max_t, np.max( latency_matrix.ravel() ) )
        dotd.S_max[dotd.t] = S_max_t
        dotd.L_max[dotd.t] = L_max_t


        A_k_x_t = np.zeros((M, M))
        for k in range(M):
            S_bar_k_x_t     = is_connectd_matrix[k, :] * throughput_matrix[k, :] / S_max_t
            L_bar_k_x_t     = is_connectd_matrix[k, :] * (1-latency_matrix[k, :]/L_max_t)
            phi_bar_i_j_tm1 = is_connectd_matrix[k, :] * dotd.phi_i_j_t[k, :, dotd.t-1] / 4
            A_k_x_t[k, :] = dotd.w1 * S_bar_k_x_t + dotd.w2*(L_bar_k_x_t) + dotd.w3*phi_bar_i_j_tm1 + is_connectd_matrix[k, :] * dotd.PI_i_t[:, dotd.t-1]

        phi_star_i_j = np.zeros((M, M))
        total_score = A_k_x_t * is_connectd_matrix
        while(np.sum(np.sum(total_score)) != 0):

            # Arg_value = np.max(total_score, 0)
            # Arg_ID    = np.argmax(total_score, 0)
            # sate_selection = np.argmax(Arg_value)
            # satellite = Arg_ID(sate_selection)
            sate_selection = np.argmax( np.max(total_score, 0) )
            satellite = np.argmax(total_score, 0)[sate_selection]
            if np.sum(phi_star_i_j[satellite,:]) < 4 and np.sum(phi_star_i_j[sate_selection,:]) < 4:
                phi_star_i_j[satellite,sate_selection] = is_connectd_matrix[satellite,sate_selection]
                phi_star_i_j[sate_selection,satellite] = is_connectd_matrix[sate_selection,satellite]

            total_score[satellite,sate_selection] = 0
            total_score[sate_selection,satellite] = 0
    
        for k in range(M):
            Obj_cap = phi_star_i_j[k,:]*throughput_matrix[k,:]/S_max_t
            Obj_Laten = phi_star_i_j[k,:]*(1-(latency_matrix[k,:]/L_max_t))
            Obj_LChurn = phi_star_i_j[k,:]*dotd.phi_i_j_t[k,:, dotd.t-1]/4
            dotd.PI_i_t[k, dotd.t] = np.sum(dotd.w1*Obj_cap+dotd.w2*Obj_Laten+dotd.w3*Obj_LChurn+phi_star_i_j[k,:]*dotd.PI_i_t[k, dotd.t-1])/4


        dotd.phi_i_j_t[:,:, dotd.t] = phi_star_i_j

        for i in range(M):
            for j in range(M):
                connectivity_matrix[i][j] = phi_star_i_j[i, j]

    elif isl_config == "CROSS_GRID":
        for i in range(len(satellites_sorted_in_orbits)):
            # Get the number of satellites in the current orbit
            n_sats_per_orbit = len(satellites_sorted_in_orbits[i])

            # Iterate through each satellite in the current orbit
            for j in range(n_sats_per_orbit):
                # Determine the index of the current satellite
                sat = total_sat_now + j
                
                # Get information about the current satellite
                current_sat = satellites_by_index[sat]
                current_sat = satellites_by_name[str(current_sat)]

                # Find satellites in adjacent orbits
                if not is_already_orbit_connected(sat, (i+1)%n_orbits, connectivity_matrix, satellites_sorted_in_orbits, satellites_by_index):
                    sat_adjacent_orbit_1 = find_adjacent_orbit_sat(current_sat, (i+1)%n_orbits, satellites_sorted_in_orbits, t)
                    sat_adjacent_orbit_1 = list(satellites_by_index.keys())[list(satellites_by_index.values()).index(str(sat_adjacent_orbit_1))]

                if not is_already_orbit_connected(sat, (i-1)%n_orbits, connectivity_matrix, satellites_sorted_in_orbits, satellites_by_index):
                    sat_adjacent_orbit_2 = find_adjacent_orbit_sat(current_sat, (i-1)%n_orbits, satellites_sorted_in_orbits, t)
                    sat_adjacent_orbit_2 = list(satellites_by_index.keys())[list(satellites_by_index.values()).index(str(sat_adjacent_orbit_2))]

                # Find satellites in adjacent orbits
                if not is_already_orbit_connected(sat, (i+2)%n_orbits, connectivity_matrix, satellites_sorted_in_orbits, satellites_by_index):
                    sat_adjacent_orbit_3 = find_adjacent_orbit_sat(current_sat, (i+2)%n_orbits, satellites_sorted_in_orbits, t)
                    sat_adjacent_orbit_3 = list(satellites_by_index.keys())[list(satellites_by_index.values()).index(str(sat_adjacent_orbit_3))]

                if not is_already_orbit_connected(sat, (i-2)%n_orbits, connectivity_matrix, satellites_sorted_in_orbits, satellites_by_index):
                    sat_adjacent_orbit_4 = find_adjacent_orbit_sat(current_sat, (i-2)%n_orbits, satellites_sorted_in_orbits, t)
                    sat_adjacent_orbit_4 = list(satellites_by_index.keys())[list(satellites_by_index.values()).index(str(sat_adjacent_orbit_4))]

                # Establish ISLs with satellites in adjacent orbits
                connectivity_matrix[sat][sat_adjacent_orbit_1] = 1
                connectivity_matrix[sat_adjacent_orbit_1][sat] = 1
                connectivity_matrix[sat][sat_adjacent_orbit_2] = 1
                connectivity_matrix[sat_adjacent_orbit_2][sat] = 1

                connectivity_matrix[sat][sat_adjacent_orbit_3] = 1
                connectivity_matrix[sat_adjacent_orbit_3][sat] = 1
                connectivity_matrix[sat][sat_adjacent_orbit_4] = 1
                connectivity_matrix[sat_adjacent_orbit_4][sat] = 1

            # Update the total number of satellites
            total_sat_now += n_sats_per_orbit

    elif isl_config == "MOTIF":
        
        M_seq = M[n_seq]
        first_ = M_seq[0]
        second_ = M_seq[1]
        first_node_orbit_offset, first_node_sat_offset = first_[0] - e[0], first_[1] - e[1]
        second_node_orbit_offset, second_node_sat_offset = second_[0] - e[0], second_[1] - e[1]

        while (first_node_orbit_offset == 0 and second_node_orbit_offset == 0) or (first_node_sat_offset == 0 and second_node_sat_offset == 0):
            n_seq += 1
            M_seq = M[n_seq]
            first_ = M_seq[0]
            second_ = M_seq[1]
            first_node_orbit_offset, first_node_sat_offset = first_[0] - e[0], first_[1] - e[1]
            second_node_orbit_offset, second_node_sat_offset = second_[0] - e[0], second_[1] - e[1]
        
        total_sat_now = 0
        for i in range(len(satellites_sorted_in_orbits)):

            n_sats_per_orbit = len(satellites_sorted_in_orbits[i])
            
            for j in range(n_sats_per_orbit):
            
                sat = total_sat_now + j
                current_sat = satellites_by_index[sat]
                current_sat = satellites_by_name[str(current_sat)]

                sat_adj_1_orbit = (i + first_node_orbit_offset) % n_orbits
                sat_adj_1_sat   = (j + first_node_sat_offset) % len(satellites_sorted_in_orbits[sat_adj_1_orbit])
                sat_adj_1 = satellites_sorted_in_orbits[sat_adj_1_orbit][sat_adj_1_sat]
                sat_adj_1 = list(satellites_by_index.keys())[list(satellites_by_index.values()).index(str(sat_adj_1.name))]

                sat_adj_2_orbit = (i + second_node_orbit_offset) % n_orbits
                sat_adj_2_sat   = (j + second_node_sat_offset) % len(satellites_sorted_in_orbits[sat_adj_2_orbit])
                sat_adj_2 = satellites_sorted_in_orbits[sat_adj_2_orbit][sat_adj_2_sat]
                sat_adj_2 = list(satellites_by_index.keys())[list(satellites_by_index.values()).index(str(sat_adj_2.name))]

                connectivity_matrix[sat][sat_adj_1] = 1
                connectivity_matrix[sat_adj_1][sat] = 1
                connectivity_matrix[sat][sat_adj_2] = 1
                connectivity_matrix[sat_adj_2][sat] = 1

            total_sat_now += n_sats_per_orbit

    # Return the updated connectivity matrix
    return connectivity_matrix

def is_already_orbit_connected(sat, dst_orbit, connectivity_matrix, satellites_sorted_in_orbits, satellites_by_index):

    n_sats_in_orbit = len(satellites_sorted_in_orbits[dst_orbit])
    for j in range(n_sats_in_orbit):
        sat_x = satellites_sorted_in_orbits[dst_orbit][j]
        sat_x = list(satellites_by_index.keys())[list(satellites_by_index.values()).index(sat_x.name.split(' ')[0])]
        if connectivity_matrix[sat][sat_x] == 1 or connectivity_matrix[sat_x][sat] == 1:
            return True
        
    return False

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

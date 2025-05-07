''''

SpaceNet Read Live TLEs/

AUTHOR:         Mohamed M. Kassem, Ph.D.
                University of Surrey

EDITOR:         Bruce Barbour
                Virginia Tech

DESCRIPTION:    This Python script supplies the utility functions for extracting and sorting constellation information, including info on its orbits and sorting of satellites 
                in their particular orbits. (FUTURE INTENTIONS - Identifying different shells of a constellation and sort within each shell)


CONTENTS:       TLE CONSTELLATION FUNCTIONS/ (STARTS AT 45)
                    get_orbital_planes(tle_filename, shell_num)
                    get_orbital_planes_classifications(tle_filename, constellation, number_of_orbits, number_of_sats_per_orbits, orbits_inclination)
                    sort_satellites_in_orbit(satellites_in_orbit, t)
'''

# =================================================================================== #
# ------------------------------- IMPORT PACKAGES ----------------------------------- #
# =================================================================================== #

import numpy as np
import jenkspy
from .lunar_dyn_utils import *
import mobility.lunar_dyn_utils as lunar_dyn
from .mobility_utils import *
import mobility.mobility_utils as mobil_utl

# =================================================================================== #
# -------------------------- TLE CONSTELLATION FUNCTIONS ---------------------------- #
# =================================================================================== #
def get_orbital_planes_classifications(
                                        tle_filename                : str, 
                                        constellation               : str, 
                                        number_of_orbits            : int, 
                                        shell_num                   : int, 
                                        orbits_inclination          : float,
                                        orbits_altitude             : float
                                      ) -> dict:
    """
    Arrange satellite orbital information from a TLE (Two-Line Element) file using Jenks Natural Breaks algorithm.

    Args:
        tle_filename (str):                 Path to the TLE file
        constellation (str):                Name of the satellite constellation
        number_of_orbits (int):             Number of orbital planes
        shell_num (int):                    Current shell number (used to classify orbit number alphabetically)
        orbits_inclination (float):         Inclination of the orbital planes
        orbits_altitude (float):            Altitude of the orbital planes

    Returns:
        dict:                               Dictionary containing arranged orbital information
                                            Key: Satellite name, Value: Tuple (Orbital number, Epoch, Inclination, RAAN, Eccentricity, Argument of Perigee, Mean Anomaly, Mean Motion, Shell number)
    """
    # tle_filename = "/home/suryaryan/t2t-plotting/dynamic-topology-generator/utils/starlink_tles/starlink_1727718467"
    # Initialize dictionaries as empty
    data_orbits                 = {}
    dump_orbital_data           = {"Epoch": [], "Satellites": [], "Inclination": [], "RAAN": [], "Mean anomaly": [], "ecc": [], "aop": [], "Mean motion": []}

    # Open TLE file in read mode
    tle_file = open(tle_filename, 'r')

    # Extract the contents of TLE file
    Lines = tle_file.readlines()

    # Defining thresholds
    thresh1, thresh2, thresh3, thresh4 = real_tle_filter(constellation, orbits_inclination, orbits_altitude)
    
    # First, we dump the TLE files into the dump_orbital_data variable; we read the three lines by three lines, and save satellite names, inclination and RAAN
    for i in range(0, len(Lines), 3):

        # First line TLE
        tle_first_line  = list([_f for _f in Lines[i+1].strip("\n").split(" ") if _f])

        # Second line TLE
        tle_second_line = list([_f for _f in Lines[i+2].strip("\n").split(" ") if _f])

        # Compute orbiting altitude
        tle_n           = float(tle_second_line[7]) * 2 * np.pi / 86400            # rad/s
        if constellation=='starlink':
            GM = 398600.435507   #km^3/s^2
            radius = 6378.137    #km
        elif constellation=='lunar':
            GM = get_value("GM")
            radius = get_value("radius")
        tle_a           = (GM / (tle_n ** 2)) ** (1. / 3.) - radius     # (altitude in km)

        # Inclination of constellation shell
        if  float(tle_second_line[2]) < (orbits_inclination + thresh1) and float(tle_second_line[2]) >= (orbits_inclination + thresh2) \
            and tle_a < (orbits_altitude + thresh3) and tle_a > (orbits_altitude + thresh4): 

            # Store TLE data in dump_orbital_data
            dump_orbital_data["Epoch"].append(tle_first_line[3])
            dump_orbital_data["Satellites"].append(Lines[i].strip())
            dump_orbital_data["Inclination"].append(tle_second_line[2])
            dump_orbital_data["RAAN"].append(float(tle_second_line[3]))
            dump_orbital_data["ecc"].append(tle_second_line[4])
            dump_orbital_data["aop"].append(tle_second_line[5])
            dump_orbital_data["Mean anomaly"].append(tle_second_line[6])
            dump_orbital_data["Mean motion"].append(tle_second_line[7])

    # Collect RAAN values in data dump
    list_of_values = [-1 for _ in range(len(dump_orbital_data["RAAN"]))]

    # Extract RAAN values for classification
    for i in range(0, len(dump_orbital_data["RAAN"])):
        list_of_values[i] = float(dump_orbital_data["RAAN"][i])
    
    # Use Jenks Natural Breaks classification to determine orbital planes
    breaks = jenkspy.jenks_breaks(list_of_values, n_classes=number_of_orbits)
    totalsatellites = 0
    # if shell_num==1:
    #     print(breaks)
    # Iterate over each determined natural break
    for b in range(1, len(breaks)):
        
        # Define the bounds from Natural Breaks
        upperBound_of_class = float(breaks[b])
        lowerBound_of_class = float(breaks[b-1])

        # Initialize variables
        class_num = b-1
        count_sats_per_orbit = 0

        # Iterate the satellite data and breaks in RAAN to arrange the satellites in their respective orbits
        for i, j in zip(list(range(len(dump_orbital_data["Satellites"]))), list(range(len(dump_orbital_data["RAAN"])))):
            #print(dump_orbital_data["Inclination"][j])
            # Only for the first break
            if b == 1:

                # Check if RAAN falls within the specified bounds
                if float(dump_orbital_data["RAAN"][j]) <= upperBound_of_class and float(dump_orbital_data["RAAN"][j]) >= lowerBound_of_class:
                    
                    # Store orbital information in data_orbits dictionary
                    data_orbits[dump_orbital_data["Satellites"][i]] = (class_num, dump_orbital_data["Epoch"][j], dump_orbital_data["Inclination"][j], dump_orbital_data["RAAN"][j], dump_orbital_data["ecc"][j], dump_orbital_data["aop"][j], dump_orbital_data["Mean anomaly"][j], dump_orbital_data["Mean motion"][j], shell_num+1)
                    
                    # Count satellites in orbit
                    count_sats_per_orbit += 1

                    # print(dump_orbital_data["Satellites"][i], dump_orbital_data["RAAN"][j])

            else:

                # Check if RAAN falls within the specified bounds
                if float(dump_orbital_data["RAAN"][j]) <= upperBound_of_class and float(dump_orbital_data["RAAN"][j]) > lowerBound_of_class:
                    
                    # Store orbital information in data_orbits dictionary
                    data_orbits[dump_orbital_data["Satellites"][i]] = (class_num, dump_orbital_data["Epoch"][j], dump_orbital_data["Inclination"][j], dump_orbital_data["RAAN"][j], dump_orbital_data["ecc"][j], dump_orbital_data["aop"][j], dump_orbital_data["Mean anomaly"][j], dump_orbital_data["Mean motion"][j], shell_num+1)
                    
                    # Count satellites in orbit
                    count_sats_per_orbit += 1

        # print("Num of Sats ----------------", count_sats_per_orbit)
                    
        # Count the total number of satellites
        totalsatellites += count_sats_per_orbit

    # Return the collected orbital information separated by orbit
    return data_orbits


def sort_satellites_in_orbit(
                                satellites_in_orbit : list, 
                                t                   : float
                            ) -> list:
    """
    Sorts a list of satellites in an orbit based on their distance from the first indexed satellite.

    Args:
        satellites_in_orbit (list): List of satellites in orbit
        t (float):                  Current simulation time

    Returns:
        list:                       List of sorted satellites based on their distance from the first indexed satellite
    """

    # Initialize lists as empty
    visited_sats = []
    sorted_sats = []

    # Select the first satellite in the orbit as the starting point
    first_sat = satellites_in_orbit[0]

    # Add the first satellite to the lists
    sorted_sats.append(first_sat)
    visited_sats.append(first_sat.name)

    # Change epoch type based on main_body
    if get_main_body_str(first_sat) != 'Earth':
        t = first_sat.epoch   #changes the type to astropy Time object
        distance_between_two_satellites = lunar_dyn.distance_between_two_satellites
    else:
        distance_between_two_satellites = mobil_utl.distance_between_two_satellites

    # Iterate through the satellites in the orbit and find the next corresponding satellite with the minimum distance
    for _ in range(len(satellites_in_orbit)):
        
        # Prevent duplicates
        next_hop = -1

        # Minimum distance set to infinity for checking
        min_distance = float('inf')
        
        # Find the next satellite with the minimum distance that hasn't been visited
        for sat in satellites_in_orbit:
            if sat.name not in visited_sats:
                distance = distance_between_two_satellites(first_sat, sat, t)
                if distance < min_distance:
                    next_hop = sat
                    min_distance = distance

        # If a valid next satellite is found, update the current satellite and add it to the sorted list
        if next_hop != -1:
            first_sat = next_hop
            visited_sats.append(first_sat.name)
            sorted_sats.append(first_sat)

    # Return the sorted list of satellites in orbit
    return sorted_sats


def real_tle_filter(operator_name, orbits_inclination, orbits_altitude):
    """
    INPUT:  operator_name (str)        : Name of the constellation (SUPPORTS: starlink, lunar)
            orbits_inclination (float) : Mean inclination of the shell (in degrees)
            orbits_altitude (float)    : Mean altitude of the shell (in kms)


    OUTPUT:  thresh1 : Inclination lower bound
             thresh2 : Inclination upper bound
             thresh3 : ALtitude lower bound
             thresh4 : Altitude upper bound

    """

    if operator_name=='starlink':
        if orbits_inclination == 53.2 and orbits_altitude == 540:   # ref Starlink FCC
            thresh1 = 0.1
            thresh2 = -0.9
            thresh3 = 7.1524
            thresh4 = -9.0524
        elif orbits_inclination == 97.6 and orbits_altitude == 560:   # ref Starlink FCC
            thresh1 = 0.1
            thresh2 = -0.9
            thresh3 = 3.0624
            thresh4 = 2.0524
        else:
            thresh1 = 0.1
            thresh2 = -0.1
            thresh3 = 1
            thresh4 = -1


    elif operator_name=='lunar':
        thresh1 = 1
        thresh2 = -1
        thresh3 = 1
        thresh4 = -1
    
    else:
        thresh1 = 0.1
        thresh2 = -0.1
        thresh3 = 1
        thresh4 = -1

    return thresh1, thresh2, thresh3, thresh4

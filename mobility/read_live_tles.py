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
from .mobility_utils import *

# =================================================================================== #
# -------------------------- TLE CONSTELLATION FUNCTIONS ---------------------------- #
# =================================================================================== #
def get_orbital_planes_classifications(
                                        tle_filename                : str, 
                                        constellation               : str, 
                                        number_of_orbits            : int, 
                                        number_of_sats_per_orbits   : int, 
                                        orbits_inclination          : float,
                                        orbits_altitude             : float
                                      ) -> dict:
    """
    Arrange satellite orbital information from a TLE (Two-Line Element) file using Jenks Natural Breaks algorithm.

    Args:
        tle_filename (str):                 Path to the TLE file
        constellation (str):                Name of the satellite constellation
        number_of_orbits (int):             Number of orbital planes
        number_of_sats_per_orbits (int):    Number of satellites per orbital plane (unused)
        orbits_inclination (float):         Inclination of the orbital planes
        orbits_altitude (float):            Altitude of the orbital planes

    Returns:
        dict:                               Dictionary containing arranged orbital information
                                            Key: Satellite name, Value: Tuple (Orbital number, Epoch, Inclination, RAAN, Eccentricity, Argument of Perigee, Mean Anomaly, Mean Motion)
    """
    # tle_filename = "/home/suryaryan/t2t-plotting/dynamic-topology-generator/utils/starlink_tles/starlink_1727718467"
    # Initialize dictionaries as empty
    data_orbits                 = {}
    dump_orbital_data           = {"Epoch": [], "Satellites": [], "Inclination": [], "RAAN": [], "Mean anomaly": [], "ecc": [], "aop": [], "Mean motion": []}

    # Open TLE file in read mode
    tle_file = open(tle_filename, 'r')

    # Open Output TLE save file in write mode
    path_segments = tle_filename.split('/')
    tle_savefilename = "/home/suryaryan/t2t-plotting/dynamic-topology-generator/utils/analysis/extracted_"+path_segments[-1]
    tle_savefile = open(tle_savefilename, 'w')

    # Extract the contents of TLE file
    Lines = tle_file.readlines()

    # First, we dump the TLE files into the dump_orbital_data variable; we read the three lines by three lines, and save satellite names, inclination and RAAN
    for i in range(0, len(Lines), 3):

        # First line TLE
        tle_first_line  = list([_f for _f in Lines[i+1].strip("\n").split(" ") if _f])

        # Second line TLE
        tle_second_line = list([_f for _f in Lines[i+2].strip("\n").split(" ") if _f])

        # Compute orbiting altitude
        tle_n           = float(tle_second_line[7]) * 2 * np.pi / 86400            # rad/s
        tle_a           = (398600.435507 / (tle_n ** 2)) ** (1. / 3.) - 6378.137     # (altitude in km)

        # Inclination of constellation shell
        if  float(tle_second_line[2]) < (orbits_inclination + 0.1) and float(tle_second_line[2]) >= (orbits_inclination - 0.9) \
            and tle_a < (orbits_altitude + 7.1524) and tle_a > (orbits_altitude - 9.0524): 
        #if (tle_a > (orbits_altitude + 7.1524)) and float(tle_second_line[2]) < (orbits_inclination + 8):    # tle : (1727475306, 1727718467, 1730840419)  threshold : (7.1524, 7.0330, 7.0080)  alt-540            ":

            # Store TLE data in dump_orbital_data
            dump_orbital_data["Epoch"].append(tle_first_line[3])
            dump_orbital_data["Satellites"].append(Lines[i].strip())
            dump_orbital_data["Inclination"].append(tle_second_line[2])
            dump_orbital_data["RAAN"].append(float(tle_second_line[3]))
            dump_orbital_data["ecc"].append(tle_second_line[4])
            dump_orbital_data["aop"].append(tle_second_line[5])
            dump_orbital_data["Mean anomaly"].append(tle_second_line[6])
            dump_orbital_data["Mean motion"].append(tle_second_line[7])

            # Storing the TLEs of all the selected sats in a file
            tle_savefile.writelines([Lines[i], Lines[i+1], Lines[i+2]]) 

    # Collect RAAN values in data dump
    list_of_values = [-1 for _ in range(len(dump_orbital_data["RAAN"]))]

    # Extract RAAN values for classification
    for i in range(0, len(dump_orbital_data["RAAN"])):
        list_of_values[i] = float(dump_orbital_data["RAAN"][i])
    
    # Use Jenks Natural Breaks classification to determine orbital planes
    breaks = jenkspy.jenks_breaks(list_of_values, n_classes=number_of_orbits)
    totalsatellites = 0

    # Iterate over each determined natural break
    for b in range(1, len(breaks)):
        
        # Define the bounds from Natural Breaks
        upperBound_of_class = float(breaks[b])
        lowerBound_of_class = float(breaks[b-1])

        # Initialize variables
        class_num = b-1
        count_sats_per_orbit = 0
        # print("Class -------------------- "+str(class_num))

        # Iterate the satellite data and breaks in RAAN to arrange the satellites in their respective orbits
        for i, j in zip(list(range(len(dump_orbital_data["Satellites"]))), list(range(len(dump_orbital_data["RAAN"])))):
           
            # Only for the first break
            if b == 1:

                # Check if RAAN falls within the specified bounds
                if float(dump_orbital_data["RAAN"][j]) <= upperBound_of_class and float(dump_orbital_data["RAAN"][j]) >= lowerBound_of_class:
                    
                    # Store orbital information in data_orbits dictionary
                    data_orbits[dump_orbital_data["Satellites"][i]] = (class_num, dump_orbital_data["Epoch"][j], dump_orbital_data["Inclination"][j], dump_orbital_data["RAAN"][j], dump_orbital_data["ecc"][j], dump_orbital_data["aop"][j], dump_orbital_data["Mean anomaly"][j], dump_orbital_data["Mean motion"][j])#Satellite name: (Inclination, RAAN, orbital number)
                    
                    # Count satellites in orbit
                    count_sats_per_orbit += 1

                    # print(dump_orbital_data["Satellites"][i], dump_orbital_data["RAAN"][j])

            else:

                # Check if RAAN falls within the specified bounds
                if float(dump_orbital_data["RAAN"][j]) <= upperBound_of_class and float(dump_orbital_data["RAAN"][j]) > lowerBound_of_class:
                    
                    # Store orbital information in data_orbits dictionary
                    data_orbits[dump_orbital_data["Satellites"][i]] = (class_num, dump_orbital_data["Epoch"][j], dump_orbital_data["Inclination"][j], dump_orbital_data["RAAN"][j], dump_orbital_data["ecc"][j], dump_orbital_data["aop"][j], dump_orbital_data["Mean anomaly"][j], dump_orbital_data["Mean motion"][j])#Satellite name: (Inclination, RAAN, orbital number)
                    
                    # Count satellites in orbit
                    count_sats_per_orbit += 1

        # print("Num of Sats ----------------", count_sats_per_orbit)
                    
        # Count the total number of satellites
        totalsatellites += count_sats_per_orbit

    """Debugging Zone"""

    # # Print total satellites for checking before completing sim
    # #ggs = [((398600.435507 / ((float(dump_orbital_data["Mean motion"][i]) * 2 * np.pi / 86400) ** 2)) ** (1. / 3.) - 6378.137) for i, data in enumerate(dump_orbital_data["Inclination"]) if float(data) > 40 and float(data) < 45 and ((398600.435507 / ((float(dump_orbital_data["Mean motion"][i]) * 2 * np.pi / 86400) ** 2)) ** (1. / 3.) - 6378.137)>0]
    # ggs = [float(data) for i, data in enumerate(dump_orbital_data["Inclination"]) if float(data) > 40 and float(data) < 45 and ((398600.435507 / ((float(dump_orbital_data["Mean motion"][i]) * 2 * np.pi / 86400) ** 2)) ** (1. / 3.) - 6378.137)>0]
    # print(".......... No. of Sat Nodes: ", totalsatellites, "  Total number of sats in desired inclination: ", len(ggs))
    
    # raan_mean = range(0,360,5)
    # RAAN_sats = {i:[] for i in raan_mean}
    # anom_sats = {i:[] for i in raan_mean}
    # for raan in raan_mean:
    #     anomalies = []
    #     for k, data in enumerate(dump_orbital_data["RAAN"]):
    #         if float(data) > raan-1.7 and float(data) < raan+1.7:
    #             RAAN_sats[raan].append(float(data))
    #             anomalies.append(float(dump_orbital_data["Mean anomaly"][k]))

    #     sat_list = sorted(anomalies)
    #     anom_mean = range(int(sat_list[0]),360+int(sat_list[0]),20)
    #     for anom_data in anom_mean:
    #         for sat_anoms in sat_list:
    #             if float(sat_anoms) > anom_data-4.5 and float(sat_anoms) < anom_data+4.5:
    #                 anom_sats[raan].append(float(sat_anoms))

    # print(".......... No. of Sat Nodes: ", totalsatellites, "  Total number of sats after anomaly filter: ", sum([len(R_sats) for R_sats in RAAN_sats.values()]))
    # print(anom_sats[raan_mean[0]])

    # import matplotlib.pyplot as plt
    # plt.hist(ggs)
    # plt.title("Number of satellites from tle: " + str(totalsatellites))
    # plt.show()

    """Debugging Zone Ends!"""

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
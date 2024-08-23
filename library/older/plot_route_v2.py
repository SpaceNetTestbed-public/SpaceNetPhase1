import matplotlib.pyplot as plt
from skyfield.api import load, EarthSatellite
import re
import numpy as np
from datetime import datetime, timezone
from mpl_toolkits.basemap import Basemap







# ================================================================================================
# SCRIPT CONTROL
# ================================================================================================
time_index = 0
plot_only_optimal = False
gs_filepath = open('/home/barbourbruce/dynamic-topology-generator/utils/gs_files/gs_default.txt', 'r')
tle_file = open('/home/barbourbruce/dynamic-topology-generator/utils/starlink_tles/starlink_1720559136', 'r')
optimal_route_filepath = '/home/barbourbruce/dynamic-topology-generator/output/optimal_routes/starlink/best_path_2024_07_09.txt'
node_indices_filepath = '/home/barbourbruce/dynamic-topology-generator/output/node_indices/starlink/nodeindex_1720559136.txt'















# ================================================================================================
# INITIALIZER
# ================================================================================================
ts = load.timescale()
sats_from_tle_dict                  = {}
gs_dict                             = {}
plot_sat_alias_dict                 = {}
plot_sat_alias_to_index_dict        = {}
node_alias_to_index_topology_dict   = {}
node_index_to_alias_topology_dict   = {}
node_info_topology_at_t             = {}
optimal_routes                      = []
dt_hist                             = []
gs_alias_list                       = ["CT", "GS", "GW", "IE"]
total_num_sat                       = 0
total_num_gs                        = 0

# ================================================================================================
# FILE PARSING - TLE DICT
# ================================================================================================
tle_lines = tle_file.readlines()
for i in range(0, len(tle_lines), 3):
    tle_sat_name                    = tle_lines[i]
    tle_first_line                  = tle_lines[i+1]
    tle_second_line                 = tle_lines[i+2]
    tle_sat_obj                     = EarthSatellite(tle_first_line, tle_second_line)
    tle_sat_obj.name                = re.match(r"STARLINK-\d+", tle_sat_name).group(0)
    sats_from_tle_dict[tle_sat_obj.name] = tle_sat_obj

# ================================================================================================
# FILE PARSING - GROUND STATION
# ================================================================================================
gs_lines = gs_filepath.readlines()
total_num_gs_in_gsfile = len(gs_lines)
for i in range(0, len(gs_lines), 1):
    gs_info = gs_lines[i].split(",")
    gs_dict['GS-'+gs_info[0]] = (gs_info[1], float(gs_info[3]), float(gs_info[2]))

# ================================================================================================
# FILE PARSING - NODE INDEXING DICT
# ================================================================================================
with open(node_indices_filepath, 'r') as node_indices_file:
    for node_assignment in node_indices_file:
        node_index, node_alias  = node_assignment.split(":")
        node_index              = int(node_index)
        node_alias              = node_alias[:-1]
        if "GS" in node_alias:
            node_alias_to_index_topology_dict[node_alias] = node_index
            node_index_to_alias_topology_dict[node_index] = node_alias
            total_num_gs += 1
        else:
            node_alias_to_index_topology_dict[node_alias] = node_index
            node_index_to_alias_topology_dict[node_index] = node_alias
            total_num_sat += 1
if total_num_gs > total_num_gs_in_gsfile:
    raise ValueError("The number of ground stations observed in the topology is greater than in the provided ground station file!")

# ================================================================================================
# FILE PARSING - OPTIMAL PATH ASSIGNMENT
# ================================================================================================
with open(optimal_route_filepath, 'r') as optimal_path_file:
    for route_line in optimal_path_file:

        # Separate datetime and route information
        dt, route_info          = route_line.split(": ", 1)

        # Build time history
        dt_info                 = dt.strip('()').split("_")
        yr, mon, day, hr, min   = map(int, dt_info[:5])
        sec                     = float(dt_info[5])
        dt_hist.append(datetime(yr, mon, day, hr, min, int(sec), int((sec - int(sec)) * 1000000), tzinfo=timezone.utc))

        # Routing information
        route_node_indices      = route_info.split(", ")
        route_node_indices[-1]  = route_node_indices[-1][:-1]   # Removes the next line command
        
        # Assign the corresponding satellite alias with the route node indices
        optimal_route_at_epoch  = []
        for route_node_index in route_node_indices:
            optimal_route_at_epoch.append(node_index_to_alias_topology_dict[int(route_node_index)])
        
        # Append to complete list
        optimal_routes.append(optimal_route_at_epoch)

# ================================================================================================
# DEFINE CURRENT TIME INDEX
# ================================================================================================
t_current   = ts.from_datetime(dt_hist[time_index])

# ================================================================================================
# CREATE DICTIONARY WITH SATELLITE GEODETIC POSITION
# ================================================================================================
for topology_node_alias, topology_node_index in node_alias_to_index_topology_dict.items():

    # Check that it's a satellite
    if not any(gs_type in topology_node_alias for gs_type in gs_alias_list):

        # Extract satellite object based on node alias
        sat_node            = sats_from_tle_dict[topology_node_alias]

        # Propagate satellite object to t
        sat_node_at_t       = sat_node.at(t_current)

        # Compute longitude/latitude of satellite object at t
        lon_at_t, lat_at_t  = sat_node_at_t.subpoint().longitude.degrees, sat_node_at_t.subpoint().latitude.degrees
        
        # Add to dictionary
        node_info_topology_at_t[sat_node.name] = (node_alias_to_index_topology_dict[sat_node.name], lon_at_t, lat_at_t)

    # And if it's a ground station
    else:

        # Add to dictionary
        node_info_topology_at_t[topology_node_alias] = gs_dict[topology_node_alias]
        
# ================================================================================================
# PLOTTING
# ================================================================================================
plt.figure(figsize=(10, 5))

# PLOT BASEMAP
m = Basemap(projection='cyl', llcrnrlat=-60, urcrnrlat=60, llcrnrlon=-180, urcrnrlon=180, resolution='c')
m.drawcoastlines()
m.drawcountries()
m.fillcontinents(color='lightgray')

# PLOT ALL SATELLITE NODES IN TOPOLOGY
if not plot_only_optimal:
    for node_alias, node_info in node_info_topology_at_t.items():

        # Extract information
        node_assigned_alias     = node_info[0]
        node_lon, node_lat      = node_info[1:]
        
        # Plot satellite node as a regular scatter point with label
        if not "GS" in node_alias:
            plt.scatter(node_lon, node_lat, s=1, marker="o", facecolors='none', edgecolors='black', zorder=20)
            plt.text(node_lon, node_lat-0.5, node_assigned_alias, fontsize=7, zorder=100)
    
# PLOT GROUND STATIONS OF INTEREST
optimal_route_at_t          = optimal_routes[time_index]
gs0                         = node_info_topology_at_t[optimal_route_at_t[0]]
gs1                         = node_info_topology_at_t[optimal_route_at_t[-1]]
plt.scatter(gs0[1], gs0[2], s=100, marker='x', linewidth=2, c='g', zorder=20, label=gs0[0])
plt.scatter(gs1[1], gs1[2], s=100, marker='x', linewidth=2, c='b', zorder=20, label=gs1[0])
plt.text(gs0[1], gs0[2]-0.5, gs0[0], fontsize=7, zorder=100)
plt.text(gs1[1], gs1[2]-0.5, gs1[0], fontsize=7, zorder=100)

# PLOT OPTIMAL ROUTE
optimal_lon = np.array([0., ] * len(optimal_route_at_t))
optimal_lat = np.array([0., ] * len(optimal_route_at_t))
for indx, optimal_node in enumerate(optimal_route_at_t):
    optimal_node_info   = node_info_topology_at_t[optimal_node]
    optimal_lon[indx]   = optimal_node_info[1]
    optimal_lat[indx]   = optimal_node_info[2]
    if plot_only_optimal and not "GS" in optimal_node:
        plt.text(optimal_lon[indx] , optimal_lat[indx]-0.5, optimal_node_info[0], fontsize=6, zorder=100)
plt.plot(optimal_lon[:], optimal_lat[:], '--*', linewidth=1.5, c='r', zorder=12)

# PLOT INFORMATION
plt.title('FW Algorithm: '+str(total_num_sat)+' nodes (time: '+str(dt_hist[time_index])+') (hops='+str(len(optimal_route_at_t)-1)+')')
plt.legend(loc='upper left')
plt.xlabel('Longitude')
plt.ylabel('Latitude')
# plt.show()
plt.savefig('/home/barbourbruce/dynamic-topology-generator/test_plot.pdf')
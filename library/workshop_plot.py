import matplotlib
#matplotlib.use('tkagg')
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from skyfield.api import load, EarthSatellite
import re
import numpy as np
from datetime import datetime, timezone
from mpl_toolkits.basemap import Basemap






# ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

# ================================================================================================
# >>> SCRIPT CONTROL - EDIT HERE <<<
# ================================================================================================
time_index              = 0
plot_only_optimal       = False
plot_in_3D              = True
lon0_3d                 = -40
lat0_3d                 = 0
gs_filepath             = open('/home/spacenet/t2t-plotting/dynamic-topology-generator/output/real/NYC_LA/terrestrial_info/terrestrial_1721256111.txt', 'r')
tle_file                = open('/home/spacenet/t2t-plotting/dynamic-topology-generator/utils/starlink_tles/starlink_1721256111', 'r')
optimal_route_filepath  = '/home/spacenet/t2t-plotting/dynamic-topology-generator/output/real/NYC_LA/optimal_routes/starlink/best_path_2024_07_17_22_41_51.0.txt'
node_indices_filepath   = '/home/spacenet/t2t-plotting/dynamic-topology-generator/output/real/NYC_LA/node_indices/starlink/nodeindex_1721256111.txt'

# ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

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
    gs_info = gs_lines[i][:-1].split(":")
    gs_dict[gs_info[1]] = (int(gs_info[0]), gs_info[2], float(gs_info[4]), float(gs_info[3]))

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
        node_info_topology_at_t[sat_node.name] = (topology_node_alias, lon_at_t, lat_at_t)

    # And if it's a ground station
    else:

        # Add to dictionary
        node_info_topology_at_t[topology_node_alias] = gs_dict[topology_node_alias][1:]
        
# ================================================================================================
# PLOTTING
# ================================================================================================
fig = plt.figure()
font = {'family' : 'monospace', 
        'size' : 12}
plt.rc('font', **font)

# PLOT BASEMAP
if not plot_in_3D:
    m = Basemap(projection='cyl', llcrnrlat=-60, urcrnrlat=60, llcrnrlon=-180, urcrnrlon=180, resolution='c')
else:
    m0 = Basemap(projection='ortho', lat_0=lat0_3d, lon_0=lon0_3d, resolution=None)
    m = Basemap(projection='ortho', lat_0=lat0_3d, lon_0=lon0_3d, llcrnrx=-m0.urcrnrx/1.75, llcrnry=-m0.urcrnry/1.75, urcrnrx=m0.urcrnrx/1.75, urcrnry=m0.urcrnry/1.75, resolution='c')
    #m = Basemap(projection='ortho', lat_0=lat0_3d, lon_0=lon0_3d, llcrnrx=-m0.urcrnrx/1.75, llcrnry=-m0.urcrnry/10.75, urcrnrx=m0.urcrnrx/1.75, urcrnry=m0.urcrnry/1.75, resolution='c')
m.drawcoastlines()
m.drawcountries()
m.fillcontinents(color='lightgray', lake_color='white')
m.drawmapboundary(fill_color='white')

# PLOT ALL SATELLITE NODES IN TOPOLOGY
if not plot_only_optimal:
    for node_alias, node_info in node_info_topology_at_t.items():

        # Extract information
        node_assigned_alias     = node_info[0]
        node_lon, node_lat      = node_info[1:]
        
        # Plot satellite node as a regular scatter point with label
        if not any(gs_type in node_alias for gs_type in gs_alias_list):
            x, y = m(node_lon, node_lat)
            #plt.scatter(x, y, s=20, marker="o", facecolors='none', edgecolors='black', zorder=20)
            #plt.text(x, y-0.5, node_assigned_alias, fontsize=7, zorder=100)
    
# PLOT ALL GROUND STATIONS
optimal_route_at_t          = optimal_routes[time_index]
gs0                         = node_info_topology_at_t[optimal_route_at_t[0]]
gs1                         = node_info_topology_at_t[optimal_route_at_t[-1]]
optimal_endpoints           = [gs0, gs1]
x1, y1 = m(gs0[1], gs0[2])
x2, y2 = m(gs1[1], gs1[2])
#plt.scatter(x1, y1, s=100, marker='^', linewidth=1.5, edgecolors='k', facecolors='none', zorder=3, label="Source ("+gs0[0]+")")
#plt.scatter(x2, y2, s=100, marker='s', linewidth=1.5, edgecolors='k', facecolors='none', zorder=3, label="Destination ("+gs1[0]+")")
# plt.text(x1, y1-0.5, gs0[0], fontsize=7, zorder=100)
# plt.text(x2, y2-0.5, gs1[0], fontsize=7, zorder=100)
for node_alias, node_info in node_info_topology_at_t.items():
    if any(gs_type in node_alias for gs_type in gs_alias_list) and node_alias not in optimal_endpoints: # Rest of ground stations
        node_assigned_alias = node_info[0]
        if node_assigned_alias == 'London' or node_assigned_alias == 'Union-WVa' or node_assigned_alias == 'Sao Paulo':
            node_lon, node_lat = node_info[1:]
            x, y = m(node_lon, node_lat)
            plt.scatter(x, y, s=80, marker='s', facecolors='None', edgecolors='red', zorder=4, linewidth=2)
            if node_assigned_alias == 'London':
                x_text = x-2000000
            else:
                x_text = x
            plt.text(x_text, y-700000, node_assigned_alias, fontsize=10, color='red', zorder=75)
handles, labels = plt.gca().get_legend_handles_labels()
sat_marker = mlines.Line2D([], [], c='black', markerfacecolor='none', markersize=6, label='Satellite', marker='o', linestyle='None')
gs_marker = mlines.Line2D([], [], c='purple', markerfacecolor='none', markersize=6, label='Ground Station (GW, CT, IE)', marker='p', linestyle='None')

# PLOT OPTIMAL ROUTE
optimal_lon = np.array([0., ] * len(optimal_route_at_t))
optimal_lat = np.array([0., ] * len(optimal_route_at_t))
for indx, optimal_node in enumerate(optimal_route_at_t):
    optimal_node_info   = node_info_topology_at_t[optimal_node]
    optimal_lon[indx]   = optimal_node_info[1]
    optimal_lat[indx]   = optimal_node_info[2]
    x, y = m(optimal_lon[indx], optimal_lat[indx])
    #plt.text(x, y+0.3, optimal_node_info[0], fontsize=7, zorder=4)

# Define colors for different types of connections
colors = {'sat-sat': 'blue', 'sat-gs': 'green', 'gs-gs': 'red'}
blue_line = mlines.Line2D([], [], color=colors['sat-sat'], markersize=5, label='Sat-Sat', linestyle='--')
green_line = mlines.Line2D([], [], color=colors['sat-gs'], markersize=5, label='GS-Sat', linestyle='--')
red_line = mlines.Line2D([], [], color=colors['gs-gs'], markersize=5, label='GS-GS', linestyle='--')
handles.extend([sat_marker, gs_marker, blue_line, green_line, red_line])
labels.extend([sat_marker.get_label(), gs_marker.get_label(), blue_line.get_label(), green_line.get_label(), red_line.get_label()])

# Iterate over pairs of nodes in the optimal route
for i in range(len(optimal_route_at_t) - 1):
    
    # Assign node for comparison
    node1 = optimal_route_at_t[i]
    node2 = optimal_route_at_t[i+1]

    # Check if nodes are terrestrial nodes
    node1_gs = any(gs_type in node1 for gs_type in gs_alias_list)
    node2_gs = any(gs_type in node2 for gs_type in gs_alias_list)

    # Determine the type of connection
    if node1_gs and node2_gs:
        color = colors['gs-gs']
    elif node1_gs or node2_gs:
        color = colors['sat-gs']
    else:
        color = colors['sat-sat']

    # Plot the line with the chosen color
    x, y = m([optimal_lon[i], optimal_lon[i+1]], [optimal_lat[i], optimal_lat[i+1]])
    #plt.plot(x, y, '--', linewidth=1.5, c=color, zorder=1)

# PLOT INFORMATION
#plt.title('FW Algorithm: '+str(total_num_sat)+' nodes (time: '+str(dt_hist[time_index])+') (hops='+str(len(optimal_route_at_t)-1)+')')
print('FW Algorithm: '+str(total_num_sat)+' nodes (time: '+str(dt_hist[time_index])+') (hops='+str(len(optimal_route_at_t)-1)+')')
# plt.xlabel('Longitude')
# plt.ylabel('Latitude')
#plt.legend(fancybox=True, framealpha=1, handles=handles, labels=labels, loc='upper left').set_zorder(100)
plt.tight_layout()
plt.show()
# plt.savefig('/home/barbourbruce/dynamic-topology-generator/test_plot.pdf')

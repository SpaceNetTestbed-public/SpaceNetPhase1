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
plot_GSs                = False
plot_only_optimal       = False
plot_in_3D              = False
plot_debug              = False   #Plots linked sats to the given sat index and their respective orbits
plot_optimal_orbits     = False   #Plots all the orbits involved in the optimal path
lon0_3d                 = 0  #-35
lat0_3d                 = -40
timestamp               = "2024_07_21_09_00_0"
gs_filepath             = open('/home/spacenet/simulator/dynamic-topology-generator/output/terrestrial_info/terrestrial_1721552400.txt', 'r')
tle_file                = open('/home/spacenet/Desktop/jacktles/alt500km/inc50/starlink_tles/starlink_1721552400', 'r')
optimal_route_filepath  = '/home/spacenet/simulator/dynamic-topology-generator/output/optimal_routes/starlink/best_path_'+timestamp+'.0.txt'
conn_filepath           = '/home/spacenet/simulator/dynamic-topology-generator/output/connectivity_matrix/starlink/topology_'+timestamp+'.0.txt'
node_indices_filepath   = '/home/spacenet/simulator/dynamic-topology-generator/output/node_indices/starlink/nodeindex_1721552400.txt'
orb_sat_txt             = '/home/spacenet/simulator/dynamic-topology-generator/output/satellites_orbits/orbits_satellites.txt'
optimal_weight = 24349
number_of_orbits = 72  #$
num_plotorbit = range(1)  #$ Number of orbits to plot
gs0 = (-74.003663, 40.717042) # NYC
gs1 = (103.850070, 1.289670) # Singapore
ref = 558  # index of satellite to be debugged for ISLs 

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
plotted_sat_index                   = {}
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
    sats_from_tle_dict[tle_sat_obj.name] = tle_sat_obj  # Example: sats_from_tle_dict['STARLINK-1234'] = EarthSatellte type

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
        node_index              = int(node_index)  # example: 100
        node_alias              = node_alias[:-1]  # example: STARLINK-1099
        node_type               = node_alias.split('-')
        node_type               = node_type[0]
        node_alias_to_index_topology_dict[node_alias] = node_index  # example: STARLINK-1111 --> 1112
        node_index_to_alias_topology_dict[node_index] = node_alias  # example: 100  --> STARLINK-1099
        if node_type in gs_alias_list:
            total_num_gs += 1
        else:
            total_num_sat += 1
if total_num_gs > total_num_gs_in_gsfile:
    raise ValueError("The number of ground stations observed in the topology is greater than in the provided ground station file!")

#$ Step 1.5: Get sat index sorted for each orbit (Dictionary which gives all sat IDs (NOT sat names (STARLINK-####)) sorted for each orbit)
sat_orbit_index = [[] for i in range(number_of_orbits)]  #index:orbit_number | value:list of sats in that orbit
sat_index_orbit = [[] for i in range(total_num_sat)]  #index:sat index | value:orbit number
with open(orb_sat_txt, 'r') as orbsat_file:
    for i, sat_index in enumerate(orbsat_file):
        if i < total_num_sat:
            line = sat_index.split("\n")
            IDs = line[0].split()   # IDs[0] --> orbit index   IDs[1] --> sat index
            sat_orbit_index[int(IDs[0])-1].append(int(IDs[1])) # Actually a nested list  ||||| Input:orbit_index  Output:sat_index list
            sat_index_orbit[int(IDs[1])].append(int(IDs[0])-1) # Actually a nested list  ||||| Input:sat_index  Output:orbit_index

# ================================================================================================
# FILE PARSING - CONNECTIVITY FILE (Checking only sat links)
# ================================================================================================
conn_mat = {}
with open(conn_filepath, 'r') as conn_file:
    for i, conn_index in enumerate(conn_file):
        line = conn_index.split(",")
        if int(line[0]) < total_num_sat and int(line[1]) < total_num_sat:       ###### Both endpoints should be satellites (exclude all the GSL links) 
            if line[0] not in conn_mat.keys():
                conn_mat[line[0]] = [int(line[1])]
            else:
                conn_mat[line[0]].append(int(line[1])) # Actually a nested list
num_links = []
count = 0
for i, js in conn_mat.items():
    num_links.append(len(js))
    if len(js)>4:
        count += 1

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
        
        #$ Adding only sat node optimal route data
        optimal_route_satonly_at_epoch = []   # Indices of all the sats in optimal route
        for i, nodes in enumerate(optimal_route_at_epoch):
            a_node = nodes.split("-")
            if a_node[0] == 'STARLINK':
                optimal_route_satonly_at_epoch.append(route_node_indices[i])

        # Append to complete list
        optimal_routes.append(optimal_route_at_epoch)

        # (Debugging purpose) List of orbits used for optimal path 
        optimal_orbits = [sat_index_orbit[int(k)][0] for k in optimal_route_satonly_at_epoch]
        optimal_orbits = np.unique(optimal_orbits)
# ================================================================================================
# DEFINE CURRENT TIME INDEX
# ================================================================================================
t_current   = ts.from_datetime(dt_hist[time_index])

#$ Filter sat's (lat,long) orbit-wise as given in nodeindices (rearranged from arranged sats in nodeindex)
lats, lons = [], []
for k in range(max(optimal_orbits)+1):
    for sat_in_orb in sat_orbit_index[k]:
        aliass = node_index_to_alias_topology_dict[sat_in_orb]
        sat_at_t    = sats_from_tle_dict[aliass].at(t_current)
        lat, lon    = sat_at_t.subpoint().latitude.degrees, sat_at_t.subpoint().longitude.degrees
        plotted_sat_index[sat_in_orb] = (lat, lon)
        lats.append(lat)
        lons.append(lon)

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
    #m = Basemap(projection='cyl', llcrnrlat=0, urcrnrlat=80, llcrnrlon=-140, urcrnrlon=-40, resolution='c')
else:
    m0 = Basemap(projection='ortho', lat_0=lat0_3d, lon_0=lon0_3d, resolution=None)
    #m = Basemap(projection='ortho', lat_0=lat0_3d, lon_0=lon0_3d, llcrnrx=-m0.urcrnrx/1.75, llcrnry=0, urcrnrx=m0.urcrnrx/1.75, urcrnry=m0.urcrnry/1.75, resolution='c')
    m = Basemap(projection='ortho', lat_0=lat0_3d, lon_0=lon0_3d, llcrnrx=-m0.urcrnrx/1.75, llcrnry=-m0.urcrnry/10.75, urcrnrx=m0.urcrnrx/1.75, urcrnry=m0.urcrnry/1.75, resolution='c')
m.drawcoastlines()
m.drawcountries()
m.fillcontinents(color='lightgray', lake_color='white')
m.drawmapboundary(fill_color='white')

# PLOT ALL SATELLITE NODES IN TOPOLOGY
if not plot_only_optimal:

    #$ Custom plotting individual orbits
    if plot_optimal_orbits:
        for i in optimal_orbits:
            X = []
            Y = []
            col = np.random.rand(3,)
            #col = [1, 0, 0]
            for j in sat_orbit_index[i]:
                x, y = m(lons[j], lats[j])
                X.append(x)
                Y.append(y)
                #plt.scatter(x, y, s=13, marker="o", c=col, edgecolors=col, facecolors='none', zorder=5)
            X.append(X[0])  #completing the orbit
            Y.append(Y[0])  #completing the orbit
            plt.plot(X, Y, color=col, marker='o')

    for node_alias, node_info in node_info_topology_at_t.items():

        # Extract information
        node_assigned_alias     = node_info[0]
        node_lon, node_lat      = node_info[1:]
        
        # Plot satellite node as a regular scatter point with label
        if not any(gs_type in node_alias for gs_type in gs_alias_list):
            x, y = m(node_lon, node_lat)
            if node_assigned_alias == node_index_to_alias_topology_dict[ref] and plot_debug:
                plt.scatter(x, y, s=50, marker="o", facecolors='none', edgecolors='red', zorder=20)
                plt.text(x, y-0.5, ref, fontsize=15, color='red', zorder=100)
            elif num_links[node_alias_to_index_topology_dict[node_assigned_alias]]>4:
                plt.scatter(x, y, s=20, marker="o", facecolors='darkcyan', edgecolors='darkcyan', zorder=20)
                plt.text(x, y-0.5, node_alias_to_index_topology_dict[node_assigned_alias], fontsize=7, zorder=100)
            else:
                plt.scatter(x, y, s=20, marker="o", facecolors='none', edgecolors='black', zorder=20)
            #plt.text(x, y-0.5, node_assigned_alias, fontsize=7, zorder=100)

####################  #$ PLOT SATS AND ITS LINKS WITH HIGHLIGHTED ORBITS (DEBUGGING ZONE STARTS) ###########################
if plot_debug:
    print(ref, node_index_to_alias_topology_dict[ref])

    ###### Get unique list of relevant orbits
    debug_orbits = [sat_index_orbit[int(k)][0] for k in conn_mat[str(ref)]]
    debug_orbits = np.unique(debug_orbits)

    ###### Lat-Lon of all satellites in the relevant orbits
    lats, lons = [], []
    for k in range(max(debug_orbits)+1):
        for sat_in_orb in sat_orbit_index[k]:
            aliass = node_index_to_alias_topology_dict[sat_in_orb]
            sat_at_t    = sats_from_tle_dict[aliass].at(t_current)
            lat, lon    = sat_at_t.subpoint().latitude.degrees, sat_at_t.subpoint().longitude.degrees
            plotted_sat_index[sat_in_orb] = (lat, lon)
            lats.append(lat)  # lats till the max index in the debug_orbit list
            lons.append(lon)  # lons till the max index in the debug_orbit list

    for i in debug_orbits:
        X, Y = [], []
        for j in sat_orbit_index[i]:
            x, y = m(lons[j], lats[j])
            X.append(x)
            Y.append(y)
        X.append(X[0])  #completing the orbit
        Y.append(Y[0])  #completing the orbit
        if sat_index_orbit[j] == sat_index_orbit[ref]:
            plt.plot(X, Y, color='green', marker=',')   ##### Plots current orbit
        else:
            plt.plot(X, Y, color='blue', marker=',')   ##### Plots neighbouring orbits 

    for neighbours in conn_mat[str(ref)]:
        print(neighbours, node_index_to_alias_topology_dict[neighbours])
        info = node_info_topology_at_t[node_index_to_alias_topology_dict[neighbours]]
        neighbour_lon, neighbour_lat = info[1:]
        x, y = m(neighbour_lon, neighbour_lat)
        if sat_index_orbit[neighbours] == sat_index_orbit[ref]:
            plt.scatter(x, y, s=30, marker="o", facecolors='green', edgecolors='green', zorder=20)
            plt.text(x, y-0.5, neighbours, fontsize=15, zorder=100)
        else:
            plt.scatter(x, y, s=30, marker="o", facecolors='blue', edgecolors='blue', zorder=20)
            plt.text(x, y-0.5, neighbours, fontsize=15, zorder=100)

######################################################### DEBUGGING ZONE ENDS ###################################

# PLOT ALL GROUND STATIONS
optimal_route_at_t          = optimal_routes[time_index]
gs0                         = node_info_topology_at_t[optimal_route_at_t[0]]
gs1                         = node_info_topology_at_t[optimal_route_at_t[-1]]
optimal_endpoints           = [gs0, gs1]
x1, y1 = m(gs0[1], gs0[2])
x2, y2 = m(gs1[1], gs1[2])
plt.scatter(x1, y1, s=150, marker='^', linewidth=1.5, edgecolors='r', facecolors='none', zorder=3, label="Source ("+gs0[0]+")")
plt.scatter(x2, y2, s=150, marker='s', linewidth=1.5, edgecolors='r', facecolors='none', zorder=3, label="Destination ("+gs1[0]+")")
# plt.text(x1, y1-0.5, gs0[0], fontsize=7, zorder=100)
# plt.text(x2, y2-0.5, gs1[0], fontsize=7, zorder=100)
if plot_GSs:
    for node_alias, node_info in node_info_topology_at_t.items():
        if any(gs_type in node_alias for gs_type in gs_alias_list) and node_alias not in optimal_endpoints: # Rest of ground stations
            node_assigned_alias = node_info[0]
            node_lon, node_lat = node_info[1:]
            x, y = m(node_lon, node_lat)
            plt.scatter(x, y, s=20, marker='o', facecolors='None', edgecolors='purple', zorder=4, linewidth=2)
            #plt.text(x, y-700000, node_assigned_alias, fontsize=10, color='red', zorder=75)
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
    plt.scatter(x, y, s=20, marker='X', facecolors='k', edgecolors='k', zorder=4, linewidth=2)
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
    plt.plot(x, y, '--', linewidth=4.5, c=color, zorder=1)

# PLOT INFORMATION
#plt.title('FW Algorithm: '+str(total_num_sat)+' nodes (time: '+str(dt_hist[time_index])+') (hops='+str(len(optimal_route_at_t)-1)+')')
print('FW Algorithm: '+str(total_num_sat)+' nodes (time: '+str(dt_hist[time_index])+') (hops='+str(len(optimal_route_at_t)-1)+')')
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.title('Timestep: ' + timestamp + ' |  # of Hops: ' + str(len(optimal_route_at_epoch)-1) + ' | # of involved orbits: ' + str(len(optimal_orbits)) + ' | Non "+ grid" sats: ' + str(count))
#plt.legend(fancybox=True, framealpha=1, handles=handles, labels=labels, loc='upper left').set_zorder(100)
plt.tight_layout()
plt.show()
# plt.savefig('/home/barbourbruce/dynamic-topology-generator/test_plot.pdf')

import matplotlib.pyplot as plt
from skyfield.api import load, EarthSatellite
import re
from datetime import datetime, timezone
from mpl_toolkits.basemap import Basemap



# Control Var:
time_index = 0
tle_file = open('../utils/starlink_tles/starlink_1705701600', 'r')
optimal_route_filepath = '../output/optimal_routes/starlink/best_path_2024_01_19.txt'
node_indices_filepath = '../output/node_indices/starlink/nodeindex_1705701600.txt'
total_sat_num = 1440
gs0 = (-74.003663, 40.717042, "NYC") # NYC
gs1 = (-0.1278, 51.5074, "London")  # London

# Step 0: Initialize objects
ts = load.timescale()
sats = []
plotted_sat_alias = {}
plotted_alias_index = {}

# Step 1: Read TLE from the file
lines = tle_file.readlines()
for i in range(0, len(lines), 3):
    tle_name            = lines[i]
    tle_first_line      = lines[i+1]
    tle_second_line     = lines[i+2]
    current_sat         = EarthSatellite(tle_first_line, tle_second_line)
    current_sat.name    = re.match(r"STARLINK-\d+", tle_name).group(0)
    sats.append(current_sat)

# Step 2: Create node index and ID catalog
node_alias_catalog = {}
node_index_catalog = {}
gs_index_catalog = {}
with open(node_indices_filepath, 'r') as node_indices_file:
    for i, node_index in enumerate(node_indices_file):
        if i < total_sat_num:
            index, nodeID = node_index.split(":")
            node_alias_catalog[nodeID[:-1]] = int(index)
            node_index_catalog[int(index)] = nodeID[:-1]
        elif i >= total_sat_num:
            gs_index, gsID = node_index.split(":")
            gs_index_catalog[gsID[:-1]] = int(gs_index)

# Step 3: Create optimal path list and assign local index
optimal_routes = []
dt_hist = []
with open(optimal_route_filepath, 'r') as optimal_file:
    for route_line in optimal_file:
        dt, numbers_part    = route_line.split(":", 1)

        # Time history
        dt_list             = dt.strip('()').split("_")
        y, mon, d, hr, m    = map(int, dt_list[:5])
        sec                 = float(dt_list[5])
        dt_hist.append(datetime(y, mon, d, hr, m, int(sec), int((sec - int(sec)) * 1000000), tzinfo=timezone.utc))

        # Route history
        route               = numbers_part.split(", ")[1:-1]
        route_alias         = []

        for node in route:
            route_alias.append(node_index_catalog[int(node)])
        optimal_routes.append(route_alias)

# Define time range
t = ts.from_datetime(dt_hist[time_index])

# Define longitude and latitude bounds (example bounds)
lon_min, lon_max = min([gs0[0], gs1[0]])-5, max([gs0[0], gs1[0]])+5
lat_min, lat_max = min([gs0[1], gs1[1]])-5, max([gs0[1], gs1[1]])+5

# Step 4: Compute and filter positions
lats, lons = [], []
for i, sat in enumerate(sats):

    # Check if sat in node alias
    if sat.name in node_alias_catalog:
        sat_at_t    = sat.at(t)
        lat, lon    = sat_at_t.subpoint().latitude.degrees, sat.at(t).subpoint().longitude.degrees
        plotted_sat_alias[sat.name] = (lat, lon)
        plotted_alias_index[sat.name] = node_alias_catalog[sat.name]
        lats.append(lat)
        lons.append(lon)
    # if lon_min <= lon <= lon_max and lat_min <= lat <= lat_max:
    #     lats.append(lat)
    #     lons.append(lon)

# Step 5: Plot the positions (continued)
plt.figure(figsize=(10, 5))
if lats and lons:  # Check if lists are not empty
    
    plt.scatter(lons, lats, s=1, marker="o", facecolors='none', edgecolors='black', zorder=20)
    for lon, lat, label in zip(lons, lats, plotted_alias_index.values()):
        plt.text(lon, lat - 0.5, label, ha='center')  # Adjust 0.5 as needed to position the text


    plt.scatter(gs0[0], gs0[1], s=100, marker='x', linewidth=2, c='g', zorder=20, label=gs0[2])
    plt.scatter(gs1[0], gs1[1], s=100, marker='x', linewidth=2, c='b', zorder=20, label=gs1[2])

    # Plot the routes
    route_points = [plotted_sat_alias[node_alias] for node_alias in optimal_routes[time_index]]
    route_lats, route_lons = zip(*route_points)
    route_lats = (gs0[1],) + route_lats + (gs1[1],)
    route_lons = (gs0[0],) + route_lons + (gs1[0],)

    plt.plot(route_lons, route_lats, '--*', linewidth=1.5, c='r', zorder=12)  # Plot with lines and markers

    m = Basemap(projection='cyl', llcrnrlat=-60, urcrnrlat=60, llcrnrlon=-180, urcrnrlon=180, resolution='c')
    m.drawcoastlines()
    m.drawcountries()
    m.fillcontinents(color='lightgray')

    plt.title('FW Algorithm: '+str(total_sat_num)+' nodes (time: '+str(dt_hist[time_index])+') (hops='+str(len(optimal_routes[time_index])+1)+')')
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    # plt.xlim((lon_min, lon_max))
    # plt.ylim((lat_min, lat_max))
    plt.legend(loc='upper left')
    plt.show()
    plt.savefig('../output/Figures/spacenet_node'+str(total_sat_num)+'_epoch'+str(time_index)+'.svg', dpi=300, bbox_inches='tight')
else:
    print("No data points to plot.")
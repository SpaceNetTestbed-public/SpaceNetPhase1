''''

SpaceNet: Generalized Plotting for Walker type constellation

AUTHOR:         Suryansh Aryan, 2025
                Virginia Tech


This Python script acts as a generalized plotting script that can incorporate CustomSatellite and Skyfield Earthsatellite type satellites in the 
constellation. 
Note: It can plot only one type at a time (If the simulation contains both Earthsatellites and Customsatellites type constellation satellites then it fails!)
'''

import os
import matplotlib
import gif_utils
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from astropy import units as u
from astropy.coordinates import spherical_to_cartesian, cartesian_to_spherical
from skyfield.api import load, EarthSatellite
from poliastro.twobody import Orbit
import poliastro.bodies as pbodies
from poliastro.util import Time
import plotly.graph_objects as go

import sys
pwd = os.getcwd()
sys.path.append("dynamic-topology-generator/")
from mobility.lunar_dyn_utils import *
os.chdir(pwd)

import yaml

import re
import numpy as np
from datetime import datetime, timezone
from mpl_toolkits.basemap import Basemap



# ================================================================================================
# >>> SCRIPT CONTROL - EDIT HERE <<<
# ================================================================================================
time_index              = 0
plot_GSs                = False
plot_only_optimal       = False
plot_in_3D              = False
plot_debug              = False   #Plots linked sats to the given sat index and their respective orbits AND also darkcyan satellites that have more than 4 ISLs
plot_optimal_orbits     = False   #Plots all the orbits involved in the optimal path
make_gif                = True   # Makes gif of all the timestep plots in this file existyin in current output path  (REQUIREMENTS: CONNECTIVITY FILES AND OPTIMAL_PATH FILES SHOULD BE EXISTING AND SEPERATE FILES FOR EACH TIMESTEP | line 153 hardcode should be rechecked)
lon0_3d                 = 50 #30  
lat0_3d                 = -15 #-35      
ll                      = [0.5, 0.5]   #scaling of the 3D plot (non-negetive) (lower left point) [0.7, 0.7]
ur                      = [0.5, 0.5]   #scaling of the 3D plot (non-negetive) (upper right point) [0.7, 0.7]
timestamp               = "2024_9_27_22_15_6"
timestamp2               = "2024_09_27_22_15_6"
tle_unix_timestamp      = "1727475306"
outputfolder_path       = "./users/johndoe/jd1/"
operator_name           = "starlink"
gif_name                = 'output_gif'
number_of_orbits        = 72+5
_timespan               = 120 #Important to change if doesnt match Phase1 settings
# gs_filepath             = open(outputfolder_path + 'output/terrestrial_info/terrestrial_'+timestamp+'.0.txt', 'r')
# tle_file                = open('/home/netsatlab/backend_test/dynamic-topology-generator/utils/'+operator_name+'_tles/'+operator_name+'_'+tle_unix_timestamp, 'r')
optimal_route_filepath  = outputfolder_path + 'output/optimal_routes/'+operator_name+'/best_path_'+timestamp2+'.0.txt'
conn_filepath           = outputfolder_path + 'output/connectivity/'+operator_name+'/topology_'+timestamp2+'.0.txt'
node_indices_filepath   = outputfolder_path + 'output/node_indices/'+operator_name+'/nodeindex_'+timestamp+'.0.txt'
orb_sat_txt             = outputfolder_path + 'output/satellites_orbits/orbits_satellites.txt'
################### Folder Paths ######################
conn_folder             = outputfolder_path + 'output/connectivity/'+operator_name+'/'
opt_route_folder        = outputfolder_path + 'output/optimal_routes/'+operator_name+'/'
gif_path                = outputfolder_path + 'GIF/'
#######################################################
Main_body = 'Moon'
Third_body = 'Earth'
ref = 1 #1263  # index of satellite to be debugged for ISLs (Only use with plot_debug=True)
shell_color = {1584:"orange",1814:"green"}
#shell_color = {440:"orange"}
default_projection_for_moon = "ortho"


ts = load.timescale()
sats_from_tle_dict                  = {}
gs_dict                             = {}
plot_sat_alias_dict                 = {}
plot_sat_alias_to_index_dict        = {}
node_alias_to_index_topology_dict   = {}
node_index_to_alias_topology_dict   = {}
node_info_topology_at_t             = {}
plotted_sat_index                   = {}
conn_mat                            = {}
optimal_routes                      = []
optimal_orbits                      = []
num_links                           = []
dt_hist                             = []
epoch_hist                          = []
lats                                = []
lons                                = []
gs_alias_list                       = ["CT", "GS", "GW", "IE"]
total_num_sat                       = 0
total_num_gs                        = 0
t_current                           = None

PLOTLY_SAT_MARKERSIZE = 5
PLOTLY_GS_MARKERSIZE = 10
PLOTLY_PATH_LINEWIDTH = 3
PLOTLY_ORBIT_LINEWIDTH = 1
PLOTLY_SAT_ALPHA = 1
PLOTLY_GS_ALPHA = 1
PLOTLY_PATH_ALPHA = 0.6
PLOTLY_ORBIT_ALPHA = 0.5



def parsing_files(filetype, sat_operator='STARLINK'):
    """
    Reading different files from Phase1
    Input:
            filetype (str) : Name of the file to be parsed (OPTIONS: tle, gs, node_index)
    """
    global sats_from_tle_dict, tle_file, gs_filepath, node_indices_filepath, total_num_gs, total_num_sat

    if filetype == 'tle':
        if sat_operator=='STARLINK':
            tle_lines = tle_file.readlines()
            for i in range(0, len(tle_lines), 3):
                tle_sat_name                    = tle_lines[i]
                tle_first_line                  = tle_lines[i+1]
                tle_second_line                 = tle_lines[i+2]
                tle_sat_obj                     = EarthSatellite(tle_first_line, tle_second_line)
                tle_sat_obj.name                = re.match(r"STARLINK-\d+", tle_sat_name, re.IGNORECASE).group(0)
                sats_from_tle_dict[tle_sat_obj.name] = tle_sat_obj  # Example: sats_from_tle_dict['STARLINK-1234'] = EarthSatellte type
        elif sat_operator=='vdvdfv':
            tle_lines = tle_file.readlines()
            tle_name = tle_file.name.split('/')
            tle_timestamp = int(tle_name[-1].split('_')[-1])
            julian_date = ( tle_timestamp / 86400.0 ) + 2440587.5
            epoch = Time(str(julian_date),format="jd",scale="tdb")
            for i in range(0, len(tle_lines), 3):
                tle_sat_name                    = tle_lines[i].strip('\n')
                tle_first_line                  = list(line for line in tle_lines[i+1].strip("\n").split(" ") if line)
                tle_second_line                 = list(line for line in tle_lines[i+2].strip("\n").split(" ") if line)
                mean_motion                     = float(tle_second_line[7]) * 2 * np.pi / 86400
                SMA                             = ((pbodies.Moon.k.to_value(u.km**3 / u.s**2) / (mean_motion ** 2)) ** (1. / 3.))*(u.km)  #kms   
                sat_data                        = Orbit.from_classical(pbodies.Moon, SMA, float("."+tle_second_line[4])*u.one, float(tle_second_line[2])*u.deg, float(tle_second_line[3])*u.deg, float(tle_second_line[5])*u.deg, float(tle_second_line[6])*u.deg, epoch) #km, kg, s, deg
                sat_obj                         = CustomSatellites(tle_sat_name, sat_data, Main_body, Third_body)
                sat_name                        = re.match(r"vdvdfv-\d+", tle_sat_name).group(0)
                sats_from_tle_dict[sat_name] = sat_obj  # Example: sats_from_tle_dict['STARLINK-1234'] = EarthSatellte type
        elif sat_operator=='LUNAR':
            os.chdir('../')
            sats_from_tle_dict = satellites_from_tle(tle_file.name, _timespan, Main_body, Third_body)
            os.chdir('library/')

    elif filetype == 'gs':

        gs_lines = gs_filepath.readlines()
        total_num_gs_in_gsfile = len(gs_lines)
        for i in range(0, len(gs_lines), 1):
            gs_info = gs_lines[i][:-1].split(":")
            gs_dict[gs_info[1]] = (int(gs_info[0]), gs_info[2], float(gs_info[4]), float(gs_info[3]))
    
    elif filetype == 'node_index':

        # (node_alias_to_index_topology_dict{} & node_alias_to_index_topology_dict{} contains satellite and GS indices and alias)
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



def index_orbit_relation():

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

    return sat_orbit_index, sat_index_orbit


def figuring_latlon(t_current, orbits, sat_orbit_index):
    if len(orbits) == 0:
        return [], []
    
    for k in range(max(orbits)+1):
        if k >= len(sat_orbit_index): continue 
        for sat_in_orb in sat_orbit_index[k]:
            aliass = node_index_to_alias_topology_dict[sat_in_orb]
            if aliass in sats_from_tle_dict:
                if operator_name.upper()=='STARLINK':
                    sat_at_t    = sats_from_tle_dict[aliass].at(t_current)
                    lat, lon    = sat_at_t.subpoint().latitude.degrees, sat_at_t.subpoint().longitude.degrees
                    plotted_sat_index[sat_in_orb] = (lat, lon)
                else:
                    # Custom logic (assumes imports exist)
                    sat_obj     = sats_from_tle_dict[aliass]
                    sat_ephem   = sat_obj.get_ephem()
                    r, v        = sat_ephem.rv(t_current)
                    rr_PA       = frame_conversions(t_current, r, 'ICRS', 'PA')
                    _, lat, lon = frame_conversions(t_current, rr_PA, 'PA', 'MER', format='spherical')
                    plotted_sat_index[sat_in_orb] = (lat, lon)
                
                if sat_in_orb in plotted_sat_index:
                    lats.append(plotted_sat_index[sat_in_orb][0])
                    lons.append(plotted_sat_index[sat_in_orb][1])
    return lats, lons


def debugging_section(sat_orbit_index, sat_index_orbit):
    # Parse Connectivity
    with open(conn_filepath, 'r') as conn_file:
        for i, conn_index in enumerate(conn_file):
            line = conn_index.split(",")
            if int(line[0]) < total_num_sat and int(line[1]) < total_num_sat:
                if line[0] not in conn_mat.keys():
                    conn_mat[line[0]] = [int(line[1])]
                else:
                    conn_mat[line[0]].append(int(line[1]))

    # Link Counts
    count = 0
    for i, js in conn_mat.items():
        num_links.append(len(js))
        if len(js)>4: count += 1

    # Parse Optimal Route
    with open(optimal_route_filepath, 'r') as optimal_path_file:
        for route_line in optimal_path_file:
            dt, route_info = route_line.split(": ", 1)
            dt_info = dt.strip('()').split("_")
            yr, mon, day, hr, min_ = map(int, dt_info[:5])
            sec = float(dt_info[5])
            dt_hist = datetime(yr, mon, day, hr, min_, int(sec), int((sec - int(sec)) * 1000000), tzinfo=timezone.utc)
            
            route_node_indices = route_info.split(", ")
            route_node_indices[-1] = route_node_indices[-1].strip()

            optimal_route_at_epoch = []
            optimal_route_satonly_at_epoch = []
            
            for idx, route_node_index in enumerate(route_node_indices):
                alias = node_index_to_alias_topology_dict[int(route_node_index)]
                optimal_route_at_epoch.append(alias)
                # Check if satellite
                if alias.split("-")[0] == operator_name.upper():
                    optimal_route_satonly_at_epoch.append(route_node_index)

            optimal_routes.append(optimal_route_at_epoch)
            
            # Get orbits involved
            for k in optimal_route_satonly_at_epoch:
                if int(k) < len(sat_index_orbit):
                    # Safety: Ensure sat index exists in map
                    if len(sat_index_orbit[int(k)]) > 0:
                        optimal_orbits.append(sat_index_orbit[int(k)][0])
            
    # Set global time
    global t_current
    t_current = ts.from_datetime(dt_hist) if operator_name.upper() == 'STARLINK' else Time(dt_hist, format='datetime', scale='tdb')

    # Calculate Lats/Lons
    unique_orbits = np.unique(optimal_orbits).tolist()
    global lats, lons
    lats, lons = figuring_latlon(t_current, unique_orbits, sat_orbit_index)
    
    # Create Topology Info
    global node_info_topology_at_t
    node_info_topology_at_t = node_topology_creation(t_current)

    return optimal_route_at_epoch, count, unique_orbits


def basemap_settings():

    if operator_name.lower() == 'starlink':
        if not plot_in_3D:
            m = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90, llcrnrlon=-180, urcrnrlon=180, resolution='c')
            #m = Basemap(projection='cyl', llcrnrlat=0, urcrnrlat=80, llcrnrlon=-140, urcrnrlon=-40, resolution='c')
        else:
            m0 = Basemap(projection='ortho', lat_0=lat0_3d, lon_0=lon0_3d, resolution=None)
            #m = Basemap(projection='ortho', lat_0=lat0_3d, lon_0=lon0_3d, llcrnrx=-m0.urcrnrx/1.75, llcrnry=-m0.urcrnry/10.75, urcrnrx=m0.urcrnrx/1.75, urcrnry=m0.urcrnry/2.75, resolution='c')   #US zoomed
            #m = Basemap(projection='ortho', lat_0=lat0_3d, lon_0=lon0_3d, llcrnrx=-m0.urcrnrx/1.75, llcrnry=-m0.urcrnry/1.75, urcrnrx=m0.urcrnrx/1.75, urcrnry=m0.urcrnry/1.75, resolution='c') #llcrnry=-m0.urcrnry/10.75  #Whole Globe
            m = Basemap(projection='ortho', lat_0=lat0_3d, lon_0=lon0_3d, llcrnrx=-ll[0]*m0.urcrnrx, llcrnry=-ll[1]*m0.urcrnrx, urcrnrx=ur[0]*m0.urcrnrx, urcrnry=ur[1]*m0.urcrnrx, resolution='c') #llcrnry=-m0.urcrnry/10.75  #Whole Globe
        m.drawcoastlines()
        m.drawcountries()
        m.fillcontinents(color='lightgray', lake_color='white')
        m.drawmapboundary(fill_color='white')
    
    else:
        from irfpy.moon import moon_map

        map = moon_map.MoonMapSmall()
        map_on_sphere = map.gridsphere()
        blon, blat = map_on_sphere.get_bgrid()
        level = map_on_sphere.get_average()

        if default_projection_for_moon == "hammer":
            #### Currently only supporting Moon (hammer projection)
            m = Basemap(projection='hammer', lon_0=0)
            m.drawparallels([-45, 0, 45])
            m.drawmeridians([-90, 0, 90])
            m.pcolormesh(blon, blat, level, latlon=True, cmap='gray')

        elif default_projection_for_moon == "ortho":
            #### 3D plot (supports Moon) (ortho projection)
            m = Basemap(projection='ortho', lat_0=lat0_3d, lon_0=lon0_3d, resolution=None)
            m.drawparallels([-45, 0, 45])
            m.drawmeridians([-90, 0, 90])
            m.pcolor(blon, blat, level, latlon=True, cmap='gray')


    return m


def node_topology_creation(tt):

    # (node_info_topology_at_t{} --> key: SAT NAME + GS NAME (GS, CT, IE, GW)  ||||  values: (NAME (SAT NAME + GS TRUE NAME (London, Tokyo etc..)), lat, long)   )
    for topology_node_alias, topology_node_index in node_alias_to_index_topology_dict.items():

        # Check that it's a satellite
        if not any(gs_type in topology_node_alias for gs_type in gs_alias_list):

            # Extract satellite object based on node alias
            sat_node            = sats_from_tle_dict[topology_node_alias]

            if operator_name.upper()=='STARLINK': #SKYFIELD TYPE
                # Propagate satellite object to t
                sat_node_at_t       = sat_node.at(tt)

                # Compute longitude/latitude of satellite object at t
                lon_at_t, lat_at_t  = sat_node_at_t.subpoint().longitude.degrees, sat_node_at_t.subpoint().latitude.degrees

            else: #CustomSatellites Type
                # Propagate satellite object to t
                sat_ephem           = sat_node.get_ephem()
                sat_r_icrs, sat_v_icrs        = sat_ephem.rv(tt)
                rr_PA = frame_conversions(tt, sat_r_icrs, 'ICRS', 'PA')
                sat_dist, sat_lat, sat_lon = frame_conversions(tt, rr_PA, 'PA', 'MER', format='spherical')
                lon_at_t, lat_at_t  = sat_lon, sat_lat

            # Add to dictionary
            node_info_topology_at_t[sat_node.name] = (topology_node_alias, lon_at_t, lat_at_t)

        # And if it's a ground station
        else:

            # Add to dictionary
            node_info_topology_at_t[topology_node_alias] = gs_dict[topology_node_alias][1:]

    return node_info_topology_at_t


def final_plotting(optimal_route_at_epoch, count, unique_orbits_list):
    print("Generating Interactive 3D Plot with Plotly...")
    
    fig = go.Figure()

    # Lists for plotting
    sat_lons, sat_lats, sat_text, sat_colors = [], [], [], []
    gs_lons, gs_lats, gs_text = [], [], []

    # --- 1. PARSE NODES ---
    sat_shell_labels = []
    for node_alias, node_info in node_info_topology_at_t.items():
        node_assigned_alias = node_info[0]
        node_lon, node_lat = node_info[1], node_info[2]

        # Ground Stations
        if any(gs_type in node_alias for gs_type in gs_alias_list):
            gs_lons.append(node_lon)
            gs_lats.append(node_lat)
            gs_text.append(f"{node_alias} ({node_assigned_alias})")

        # Satellites
        else:
            sat_lons.append(node_lon)
            sat_lats.append(node_lat)
            sat_text.append(node_assigned_alias)

            # Color Logic
            c = 'orange'
            shell_label = 'Satellites'
            node_idx = node_alias_to_index_topology_dict.get(node_assigned_alias, 0)

            if plot_debug and node_assigned_alias == node_index_to_alias_topology_dict.get(ref, ""):
                c = 'red'
                shell_label = 'Debug: Reference Satellite'
            elif plot_debug and node_idx < len(num_links) and num_links[node_idx] > 4:
                c = 'darkcyan'
                shell_label = 'Debug: High Link Count'
            else:
                try:
                    shell_list = sorted(list(shell_color.keys()))
                    c = shell_color[shell_list[-1]]
                    shell_label = f'Shell {len(shell_list)} Satellites'
                    for shell_idx, s in enumerate(shell_list):
                        if node_idx < s:
                            c = shell_color[s]
                            shell_label = f'Shell {shell_idx + 1} Satellites'
                            break
                except:
                    pass
            sat_colors.append(c)
            sat_shell_labels.append(shell_label)

    # z-ordering: add_trace line types first and then marker types such that nodes
    # layer on top of the lines for better representation

    # --- 2. PLOT OPTIMAL ROUTE ---
    type_lons = {'sat-sat': [], 'sat-gs': [], 'gs-gs': []}
    type_lats = {'sat-sat': [], 'sat-gs': [], 'gs-gs': []}
    
    path_lons_for_center = []
    path_lats_for_center = []

    for i in range(len(optimal_route_at_epoch) - 1):
        node1 = optimal_route_at_epoch[i]
        node2 = optimal_route_at_epoch[i+1]
        info1 = node_info_topology_at_t.get(node1)
        info2 = node_info_topology_at_t.get(node2)
        if not info1 or not info2: continue

        path_lons_for_center.append(info1[1])
        path_lats_for_center.append(info1[2])

        node1_gs = any(gs_type in node1 for gs_type in gs_alias_list)
        node2_gs = any(gs_type in node2 for gs_type in gs_alias_list)
        
        key = 'sat-sat'
        if node1_gs and node2_gs: key = 'gs-gs'
        elif node1_gs or node2_gs: key = 'sat-gs'

        type_lons[key].extend([info1[1], info2[1], None])
        type_lats[key].extend([info1[2], info2[2], None])

    if len(optimal_route_at_epoch) > 0:
        last_node = optimal_route_at_epoch[-1]
        if last_node in node_info_topology_at_t:
            path_lons_for_center.append(node_info_topology_at_t[last_node][1])
            path_lats_for_center.append(node_info_topology_at_t[last_node][2])

    colors = {'sat-sat': 'blue', 'sat-gs': 'green', 'gs-gs': 'red'}
    for key in ['sat-sat', 'sat-gs', 'gs-gs']:
        if type_lons[key]:
            fig.add_trace(go.Scattergeo(
                lon=type_lons[key], lat=type_lats[key], mode='lines',
                line=dict(width=PLOTLY_PATH_LINEWIDTH, color=colors[key]),
                name=f'Link: {key.upper()}', opacity=PLOTLY_PATH_ALPHA
            ))

    # --- 3. PLOT DEBUG ORBITS ---
    orbit_lons, orbit_lats = [], []
    if plot_optimal_orbits:
        for orb_idx in unique_orbits_list:
            orb_sats = sat_orbit_index[orb_idx]
            temp_lons, temp_lats = [], []
            for s in orb_sats:
                if s in plotted_sat_index:
                    temp_lons.append(plotted_sat_index[s][1])
                    temp_lats.append(plotted_sat_index[s][0])
            if temp_lons:
                temp_lons.append(temp_lons[0]) 
                temp_lats.append(temp_lats[0])
                orbit_lons.extend(temp_lons + [None])
                orbit_lats.extend(temp_lats + [None])
        
        fig.add_trace(go.Scattergeo(
            lon=orbit_lons, lat=orbit_lats, mode='lines',
            line=dict(width=PLOTLY_ORBIT_LINEWIDTH, color='cyan', dash='dot'),
            name='Active Orbits', opacity=PLOTLY_ORBIT_ALPHA
        ))

    # --- 4. ADD TRACES (NODES) ---
    # One trace per shell so each gets its own legend entry. Plotly draws
    # one legend swatch per trace, not per color, so a single combined
    # trace could only ever show one "Satellites" entry regardless of how
    # many shells were present.
    shell_buckets = {}
    for lon, lat, text, color, label in zip(sat_lons, sat_lats, sat_text, sat_colors, sat_shell_labels):
        bucket = shell_buckets.setdefault(label, {'lon': [], 'lat': [], 'text': [], 'color': color})
        bucket['lon'].append(lon)
        bucket['lat'].append(lat)
        bucket['text'].append(text)

    for label, bucket in shell_buckets.items():
        fig.add_trace(go.Scattergeo(
            lon=bucket['lon'], lat=bucket['lat'], text=bucket['text'], mode='markers',
            marker=dict(size=PLOTLY_SAT_MARKERSIZE, color=bucket['color'], opacity=PLOTLY_SAT_ALPHA, symbol='circle'),
            name=label
        ))

    fig.add_trace(go.Scattergeo(
        lon=gs_lons, lat=gs_lats, text=gs_text, mode='markers',
        marker=dict(size=PLOTLY_GS_MARKERSIZE, color='purple', opacity=PLOTLY_GS_ALPHA, symbol='diamond', line=dict(width=1, color='white')),
        name='Ground Stations'
    ))

    # --- 6. LAYOUT CONFIGURATION ---
    fig.update_layout(
        title=f'Interactive 3D View | Time: {timestamp}',
        autosize=True,
        height=900,
        geo=dict(
            projection_type='orthographic',
            
            # SCALE: Reduced to 0.9 to ensure no clipping edges
            projection_scale=0.9,
                        
            fitbounds=False, # Disable auto-zoom
            
            showland=True, showocean=True, showcountries=True,
            landcolor='rgb(30, 30, 30)', oceancolor='rgb(10, 10, 20)',
            bgcolor='black', 
            showframe=False
        ),
        margin=dict(l=20, r=20, t=40, b=20), # Small margins to keep it clean
        paper_bgcolor='black', font=dict(color='white')
    )
    
    output_filename = gif_path + "interactive_plot.html"
    print(f"Saving plot to {output_filename}...")
    fig.write_html(output_filename)
    print("Done! Open 'interactive_plot.html' in your web browser to view.")


def gif_creator():

    global node_info_topology_at_t

    conn_mat_global = {}   # Only exists if make_gif exists
    COUNT_global = []      # Only exists if make_gif exists
    num_links_global = {}  # Only exists if make_gif exists
    conn_sorted_path = sorted(os.listdir(conn_folder))

    for itr, conn_path_iter in enumerate(conn_sorted_path):
        conn_mat = {}
        with open(conn_folder+conn_path_iter, 'r') as conn_file:
            for conn_index in conn_file:
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
        conn_mat_global[itr] = conn_mat
        COUNT_global.append(count)
        num_links_global[itr] = num_links
    #print(len(num_links_global[3]))

    optroute_global = {}                                            # Only exists if make_gif exists
    optroute_sorted_path = sorted(os.listdir(opt_route_folder))     # Only exists if make_gif exists
    node_info_topology_at_t_global = {}                             # Only exists if make_gif exists
    TIMESTAMPS = []
    # ================================================================================================
    # FILE PARSING - OPTIMAL PATH ASSIGNMENT
    # ================================================================================================
    for itr, route_path_iter in enumerate(optroute_sorted_path):
        optimal_routes = []
        with open(opt_route_folder+route_path_iter, 'r') as optimal_path_file:
            for route_line in optimal_path_file:

                # Separate datetime and route information
                dt, route_info          = route_line.split(": ", 1)

                # Build time history
                dt_info                 = dt.strip('()').split("_")
                yr, mon, day, hr, min   = map(int, dt_info[:5])
                sec                     = float(dt_info[5])
                datetime_obj = datetime(yr, mon, day, hr, min, int(sec), int((sec - int(sec)) * 1000000), tzinfo=timezone.utc)
                dt_hist = datetime_obj
                datetime_gd = Time(datetime_obj,format='datetime',scale='tdb')
                datetime_gd.format = 'jd'
                epoch_hist = datetime_gd
                timestamp               = str(yr)+"_"+str(mon)+"_"+str(day)+"_"+str(hr)+"_"+str(min)+"_"+str(int(sec))
                # Routing information
                route_node_indices      = route_info.split(", ")
                route_node_indices[-1]  = route_node_indices[-1][:-1]   # Removes the next line command
                
                # Assign the corresponding satellite alias with the route node indices
                optimal_route_at_epoch  = []
                for route_node_index in route_node_indices:
                    optimal_route_at_epoch.append(node_index_to_alias_topology_dict[int(route_node_index)])
                
                #$ Adding only sat node optimal route data
                # optimal_route_satonly_at_epoch = []   # Indices of all the sats in optimal route
                # for i, nodes in enumerate(optimal_route_at_epoch):
                #     a_node = nodes.split("-")
                #     if a_node[0] == 'STARLINK':
                #         optimal_route_satonly_at_epoch.append(route_node_indices[i])

                # Append to complete list
                optimal_routes.append(optimal_route_at_epoch)
        TIMESTAMPS.append(timestamp)
        optroute_global[itr] = optimal_routes

        # ================================================================================================
        # DEFINE CURRENT TIME INDEX
        # ================================================================================================
        if operator_name.upper() == 'STARLINK':
            t_current   = ts.from_datetime(dt_hist)
        elif operator_name.upper() == 'LUNAR':
            t_current   = epoch_hist

        node_info_topology_at_t = {}
        # ================================================================================================
        # CREATE DICTIONARY WITH SATELLITE GEODETIC POSITION
        # ================================================================================================
        node_info_topology_at_t = node_topology_creation(t_current)
        node_info_topology_at_t_global[itr] = node_info_topology_at_t
    # ================================================================================================
    # PLOTTING
    # ================================================================================================
    figures = []
    for itr in range(len(conn_sorted_path)):
        conn_mat = conn_mat_global[itr]
        num_links = num_links_global[itr]
        node_info_topology_at_t = node_info_topology_at_t_global[itr]
        optimal_routes = optroute_global[itr]
        count = COUNT_global[itr]
        timestamp = TIMESTAMPS[itr]

        fig = plt.figure(figsize=(16,9))
        font = {'family' : 'monospace', 
                'size' : 12}
        plt.rc('font', **font)

        # PLOT BASEMAP
        m = basemap_settings()

        # PLOT ALL SATELLITE NODES IN TOPOLOGY (IF not make_gif then just plot OTHERWISE make the GIF OFC!)
        if not plot_only_optimal:

            for node_alias, node_info in node_info_topology_at_t.items():

                # Extract information
                node_assigned_alias     = node_info[0]
                node_lon, node_lat      = node_info[1:]

                # Plot ONLY satellite node as a regular scatter point with label
                if not any(gs_type in node_alias for gs_type in gs_alias_list):

                    node_idx = node_alias_to_index_topology_dict[node_assigned_alias]
                    shell_span = 0
                    shell_list = list(shell_color.keys())
                    while node_idx>=shell_list[shell_span]:
                        shell_span += 1
                    coloring_shell_sats = shell_color[shell_list[shell_span]]

                    x, y = m(node_lon, node_lat)
                    if node_assigned_alias == node_index_to_alias_topology_dict[ref] and plot_debug:
                        plt.scatter(x, y, s=50, marker="o", facecolors='none', edgecolors='red', zorder=20)
                        #plt.text(x, y-0.5, ref, fontsize=15, color='red', zorder=100)
                    elif num_links[node_alias_to_index_topology_dict[node_assigned_alias]]>4 and plot_debug:
                        plt.scatter(x, y, s=20, marker="o", facecolors='darkcyan', edgecolors='darkcyan', zorder=20)
                        #plt.text(x, y-0.5, node_alias_to_index_topology_dict[node_assigned_alias], fontsize=7, zorder=100)
                    else:
                        plt.scatter(x, y, s=20, marker="o", facecolors=coloring_shell_sats, edgecolors=coloring_shell_sats, zorder=20)

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
        plt.scatter(x1, y1, s=150, marker='^', linewidth=1.5, edgecolors='yellow', facecolors='none', zorder=3, label="Source ("+gs0[0]+")")
        plt.scatter(x2, y2, s=150, marker='s', linewidth=1.5, edgecolors='yellow', facecolors='none', zorder=3, label="Destination ("+gs1[0]+")")
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
        shell_list = sorted(list(shell_color.keys()))
        sat_shell_markers = [
            mlines.Line2D([], [], c=shell_color[shell_list[i]], markerfacecolor='none', markersize=6,
                          label=f'Shell {i + 1} Satellites', marker='o', linestyle='None')
            for i in range(len(shell_list))
        ]
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
        handles.extend(sat_shell_markers + [gs_marker, blue_line, green_line, red_line])
        labels.extend([m.get_label() for m in sat_shell_markers] + [gs_marker.get_label(), blue_line.get_label(), green_line.get_label(), red_line.get_label()])
        plt.legend(handles=handles, labels=labels, loc='upper left',
                   bbox_to_anchor=(1.01, 1.0), fontsize=7, markerscale=0.8,
                   frameon=True, borderaxespad=0.0)

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
        print('FW Algorithm: '+str(total_num_sat)+' nodes (time: '+str(dt_hist)+') (hops='+str(len(optimal_route_at_t)-1)+')')
        # plt.xlabel('Longitude')
        # plt.ylabel('Latitude')
        #plt.title('Timestep: ' + timestamp + ' |  # of Hops: ' + str(len(optimal_routes[0])-1) + ' | Non "+ grid" sats: ' + str(count))
        plt.title('Timestep: ' + timestamp + ' |  # of Hops: ' + str(len(optimal_routes[0])-1))
        plt.tight_layout()
        #plt.show()

        figures.append(fig)
        print('Saved ' + str(itr+1))
        if itr<10:
            plt.savefig(gif_path+"fig0"+str(itr)+".jpg")
        else:
            plt.savefig(gif_path+"fig"+str(itr)+".jpg")
        plt.close('all')
    
    gif_utils.convert_gif(gif_path, gif_path+gif_name, 350)
    print('Exiting and saving...')

def extract_datetime(name: str):
    # Match both 2024_9_27_22_15_6.0 and 2024_09_27_22_15_06.0
    m = re.search(r'(\d{4})_(\d{1,2})_(\d{1,2})_(\d{1,2})_(\d{1,2})_(\d{1,2})', name)
    if m:
        y, mo, d, h, mi, s = map(int, m.groups())
        return datetime(y, mo, d, h, mi, s)
    return datetime.max  # fallback so invalid names sort last

def matchFilePath(folder, pattern):
    files = sorted(os.listdir(folder), key=extract_datetime)
    count = 0
    for filename in files:
        if re.match(pattern, filename) and (count == time_step_count or len(files) == 1):
            filepath = os.path.join(folder, filename)
            return filepath
        count += 1
    raise Exception("No matching pattern for folder: ", folder)


def file_to_unix(filename: str) -> int:
    # Match the datetime pattern (ignore .0 or extension)
    m = re.search(r'(\d{4}_\d{1,2}_\d{1,2}_\d{1,2}_\d{1,2}_\d{1,2})', filename)
    if not m:
        raise ValueError("No valid datetime found in filename")

    # Extract the groups
    year, month, day, hour, minute, second = map(int, m.group(1).split('_'))

    # Convert to datetime and then to timestamp
    dt = datetime(year, month, day, hour, minute, second, tzinfo=timezone.utc)
    return str(int(dt.timestamp()))

def midpoint(lat1, lon1, lat2, lon2):
    # convert to radians
    lat1, lon1 = math.radians(lat1), math.radians(lon1)
    lat2, lon2 = math.radians(lat2), math.radians(lon2)

    dlon = lon2 - lon1

    bx = math.cos(lat2) * math.cos(dlon)
    by = math.cos(lat2) * math.sin(dlon)

    lat3 = math.atan2(
        math.sin(lat1) + math.sin(lat2),
        math.sqrt((math.cos(lat1) + bx)**2 + by**2)
    )
    lon3 = lon1 + math.atan2(by, math.cos(lat1) + bx)

    return math.degrees(lat3), (math.degrees(lon3) + 540) % 360 - 180

def setConfigurations():
    if (len(sys.argv) == 2):
        gif_config_path = sys.argv[1]

    try:
        with open(gif_config_path, 'r') as stream:
            config = yaml.safe_load(stream)
    except Exception as exc:
        print(exc)
        config = None

    global plot_in_3D, outputfolder_path, gs_filepath, optimal_route_filepath, tle_file, conn_filepath, node_indices_filepath, orb_sat_txt, conn_folder, opt_route_folder, gif_path, _timespan, number_of_orbits
    global lon0_3d, lat0_3d, make_gif, shell_color, plot_GSs, plot_only_optimal, plot_debug, plot_optimal_orbits, time_step_count, timestamp

    plot_GSs                = config['plot_GSs']
    plot_only_optimal       = config['plot_only_optimal']
    plot_in_3D              = config['plot_in_3D']
    plot_debug              = config['plot_debug']   #Plots linked sats to the given sat index and their respective orbits AND also darkcyan satellites that have more than 4 ISLs
    plot_optimal_orbits     = config['plot_optimal_orbits']   #Plots all the orbits involved in the optimal path
    make_gif                = config['make_gif']   # Makes gif of all the timestep plots in this file existyin in current output path  (REQUIREMENTS: CONNECTIVITY FILES AND OPTIMAL_PATH FILES SHOULD BE EXISTING AND SEPERATE FILES FOR EACH TIMESTEP | line 153 hardcode should be rechecked)   
    ll                      = [0.5, 0.5]   #scaling of the 3D plot (non-negetive) (lower left point) [0.7, 0.7]
    ur                      = [0.5, 0.5]   #scaling of the 3D plot (non-negetive) (upper right point) [0.7, 0.7]
    outputfolder_path       = config['outputfolder_path']
    time_step               = config['time_step']
    try:
        with open(outputfolder_path + "sat_config.yaml", 'r') as stream:
            sat_config = yaml.safe_load(stream)
        with open(outputfolder_path + "main_config.yaml", 'r') as stream:
            main_config = yaml.safe_load(stream)
    except Exception as exc:
        print(exc)
        sat_config = None
        main_config = None
    operator_name           = sat_config['operator_name']
    gif_name                = 'output_gif'
    total_orbits = 0
    total_sats = 0
    for shell_name, shell in sat_config["shells"].items():
        orbit_total = shell["orbits"]
        total_orbits += orbit_total
        total_sats += shell["sat_per_orbit"] * shell["orbits"]
    number_of_orbits        = total_orbits
    groundStation1 = main_config['SourceNode']
    groundStation2 = main_config['DestNode']

    sx, sy = 0, 0
    dx, dy = 0, 0

    with open(main_config['GroundStationFile'], "r") as f:
        for line in f:
            parts = line.strip().split(",")
            value = parts[0]
            if groundStation1 == int(value):
                sy = float(parts[2])
                sx = float(parts[3])
            if groundStation2 == int(value):
                dy = float(parts[2])
                dx = float(parts[3])
    # center gif between the ground stations
    if (config['center_gif']):
        lat0_3d, lon0_3d = midpoint(sy, sx, dy, dx)
    else: 
        lon0_3d                 = config['long'] #30  
        lat0_3d                 = config['lat'] #-35   

    _timespan               = sat_config['Sim_Length']['TimeStepDuration'] * sat_config['Sim_Length']['TimeStepCount']  #Important to change if doesnt match Phase1 settings
    time_step_count         = 0
    if (time_step != 0 and _timespan % sat_config['Sim_Length']['TimeStepDuration'] == 0):
        time_step_count = time_step // sat_config['Sim_Length']['TimeStepDuration']
    optimal_route_filepath  = matchFilePath(outputfolder_path + 'output/optimal_routes/'+operator_name, r'^best_path.*.0.txt$')
    gs_filepath             = open(matchFilePath(outputfolder_path + 'output/terrestrial_info/', r'^terrestrial_.*.0.txt$'), 'r')
    node_indices_filepath   = matchFilePath(outputfolder_path + 'output/node_indices/'+operator_name, r'^nodeindex.*.0.txt')
    tle_file_name           = os.listdir(outputfolder_path + 'output/TLE/')[0]
    tle_file                = open(outputfolder_path + 'output/TLE/' + tle_file_name, 'r')
    conn_filepath           = matchFilePath(outputfolder_path + 'output/connectivity/'+operator_name, r'^topology_.*.0.txt')
    orb_sat_txt             = outputfolder_path + 'output/satellites_orbits/orbits_satellites.txt'
    conn_folder             = outputfolder_path + 'output/connectivity/'+operator_name+'/'
    opt_route_folder        = outputfolder_path + 'output/optimal_routes/'+operator_name+'/'
    gif_path                = outputfolder_path + 'gifs/' + config['gif_name'] + '/'
    os.makedirs(gif_path, exist_ok=True)
    timestamp               = re.search(r'(\d{4}_\d{1,2}_\d{1,2}_\d{1,2}_\d{1,2}_\d{1,2})', optimal_route_filepath).group(1)

    shell_color = {}
    sat_count = 0

    for shell in sat_config['shells']:
        sat_count += sat_config['shells'][shell]["sat_per_orbit"] * sat_config['shells'][shell]["orbits"]
        shell_color[sat_count] = config['shells'][shell]

if __name__ == '__main__':
    setConfigurations()
    parsing_files('tle',operator_name.upper())
    parsing_files('gs',operator_name.upper())
    parsing_files('node_index',operator_name.upper())
    sat_orbit_index, sat_index_orbit = index_orbit_relation()
    # print(sats_from_tle_dict)

    if not make_gif:
        # Run the processing and plotting
        optimal_route, count, unique_orbits = debugging_section(sat_orbit_index, sat_index_orbit)
        
        # Plot Interactive
        final_plotting(optimal_route, count, unique_orbits)
    else:
        gif_creator()

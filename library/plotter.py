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

import sys
pwd = os.getcwd()
sys.path.append("../")
from mobility.lunar_dyn_utils import *
os.chdir(pwd)

import re
import numpy as np
from datetime import datetime, timezone
from mpl_toolkits.basemap import Basemap



# ================================================================================================
# >>> SCRIPT CONTROL - EDIT HERE <<<
# ================================================================================================
time_index              = 0
plot_GSs                = True
plot_only_optimal       = False
show_optimal            = True   #If false it would not plot the optimal path (doesn't work if make_gif = True)
plot_in_3D              = True
plot_debug              = False   #Plots linked sats to the given sat index and their respective orbits AND also darkcyan satellites that have more than 4 ISLs
plot_optimal_orbits     = False   #Plots all the orbits involved in the optimal path
make_gif                = True   # Makes gif of all the timestep plots in this file existyin in current output path  (REQUIREMENTS: CONNECTIVITY FILES AND OPTIMAL_PATH FILES SHOULD BE EXISTING AND SEPERATE FILES FOR EACH TIMESTEP | line 153 hardcode should be rechecked)
lon0_3d                 = -40 #30  
lat0_3d                 = 50 #-35      
ll                      = [0.4, 0.4]   #scaling of the 3D plot (non-negetive) (lower left point) [0.7, 0.7]
ur                      = [0.4, 0.4]   #scaling of the 3D plot (non-negetive) (upper right point) [0.7, 0.7]
timestamp               = "2025_07_05_06_00_40"
tle_unix_timestamp      = "1751837425"
outputfolder_name       = "Acta/case2/congestion/test1/"
operator_name           = "starlink"
gif_name                = 'gif_test'
number_of_orbits        = 78
_timespan               = 300  #1800 #Important to change if doesnt match Phase1 settings
gs_filepath             = open('/home/spacenet/simulator/gitlab/dynamic-topology-generator/output/'+outputfolder_name+'/terrestrial_info/terrestrial_'+tle_unix_timestamp+'.txt', 'r')
tle_file                = open('/home/spacenet/simulator/gitlab/dynamic-topology-generator/utils/'+operator_name+'_tles/'+operator_name+'_'+tle_unix_timestamp, 'r')
optimal_route_filepath  = '/home/spacenet/simulator/gitlab/dynamic-topology-generator/output/'+outputfolder_name+'/optimal_routes/'+operator_name+'/best_path_'+timestamp+'.0.txt'
conn_filepath           = '/home/spacenet/simulator/gitlab/dynamic-topology-generator/output/'+outputfolder_name+'/connectivity/'+operator_name+'/topology_'+timestamp+'.0.txt'
node_indices_filepath   = '/home/spacenet/simulator/gitlab/dynamic-topology-generator/output/'+outputfolder_name+'/node_indices/'+operator_name+'/nodeindex_'+tle_unix_timestamp+'.txt'
topo_graph_filepath     = '/home/spacenet/simulator/gitlab/dynamic-topology-generator/output/'+outputfolder_name+'/topology_graph/'+operator_name+'/topology_graph_'+timestamp+'.0.txt'
orb_sat_txt             = '/home/spacenet/simulator/gitlab/dynamic-topology-generator/output/'+outputfolder_name+'/satellites_orbits/orbits_satellites.txt'
################### Folder Paths ######################
conn_folder             = '/home/spacenet/simulator/gitlab/dynamic-topology-generator/output/'+outputfolder_name+'/connectivity/'+operator_name+'/'
opt_route_folder        = '/home/spacenet/simulator/gitlab/dynamic-topology-generator/output/'+outputfolder_name+'/optimal_routes/'+operator_name+'/'
gif_path                = '/home/spacenet/GIFS/Acta/case2/congestion/test1/'
#######################################################
Main_body = 'Moon'
Third_body = 'Earth'
ref = 1 #1263  # index of satellite to be debugged for ISLs (Only use with plot_debug=True)
shell_color = {1584:"orange",1814:"green"}
#shell_color = {440:"orange"}
default_projection_for_moon = "ortho"

###### To plot only terrestrial node: plot_only_optimal=True, show_optimal=False, plot_GSs=True, make_gif=False


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
                tle_sat_obj.name                = re.match(r"STARLINK-\d+", tle_sat_name).group(0)
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

    for k in range(max(orbits)+1):
        for sat_in_orb in sat_orbit_index[k]:
            aliass = node_index_to_alias_topology_dict[sat_in_orb]
            if operator_name.upper()=='STARLINK':  #Skyfield type
                sat_at_t    = sats_from_tle_dict[aliass].at(t_current)
                lat, lon    = sat_at_t.subpoint().latitude.degrees, sat_at_t.subpoint().longitude.degrees
                plotted_sat_index[sat_in_orb] = (lat, lon)
            else:  #CustomSatellite type
                sat_obj     = sats_from_tle_dict[aliass]
                sat_ephem_at_epoch = sat_obj.get_ephem()
                sat_at_epoch_r, sat_at_epoch_v = sat_ephem_at_epoch.rv(t_current)
                rr_PA = frame_conversions(t_current, sat_at_epoch_r, 'ICRS', 'PA')
                sat_at_epoch_sph = frame_conversions(t_current, rr_PA, 'PA', 'MER', format='spherical')
                plotted_sat_index[sat_in_orb] = (sat_at_epoch_sph[1], sat_at_epoch_sph[2])
            lats.append(plotted_sat_index[sat_in_orb][0])
            lons.append(plotted_sat_index[sat_in_orb][1])
    return lats, lons


def debugging_section(sat_orbit_index, sat_index_orbit):

    #$ ================================================================================================
    # FILE PARSING - CONNECTIVITY FILE (Checking only sat links)
    #$ ================================================================================================
    with open(conn_filepath, 'r') as conn_file:
        for i, conn_index in enumerate(conn_file):
            line = conn_index.split(",")
            if int(line[0]) < total_num_sat and int(line[1]) < total_num_sat:       ###### Both endpoints should be satellites (exclude all the GSL links) 
                if line[0] not in conn_mat.keys():
                    conn_mat[line[0]] = [int(line[1])]
                else:
                    conn_mat[line[0]].append(int(line[1])) # Actually a nested list

    count = 0
    for i, js in conn_mat.items():
        num_links.append(len(js))
        if len(js)>4:
            count += 1

    # ================================================================================================
    # FILE PARSING - OPTIMAL PATH ASSIGNMENT
    # ================================================================================================
    with open(optimal_route_filepath, 'r') as optimal_path_file:
        for route_line in optimal_path_file:  #Only one line so only one iteration

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
                if a_node[0] == operator_name.upper():
                    optimal_route_satonly_at_epoch.append(route_node_indices[i])

            # Append to complete list
            optimal_routes.append(optimal_route_at_epoch)

            # (Debugging purpose) List of orbits used for optimal path 
            optimal_orbits = [sat_index_orbit[int(k)][0] for k in optimal_route_satonly_at_epoch]
            optimal_orbits = np.unique(optimal_orbits)

    # ================================================================================================
    # DEFINE CURRENT TIME INDEX
    # ================================================================================================
    global t_current
    if operator_name.upper() == 'STARLINK':
        t_current   = ts.from_datetime(dt_hist)
    elif operator_name.upper() == 'LUNAR':
        t_current   = epoch_hist

    #$ Filter sat's (lat,long) orbit-wise as given in nodeindices (rearranged from arranged sats in nodeindex)
    global lats, lons
    lats, lons = figuring_latlon(t_current, optimal_orbits, sat_orbit_index)

    # ================================================================================================
    # CREATE DICTIONARY WITH SATELLITE GEODETIC POSITION
    # ================================================================================================
    global node_info_topology_at_t
    node_info_topology_at_t = node_topology_creation(t_current)

    return optimal_route_at_epoch, count


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


def sat_color_scheme(shell_color, topo_graph_path, total_sats):
    ###### Coloring scheme only for satellites

    sat_color_map = [0 for i in range(total_sats)]

    level_vals = []
    shell_span = 0
    shell_list = list(shell_color.keys())

    with open(topo_graph_path, 'r') as topology_graph:
        for i, graph_links in enumerate(topology_graph):
            line = graph_links.split("\t\t\t\t\t\t")
            indices = line[0].split(",")
            if int(indices[0]) >= total_sats:
                break
            if float(line[-1]) not in level_vals:
                level_vals.append(float(line[-1]))
            
            #### defaulting normal color scheme for different shells
            while int(indices[0])>shell_list[shell_span]:
                shell_span += 1 
            sat_color_map[int(indices[0])] = shell_color[shell_list[shell_span]]
    
    ### Removing no congestion links
    level_vals.remove(float(1))

    if not level_vals:  #No congestion in the entire topology
        return sat_color_map
    else:
        #### Color coding
        color_gradient = [(0.5+0.5*i/len(level_vals),0,0) for i in range(len(level_vals))]
        layer_tracker = [0 for i in range(total_sats)]
        with open(topo_graph_path, 'r') as topology_graph:
            for i, graph_links in enumerate(topology_graph):
                line = graph_links.split("\t\t\t\t\t\t")
                indices = line[0].split(",")
                if int(indices[0]) >= total_sats:
                    break
                if float(line[-1]) != float(1): #Only congestion relevant sats
                    if  level_vals.index(float(line[-1])) > layer_tracker[int(indices[0])]: #Origin sat (coloring based on the link with highest congestion)
                        idx = level_vals.index(float(line[-1]))
                        layer_tracker[int(indices[0])] = idx
                        sat_color_map[int(indices[0])] = color_gradient[idx]
        return sat_color_map


def final_plotting(optimal_route_at_epoch, count):

    global lats, lons

    fig = plt.figure()
    font = {'family' : 'monospace', 
            'size' : 12}
    plt.rc('font', **font)

    # PLOT BASEMAP
    m = basemap_settings()

    coloring_shell_sats = sat_color_scheme(shell_color, topo_graph_filepath, total_num_sat)

    # PLOT ALL SATELLITE NODES IN TOPOLOGY (IF not make_gif then just plot OTHERWISE make the GIF OFC!)
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

                node_idx = node_alias_to_index_topology_dict[node_assigned_alias]

                x, y = m(node_lon, node_lat)
                if node_assigned_alias == node_index_to_alias_topology_dict[ref] and plot_debug:
                    plt.scatter(x, y, s=50, marker="o", facecolors='red', edgecolors='red', zorder=20)
                    plt.text(x, y-0.5, ref, fontsize=15, color='red', zorder=100)
                elif num_links[node_alias_to_index_topology_dict[node_assigned_alias]]>4 and plot_debug:
                    plt.scatter(x, y, s=20, marker="o", facecolors='darkcyan', edgecolors='darkcyan', zorder=20)
                    plt.text(x, y-0.5, node_alias_to_index_topology_dict[node_assigned_alias], fontsize=7, zorder=100)
                else:
                    plt.scatter(x, y, s=20, marker="o", facecolors=coloring_shell_sats[node_idx], edgecolors=coloring_shell_sats[node_idx], zorder=20)
                    plt.text(x, y-0.5, node_alias_to_index_topology_dict[node_assigned_alias], fontsize=7, zorder=100)

    ####################  #$ PLOT SATS AND ITS LINKS WITH HIGHLIGHTED ORBITS (DEBUGGING ZONE STARTS) ###########################
    if plot_debug:
        print(ref, node_index_to_alias_topology_dict[ref])

        ###### Get unique list of relevant orbits
        debug_orbits = [sat_index_orbit[int(k)][0] for k in conn_mat[str(ref)]]
        debug_orbits = np.unique(debug_orbits)

        ###### Lat-Lon of all satellites in the relevant orbits
        lats, lons = figuring_latlon(t_current, debug_orbits, sat_orbit_index)

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
                plt.text(x, y-0.5, neighbours, fontsize=15, color='green', zorder=100)
            else:
                plt.scatter(x, y, s=30, marker="o", facecolors='blue', edgecolors='blue', zorder=20)
                plt.text(x, y-0.5, neighbours, fontsize=15, color='blue', zorder=100)

    ######################################################### DEBUGGING ZONE ENDS ###################################

    # PLOT ALL GROUND STATIONS
    optimal_route_at_t          = optimal_route_at_epoch
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
                # plt.text(x-300000, y+50000, node_assigned_alias, fontsize=9, color='purple', zorder=100)
    handles, labels = plt.gca().get_legend_handles_labels()
    sat_marker = mlines.Line2D([], [], c='black', markerfacecolor='none', markersize=6, label='Satellite', marker='o', linestyle='None')
    gs_marker = mlines.Line2D([], [], c='purple', markerfacecolor='none', markersize=6, label='Ground Station (GW, CT, IE)', marker='p', linestyle='None')

    # PLOT OPTIMAL ROUTE
    if show_optimal:
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
    print('FW Algorithm: '+str(total_num_sat)+' nodes (time: '+str(dt_hist)+') (hops='+str(len(optimal_route_at_t)-1)+')')
    #plt.xlabel('Longitude')
    #plt.ylabel('Latitude')
    #plt.title('Timestep: ' + timestamp + ' |  # of Hops: ' + str(len(optimal_route_at_epoch)-1) + ' | # of involved orbits: ' + str(len(optimal_orbits)) + ' | Non "+ grid" sats: ' + str(count))
    #plt.legend(fancybox=True, framealpha=1, handles=handles, labels=labels, loc='upper left').set_zorder(100)
    plt.tight_layout()
    plt.show()


def gif_creator():

    global node_info_topology_at_t

    conn_mat_global = {}   # Only exists if make_gif exists
    COUNT_global = []      # Only exists if make_gif exists
    num_links_global = {}  # Only exists if make_gif exists
    conn_sorted_path = sorted(os.listdir(conn_folder))
    ######################################### HARDCODED FOR 2024_09_27 FILES #########################################
    # temp = conn_sorted_path[6:9]
    # temp.extend(conn_sorted_path)
    # conn_sorted_path = temp
    # conn_sorted_path.pop(-1)
    # conn_sorted_path.pop(-1)
    # conn_sorted_path.pop(-1)
    # conn_sorted_path.pop(-1)

    # conn_sorted_path.insert(0, conn_sorted_path[5])
    # conn_sorted_path.pop(6)
    # conn_sorted_path.insert(6, conn_sorted_path[-1])
    # conn_sorted_path.pop(-1)

    # conn_sorted_path.insert(6, conn_sorted_path[11])
    # conn_sorted_path.pop(12)
    # conn_sorted_path.insert(12, conn_sorted_path[17])
    # conn_sorted_path.pop(18)
    # conn_sorted_path.insert(18, conn_sorted_path[23])
    # conn_sorted_path.pop(24)
    # conn_sorted_path.insert(24, conn_sorted_path[29])
    # conn_sorted_path.pop(30)
    ##################################################################################################################

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
    ######################################### HARDCODED FOR 2024_09_27 FILES #########################################
    # temp = optroute_sorted_path[6:9]
    # temp.extend(optroute_sorted_path)
    # optroute_sorted_path = temp
    # optroute_sorted_path.pop(-1)
    # optroute_sorted_path.pop(-1)
    # optroute_sorted_path.pop(-1)
    # optroute_sorted_path.pop(-1)
    
    # optroute_sorted_path.insert(0, optroute_sorted_path[5])
    # optroute_sorted_path.pop(6)
    # optroute_sorted_path.insert(6, optroute_sorted_path[-1])
    # optroute_sorted_path.pop(-1)

    # optroute_sorted_path.insert(6, optroute_sorted_path[11])
    # optroute_sorted_path.pop(12)
    # optroute_sorted_path.insert(12, optroute_sorted_path[17])
    # optroute_sorted_path.pop(18)
    # optroute_sorted_path.insert(18, optroute_sorted_path[23])
    # optroute_sorted_path.pop(24)
    # optroute_sorted_path.insert(24, optroute_sorted_path[29])
    # optroute_sorted_path.pop(30)
    ##################################################################################################################
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

        coloring_shell_sats = sat_color_scheme(shell_color, topo_graph_filepath, total_num_sat)

        # PLOT ALL SATELLITE NODES IN TOPOLOGY (IF not make_gif then just plot OTHERWISE make the GIF OFC!)
        if not plot_only_optimal:

            for node_alias, node_info in node_info_topology_at_t.items():

                # Extract information
                node_assigned_alias     = node_info[0]
                node_lon, node_lat      = node_info[1:]

                # Plot ONLY satellite node as a regular scatter point with label
                if not any(gs_type in node_alias for gs_type in gs_alias_list):

                    node_idx = node_alias_to_index_topology_dict[node_assigned_alias]

                    x, y = m(node_lon, node_lat)
                    if node_assigned_alias == node_index_to_alias_topology_dict[ref] and plot_debug:
                        plt.scatter(x, y, s=50, marker="o", facecolors='none', edgecolors='red', zorder=20)
                        #plt.text(x, y-0.5, ref, fontsize=15, color='red', zorder=100)
                    elif num_links[node_alias_to_index_topology_dict[node_assigned_alias]]>4 and plot_debug:
                        plt.scatter(x, y, s=20, marker="o", facecolors='darkcyan', edgecolors='darkcyan', zorder=20)
                        #plt.text(x, y-0.5, node_alias_to_index_topology_dict[node_assigned_alias], fontsize=7, zorder=100)
                    else:
                        plt.scatter(x, y, s=20, marker="o", facecolors=coloring_shell_sats[node_idx], edgecolors=coloring_shell_sats[node_idx], zorder=20)

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


if __name__ == '__main__':
    parsing_files('tle',operator_name.upper())
    parsing_files('gs',operator_name.upper())
    parsing_files('node_index',operator_name.upper())
    sat_orbit_index, sat_index_orbit = index_orbit_relation()

    # print(sats_from_tle_dict)

    if not make_gif:
        optimal_route_at_epoch, count = debugging_section(sat_orbit_index, sat_index_orbit)
        final_plotting(optimal_route_at_epoch, count)
    else:
        gif_creator()
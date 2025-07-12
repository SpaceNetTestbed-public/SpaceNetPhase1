''''

SpaceNet: Generalized Mobility tracking for different celestial bodies (Originally made for Moon)

AUTHOR:         Suryansh Aryan, 2025
                Virginia Tech


This Python script acts as a mobility script for custom celestial bodies (non-Earth). This script also acts as a bridge between 
mobility_util dedicated for propagation for Earth satellites. Creates CustomSatellites class that is analogous to Skyfield Earthsatellites class!
'''

import numpy as np
from astropy import units as u
import re

from poliastro.core.propagation import func_twobody
from poliastro.core.perturbations import J2_perturbation, third_body
import poliastro.bodies as pbodies
from poliastro.bodies import Earth, Sun, Moon
from poliastro.twobody import Orbit
from poliastro.twobody.propagation import CowellPropagator

from poliastro.twobody.sampling import EpochsArray, TrueAnomalyBounds, EpochBounds
from poliastro.util import time_range, Time

from astropy.coordinates import solar_system_ephemeris, SphericalRepresentation, get_body, cartesian_to_spherical, spherical_to_cartesian
solar_system_ephemeris.set("jpl")
from astropy.time import Time, TimeDelta
from poliastro.ephem import Ephem, build_ephem_interpolant
import time

from poliastro.plotting import OrbitPlotter2D, OrbitPlotter3D
import plotly.io as pio  
import matplotlib.pyplot as plt

import pyshtools as pysh


import math
import threading
import sys
sys.path.append("../")
from link.link_utils import *

################## Put the naming terminologies and assumptions on another script so that the entire testbed can refernce that (FIX THIS) 
operator_namechart = {
                            'Earth':'starlink',
                            'Moon' :'lunar'
                            }
######################################################################################
# LIST OF GLOBAL VARIABLES IN USE: (init_epoch, body_r, perturber, pck_data, lunar_harmonics_arr, max_degrees_considered, )
body_r = None

class CustomSatellites:
    """Contains useful information about the satellite in question

    Parameters
    ----------
    name : str
        Name of the satellite.
    data : poliastro.twobody.Orbit type
        All information about the initial orbit of the satellite.
    epoch : poliastro.util.Time type
        Initial epoch of the satellite (format=jd, scale=tdb).  (Julian date is the default)
    main_body : str
        Name of the host celestial body (for Keplerian orbit) [Supports: Earth, Moon, Mars, Sun]
    third_body : str
        Name of the third body (as perturbation) [Supports: Earth, Moon, Mars, Sun]

    """

    def __init__(self, name, data, main_body='Moon', third_body='Earth'):
        self.name = name
        self.data = data
        self.epoch = data.epoch
        self.main_body = self.get_attractor_obj(main_body)
        self.third_body = self.get_perturber_obj(third_body)

        # Private variables
        self.__ephem = None
    
    def get_ephem(self):
        if not self.__ephem:
            print('Cant provide Ephermedis, since it doesnt exist! ')
        else:
            return self.__ephem
    
    def get_OEs(self, units):

        oes = self.data.classical()

        if units=='rad':
            oes[2] = oes[2]*np.pi/180
            oes[3] = oes[3]*np.pi/180
            oes[4] = oes[4]*np.pi/180
            oes[5] = oes[5]*np.pi/180
        
        return oes
    

    def get_latlon_from_ephem(self, ep):

        ephem_sat = self.get_ephem()
        rr, vv = ephem_sat.rv(ep)
        sat_sph = cartesian_to_spherical(rr[0],rr[1],rr[2])
        return sat_sph[1].to_value(u.deg), sat_sph[2].to_value(u.deg)


    def get_initial_kinematics(self, kin_type='cartesian'):
        from astropy.coordinates import SphericalRepresentation
        
        if kin_type=='cartesian':
            return self.data.rv()
        elif kin_type=='spherical':
            _, vel = self.data.rv()
            lat = self.data.represent_as(SphericalRepresentation).lat.value     #rad
            lon = self.data.represent_as(SphericalRepresentation).lon.value     #rad
            r = self.data.represent_as(SphericalRepresentation).distance.value  #same units as self.data
            return np.array(r,lat,lon), vel.value
        
    def get_attractor_obj(self, main_body):

        if main_body=='Moon':
            return pbodies.Moon
        elif main_body=='Earth':
            return pbodies.Earth
        elif main_body=='Sun':
            return pbodies.Sun
        elif main_body=='Mars':
            return pbodies.Mars

    def get_perturber_obj(self, third_body):

        if third_body=='Moon':
            return  pbodies.Moon
        elif third_body=='Earth':
            return  pbodies.Earth
        elif third_body=='Sun':
            return  pbodies.Sun
        elif third_body=='Mars':
            return  pbodies.Mars
    
    def get_body_str(self):

        if self.main_body==pbodies.Moon:
            return 'Moon'
        elif self.main_body==pbodies.Earth:
            return 'Earth'
        elif self.main_body==pbodies.Sun:
            return 'Sun'
        elif self.main_body==pbodies.Mars:
            return 'Mars'
    
    def convert_epoch(self, to_format):
        """
        Outputs required format of the Epoch

        Args:
            to_format (dict): Desirable format

        Returns:
            epoch in desirable format
        """

        if to_format=='datetime':
            date_time = []
            self.epoch.format = 'datetime'
            time_string = self.epoch.strftime('%Y-%d-%d %H:%M:%S')
            data = time_string.split(' ')
            date = data[0].split('-')
            date = list(map(int, date))
            day = data[1].split(':')
            day = list(map(int, day))
            date_time.extend(date)
            date_time.extend(day)
            return date_time
        
        elif to_format=='unix':
            self.epoch.format = "unix"
            unix = self.epoch.value
            return unix 
        
    def construct_ephem(self, tp, rtol=2e-14):
        """
        Outputs required format of the Epoch

        Args:
            tp (astropy.Time.TimeDelta object): How far does the simulation has to ppropagate 

        Returns:
            self._ephem (poliastro.ephem.Ephem object): Ephermedis data for the satellite propagating tp time
        """
        self.__ephem = self.data.to_ephem(
                            EpochsArray(self.data.epoch + tp, method=CowellPropagator(rtol=rtol, f=total_dynamics)),
                            )
            
########################################################################################################################################        
    

def check_gsl_connection(
                        current_gs,
                        current_sat,
                        time,
                        main_config
                        ):
    """
    Calculates the maximum Ground Station-to-Satellite Link (GSL) length

    Args:
        main_configurations (dict): simulation definitions from the YAML configuration file

    Returns:
        max_gsl_length_m (float): maximum gs-sat link length (in meters)
    """
    
    # Initialize return variable 
    check = False

    main_body = current_sat.main_body
    sat_vec = get_sat_MERxyz(current_sat, time)
    gs_vec = get_gs_MERxyz(current_gs, main_body)
    distance_m = np.linalg.norm(sat_vec - gs_vec)

    # Calculate gs within satellite's range
    theta = np.arccos(np.dot(gs_vec - sat_vec, -sat_vec)/(np.linalg.norm(sat_vec)*distance_m))

    # Calculate same satellite within this gs's range
    eps = np.arccos(np.dot(sat_vec - gs_vec, gs_vec)/(np.linalg.norm(gs_vec)*distance_m))

    if theta<=float(main_config["satellite_FOV"])*np.pi/180 and eps>=float(main_config["min_elevation_angle"])*np.pi/180:
        check = True
    else:
        check = False
    
    return check, distance_m


def calc_elfo_gs_sat_thread(
                                ground_stations, 
                                satellites_by_name, 
                                satellites_by_index, 
                                time_t, 
                                main_config, 
                                sat_config, 
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
            is_visible, distance_m = check_gsl_connection(gs, satellites_by_name[str(satellites_by_index[sid])], time_t, main_config)
            
            # Check if the calculated distance is within the maximum GSL length
            if is_visible:
                # If in range, append a tuple to the result list
                ground_station_satellites_in_range.append((distance_m, sid, gs["gid"]))

    # Return the list of valid ground station-satellite pairs
    return ground_station_satellites_in_range


def get_gs_MERxyz(
                    ground_station,
                    main_body
                    ):

    gs_lat = float(ground_station["latitude_degrees_str"])*np.pi/180   #rad
    gs_lon = float(ground_station["longitude_degrees_str"])*np.pi/180  #rad
    gs_R = main_body.R_mean.to_value(u.m) + ground_station["elevation_m_float"]  #m
    gs_MER_xyz = list(spherical_to_cartesian(gs_R, gs_lat, gs_lon))
    gs_MER_xyz = [ele.value for ele in gs_MER_xyz] #m
    return gs_MER_xyz


def get_sat_MERxyz(
                    satellite,
                    current_epoch
                    ):

    sat_ephem = satellite.get_ephem()
    rr, vv = sat_ephem.rv(current_epoch)
    rr = rr.to_value(u.m)   #Cartesian ICRS
    rr_PA = frame_conversions(current_epoch, rr, 'ICRS', 'PA')
    sat_MER_xyz = frame_conversions(current_epoch, rr_PA, 'PA', 'MER')
    return sat_MER_xyz


def store_sat_xyzcoords(
                    satellites_by_name,
                    satellites_by_index,
                    current_time
                    ):
    
    coord_x = []
    coord_y = []
    coord_z = []
    for sid in range(len(satellites_by_index)):
        sat_vec = get_sat_MERxyz(satellites_by_name[str(satellites_by_index[sid])], current_time)
        coord_x.append(sat_vec[0])
        coord_y.append(sat_vec[1])
        coord_z.append(sat_vec[2])
    coords = [coord_x, coord_y, coord_z]
    return coords


def distance_between_ground_station_satellite(
                                              ground_station, 
                                              satellite, 
                                              current_epoch
                                              ):
    """
    Calculates the distance between a ground station and a satellite in a suitable rotating frame (Lunar: MER frame) at a specific time

    Args:
        ground_station (dict): ground station (Moon: Location in MER frame)
        satellite (CustomSatellite object): All data about satellite   (Moon: Ephem Location in ICRS frame)
        current_epoch (Time type): Epoch corresponding to current satellite position

    Returns:
        distance (float): distance between a ground station and a satellite (in meters)
    """
    main_body = satellite.main_body
    
    gs_lat = float(ground_station["latitude_degrees_str"])*np.pi/180   #rad
    gs_lon = float(ground_station["longitude_degrees_str"])*np.pi/180  #rad
    gs_R = main_body.R_mean.to_value(u.m) + ground_station["elevation_m_float"]  #m
    gs_MER_xyz = list(spherical_to_cartesian(gs_R, gs_lat, gs_lon))
    gs_MER_xyz = [ele.value for ele in gs_MER_xyz] #m

    sat_ephem = satellite.get_ephem()
    rr, vv = sat_ephem.rv(current_epoch)
    rr = rr.to_value(u.m)   #Cartesian ICRS
    rr_PA = frame_conversions(current_epoch, rr, 'ICRS', 'PA')
    sat_MER_xyz = frame_conversions(current_epoch, rr_PA, 'PA', 'MER')

    # Euclidean Distance in Cartesian format (Makes sense for radio and optical comm) 
    distance = np.linalg.norm(sat_MER_xyz - gs_MER_xyz)

    # Return the distance between the ground station and satellite in meters
    return distance


def distance_between_two_satellites(
                                    satellite1, 
                                    satellite2, 
                                    current_epoch
                                    ):
    """
    Calculates the distance between two satellites in a suitable frame (Lunar: ICRS frame) at a specific time (epoch)

    Args:
        satellite1 (CustomSatellites object): satellite class containing orbit and constellation data   (Moon: Location in ICRS frame)
        satellite2 (CustomSatellites object): satellite class containing orbit and constellation data   (Moon: Location in ICRS frame)
        current_epoch (Time type): Epoch corresponding to current satellites position

    Returns:
        distance (float): distance between a ground station and a satellite (in meters)
    """
    
    sat1_ephem = satellite1.get_ephem()
    rr1, vv1 = sat1_ephem.rv(current_epoch)
    rr1 = rr1.to_value(u.m)

    sat2_ephem = satellite2.get_ephem()
    rr2, vv2 = sat2_ephem.rv(current_epoch)
    rr2 = rr2.to_value(u.m)

    distance = np.linalg.norm(rr2-rr1)

    return distance


def satellites_from_tle(tle_path, timespan, main_body, third_body, max_deg=50):

    # name = constellation_namechart[main_body]
    # regex_pattern = name.upper()+"-\d+"
    global init_epoch
    global perturber
    global max_degrees_considered

    max_degrees_considered = max_deg
    sats_from_tle_dict = {}
    
    tle_file = open(tle_path, 'r')
    tle_lines = tle_file.readlines()
    tle_name = tle_file.name.split('/')
    #### Figuring epoch from TLE filename
    tle_timestamp = int(tle_name[-1].split('_')[-1])
    julian_date = ( tle_timestamp / 86400.0 ) + 2440587.5
    init_epoch = Time(str(julian_date),format="jd",scale="tdb")
    #####################################
    #### Converting timespan to TimeDelta format
    tp = TimeDelta(np.linspace(0, timespan*u.second, num=1000))  ## FIX THIS (find a way to parameterize num based on timespan)
    ##################################### 

    #### Constructing CustomSatelltie objects for each satellite in TLE file
    load_pck()
    load_harmonics()
    for i in range(0, len(tle_lines), 3):
        tle_sat_name                    = tle_lines[i]
        #sat_name                        = re.match(r"LUNAR-\d+", tle_sat_name).group(0)
        sat_name                        = tle_sat_name.strip('\n')
        tle_first_line                  = list(line for line in tle_lines[i+1].strip("\n").split(" ") if line)
        tle_second_line                 = list(line for line in tle_lines[i+2].strip("\n").split(" ") if line)
        mean_motion                     = float(tle_second_line[7]) * 2 * np.pi / 86400
        SMA                             = ((pbodies.Moon.k.to_value(u.km**3 / u.s**2) / (mean_motion ** 2)) ** (1. / 3.))*(u.km)  #kms   
        sat_orbit_data                  = construct_orbit_MER(pbodies.Moon, SMA, float("."+tle_second_line[4])*u.one, float(tle_second_line[2])*u.deg, float(tle_second_line[3])*u.deg, float(tle_second_line[5])*u.deg, float(tle_second_line[6])*u.deg, init_epoch) #km, kg, s, deg
        sat_obj                         = CustomSatellites(sat_name, sat_orbit_data, main_body, third_body)
        perturber = sat_obj.third_body   #ASSUMES every sat in the TLE has same third body effect (which should be true)
        third_body_ephem(init_epoch, timespan*u.second, sat_obj.main_body)
        sat_obj.construct_ephem(tp, rtol=1e-8)
        sats_from_tle_dict[sat_name]    = sat_obj
    ######################################

    return sats_from_tle_dict


def load_pck(path='./mobility/data/moon/moon_pa_de421_1900-2050.bpc'):
    ### Only supports Lunar libration tracking via PCK (FIX THIS)

    global pck_moon

    from jplephem.pck import PCK
    pck_moon = PCK.open(path)


def load_harmonics():
    ### Only supports Lunar Harmonics (FIX THIS)

    global lunar_harmonics_arr

    hlm = pysh.SHCoeffs.from_file("./mobility/data/moon/jgl150q1.sha", format='shtools',header=True)  #header=True mandatory for files with first line as header (jgl150q1.sha - Based on LP150Q spherical Harmonics model)
    lunar_harmonics_arr = hlm.to_array()
    zonals_only = True
    if zonals_only == True:
        for i, C_i in enumerate(lunar_harmonics_arr[0]):
            for j, orders_j in enumerate(C_i):
                if j != 0:
                    lunar_harmonics_arr[0][i][j] = 0.e+00  # All Cm,n = 0 for n>0
        
        lunar_harmonics_arr[1] = np.zeros((len(lunar_harmonics_arr[1]),len(lunar_harmonics_arr[1][0])))


def get_current_states(sat, time):
    ### Gives position and velocity w.r.t Current frame (most probably MoonICRS) at timstamp 'time'

    sat_ephem = sat.get_ephem()
    pos, vel = sat_ephem.rv(time)

    return pos.to_value(u.km), vel.to_value(u.km/u.s)

def get_current_epoch(
                    t: float
                     ):

    dt = TimeDelta(t,format='sec')
    current_epoch = init_epoch + dt
    return current_epoch


def get_value(property):

    body = pbodies.Moon
    GM = body.k.to_value(u.km**3 / u.s**2)
    radius = body.R_mean.to_value(u.km)
    if property.upper()=='GM':
        return GM
    elif property.lower()=='radius':
        return radius
    elif property.lower()=='all':
        return GM, radius


def construct_orbit_MER(body, a, ecc, inc, raan, argp, nu, ep):

    orb = Orbit.from_classical(body, a, ecc, inc, raan, argp, nu, ep) #Makes satellite w.r.t MoonICRS (for Earth its EarthITES Earth Equator frame which is correct)
    rr, vv = orb.rv()
    rr_fix = frame_conversions(ep, rr, 'PA', 'ICRS')  #Trick to find the correct positon vector in MoonICRS 
    vv_fix = frame_conversions(ep, vv, 'PA', 'ICRS')  #Trick to find the correct positon vector in MoonICRS 
    del orb
    orb_new = Orbit.from_vectors(Moon, rr_fix, vv_fix, ep)
    return orb_new


def plot_sphere(ax, radius):

    # Make data
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)
    x = radius * np.outer(np.cos(u), np.sin(v))
    y = radius * np.outer(np.sin(u), np.sin(v))
    z = radius * np.outer(np.ones(np.size(u)), np.cos(v))

    # Plot the surface
    ax.plot_surface(x, y, z)


def return_OEs(ephem, main_body, times, flag, ecc_vec_plot):
    global max_degrees_considered
    
    from poliastro.core.elements import rv2coe
    k = main_body.k.to(u.km**3 / u.s**2).value
    rr, vv = ephem.rv()

    SMA = []
    ecc = []
    incl = []
    raans = []
    aops = []
    true_anoms = []
    for r, v in zip(rr.to_value(u.km), vv.to_value(u.km / u.s)):
        SMA.append(rv2coe(k, r, v)[0]/(1-rv2coe(k, r, v)[1]**2))
        ecc.append(rv2coe(k, r, v)[1])
        incl.append(rv2coe(k, r, v)[2]*180/np.pi)
        raans.append(rv2coe(k, r, v)[3]*180/np.pi)
        aops.append(rv2coe(k, r, v)[4]*180/np.pi)
        true_anoms.append(rv2coe(k, r, v)[5]*180/np.pi)

    if flag:
        fig, axs = plt.subplots(5)
        fig.suptitle('Osculating Orbital Elements with Time (SH Degree: '+ str(max_degrees_considered) + ')')
        axs[0].plot(times.value, SMA)
        axs[0].set(xlabel='Time (days)', ylabel='Semi-Major Axis')
        axs[1].plot(times.value, ecc)
        axs[1].set(xlabel='Time (days)', ylabel='Eccentricity')
        axs[2].plot(times.value, incl)
        axs[2].set(xlabel='Time (days)', ylabel='Inclination')
        axs[3].plot(times.value, raans)
        axs[3].set(xlabel='Time (days)', ylabel='RAAN')
        axs[4].plot(times.value, aops)
        axs[4].set(xlabel='Time (days)', ylabel='AOP')
        plt.show()
    
    if ecc_vec_plot:

        ECC_X = ecc*np.cos(np.multiply(aops, np.pi/180))
        ECC_Y = ecc*np.sin(np.multiply(aops, np.pi/180))
        fig2 = plt.figure(2)
        ax = fig2.gca()
        ax.set_xlim([-0.1,0.1])
        ax.set_ylim([-0.1,0.1])
        plt.plot(ECC_X, ECC_Y)
        collision_ecc = 1 - main_body.R_mean.to_value(u.km)/SMA[0]
        circ = plt.Circle((0, 0), collision_ecc, color='k', linestyle='--', fill=False)
        ax.add_patch(circ)
        plt.show()

    return SMA, ecc, incl, raans, aops, true_anoms


def eul_rotations(angle, type):

    if type == "z":
        Rot = np.array([[np.cos(angle), np.sin(angle), 0], 
                            [-np.sin(angle), np.cos(angle), 0],
                            [0, 0, 1]])
    elif type == "x":
        Rot = np.array([[1, 0, 0],
                          [0, np.cos(angle), np.sin(angle)], 
                            [0, -np.sin(angle), np.cos(angle)]])
    elif type == "y":
        Rot = np.array([[np.cos(angle), 0, -np.sin(angle)],
                        [0, 1, 0], 
                        [np.sin(angle), 0, np.cos(angle)]])
    else:
        Rot = None
    
    return Rot


def lunar_harmonics_propagator(t, state, body, init_epoch):
    ##### Every calculation in this block in Units: km, s, kg, deg
    global lunar_harmonics_arr
    global max_degrees_considered
    global pck_moon

    radius = body.R_mean.to(u.km)
    
    dt = TimeDelta(t,format='sec')
    current_epoch = init_epoch + dt
    phi, theta, psi = pck_moon.segments[0].compute(current_epoch.value,0,False)

    Rot_pa_icrs = np.matmul(np.matmul(eul_rotations(-phi,'z'), eul_rotations(-theta,'x')), eul_rotations(-psi,'z'))
    rot_pos = np.matmul(np.linalg.inv(Rot_pa_icrs), state[0:3])
    sph_rotated_pos = cartesian_to_spherical(rot_pos[0],rot_pos[1],rot_pos[2])
    acc = pysh.gravmag.MakeGravGridPoint(lunar_harmonics_arr, body.k.to(u.km**3 / u.s**2).value, radius.value, sph_rotated_pos[0].value, sph_rotated_pos[1].to(u.degree).value, sph_rotated_pos[2].to(u.degree).value, [max_degrees_considered,0,0])
    acc_cartesian_rot = spherical_to_cartesian(acc[0], acc[1], acc[2])
    acc_cartesian = np.matmul(Rot_pa_icrs, acc_cartesian_rot)
    del_state = np.array([0, 0, 0, acc_cartesian[0], acc_cartesian[1], acc_cartesian[2]]) #km/s^2

    return del_state


def solar_radiation_dynamics(t, state, body_r, init_epoch, data, perturber):
    ###### Cannonball model refered from the package Basilisk CU Boulder
    
    dt = TimeDelta(t,format='sec')
    current_epoch = init_epoch + dt

    Cr = data[0]
    Area = data[1]
    mass = data[2]
    SF_au = 1372.5398 #W/m^2
    c = 299792458 #m/s
    AU = 1.496*10**(11)  #m
    p_sr = (SF_au/c) #N/m^2

    #current_epoch = sats[0].epoch + tofs
    sun_GCRS = get_body("sun",current_epoch)
    sun_GCRS.representation_type = "cartesian"  # (x,y,z) in kms
    sun_GCRS_vector = [sun_GCRS.x.value, sun_GCRS.y.value, sun_GCRS.z.value]
    if perturber==Earth:
        central_body_to_Earth = body_r(t)
    r_to_sun = -state[:3] + central_body_to_Earth + sun_GCRS_vector  #km
    
    sun_moonICRS_vector = central_body_to_Earth + sun_GCRS_vector

    if np.dot(state[:3],-sun_moonICRS_vector)<=0:
        shadow_factor = 1
    else:
        theta = np.arccos(np.dot(state[:3],-sun_moonICRS_vector)/(np.linalg.norm(state[:3])*np.linalg.norm(-sun_moonICRS_vector)))
        normal_vec = state[:3]*np.sin(theta)
        if np.linalg.norm(normal_vec) > Moon.R_mean.value:
            shadow_factor = 1
        else:
            shadow_factor = 0

    sc_factor = -Cr*p_sr*Area*AU**2/np.linalg.norm(1000*r_to_sun)**3
    SRP_acc = shadow_factor*(sc_factor/mass)*(1000*r_to_sun)  #m/s^2

    du_srp = np.array([0, 0, 0, SRP_acc[0], SRP_acc[1], SRP_acc[2]])*10**(-3)
    return du_srp


def total_dynamics(t0, state, k):
    global init_epoch
    global body_r
    global perturber

    du_kep = func_twobody(t0, state, k)   #km/s^2
    du_harmonics = lunar_harmonics_propagator(t0,state,Moon,init_epoch)
    data = [1.3,23,14000] #[coeff of reflectivity, Surface Area m^2, Mass kg] # Lunar Gateway properties
    du_srp = solar_radiation_dynamics(t0, state, body_r, init_epoch, data, perturber)
    ax, ay, az = third_body(
        t0,
        state,
        k,
        k_third= perturber.k.to(u.km**3 / u.s**2).value,
        perturbation_body=body_r,
    )   # Accelerations are in km/s^2
    du_ad = np.array([0, 0, 0, ax, ay, az])

    return du_kep + du_harmonics + du_ad + du_srp


def ephem_calc(tofs, satellite):
    """
    Generates the Ephermedis data for the satellties revolving around custom celestial bodies

    Args:
        tofs (DeltaTime type): Time intervals for the propagation from initial epoch
        satellite (Orbit object): satellite ephem   (Moon: Location in ICRS frame)

    Returns:
        ephem (Ephem object): Stores the ephermedis data of the satellite at each time interval propagated from the initial epoch
    """
    global main_body
    global init_epoch

    if init_epoch == None:
        init_epoch = satellite.epoch

    if main_body == None:
        main_body = satellite.attractor

    ephem = satellite.to_ephem(EpochsArray(satellite.epoch + tofs, method=CowellPropagator(rtol=2e-14, f=total_dynamics)))

    return ephem


def third_body_ephem(
                    initial_epoch : Time, 
                    prop_time     : u.quantity.Quantity, 
                    main_body     : pbodies):
    """
    Generates the Ephermedis data for the satellties revolving around custom celestial bodies

    Args:
        initial_epoch   (Time type)       : Initial epoch
        prop_time       (Quantity type)   : propagation time for the entire simulation  (in seconds Astropy unit)
        perturber       (pbodies)         : Perturbing third body for the satellites

    Returns:
        body_r (Interp1d object): Stores the interpolated states of the third body w.r.t main attractor 
        
        (Note: for few custom bodies, poliastro doesnt support reference frames yet and therefore there can be errors for some small bodies taken as attractors or thord body)
    """
    global init_epoch
    global perturber
    global body_r
    
    ###### If body_r already exists then no need to execute the method
    if body_r:
        return

    if init_epoch == None:
        init_epoch = initial_epoch
    
    if main_body == 'Moon':
        main_body = pbodies.Moon
    elif main_body == 'Mars':
        main_body = pbodies.Mars
    elif main_body == 'Sun':
        main_body = pbodies.Sun

    body_r = build_ephem_interpolant(
    perturber,
    1000*prop_time.to_value(u.day) * u.day,  #Inverse Time period btw two intervals (time interval depended on the total time of propagation)
    (initial_epoch.value*u.day, initial_epoch.value*u.day + prop_time.to(u.day)),  #interpolated motion till this end time
    rtol=1e-8,  #relative tolerance
    attractor=main_body
)    ###### format: Interp1d format --> input: Float type (Epoch seconds)   Output: 1D list (XYZ coordinates in kms)


def frame_conversions(epoch, position, from_frame, to_frame, format='cartesian'):
    # This is needed for two primary reasons: First, correlation and inter-communication between lunar satellite
    # and lunar GS in a common frame (selenographic coord system with lat lon logic), Secondly correlation 
    # and inter-communication between Earth sats and Lunar sats (Probably in a common ICRF system)
    """
    Args:
        epoch (Time): Time epoch (Astropy Time type)
        ephem (Ephem): Ephermeris of the satellite (Poliastro Ephem type)
        from_frame (str): Current reference frame for the kinematics
        to_frame (str): Required reference frame for the kinematics

    Returns:
        transformed_pos (np.ndarray of floats): Position vector in the transformed frame in cartesian representation (same units as input, same units as input, same units as input)
        OR
        position_sph (np.ndarray of floats): Position vector in the transformed frame in spherical representation (same units as input, deg, deg)
    """
    global pck_moon

    # position, velocity_carte = ephem.rv(epoch) #Cartesian
    # position = position.to_value(u.km)

    #### ICRF system to PA system (reading PCK file with correct epoch is required)
    phi, theta, psi = pck_moon.segments[0].compute(epoch.value,0,False)

    Rot_pa_icrs = np.matmul(np.matmul(eul_rotations(-phi,'z'), eul_rotations(-theta,'x')), eul_rotations(-psi,'z'))
    if from_frame=='ICRS' and to_frame=='PA':
        transformed_pos = np.matmul(np.linalg.inv(Rot_pa_icrs), position)
    elif from_frame=='PA' and to_frame=='ICRS':
        transformed_pos = np.matmul(Rot_pa_icrs, position)


    #### DE440 PA frame and DE421 MER frame (GIVEN BY A CONSTANT MATRIX INDEPENDENT OF EPOCH (Source-JPL))
    phi = (-0.2785/3600)*np.pi/180
    theta = (-78.6944/3600)*np.pi/180
    psi = (-67.8526/3600)*np.pi/180
    Rot_pa_mer = np.matmul(np.matmul(eul_rotations(phi, 'x'), eul_rotations(theta, 'y')), eul_rotations(psi, 'z'))
    if from_frame == 'PA' and to_frame == 'MER':
        transformed_pos = np.matmul(Rot_pa_mer, position)
    elif from_frame == 'MER' and to_frame == 'PA':
        transformed_pos = np.matmul(np.linalg.inv(Rot_pa_mer), position)

    if format=='spherical':
        position_sph = cartesian_to_spherical(transformed_pos[0], transformed_pos[1], transformed_pos[2]) #Spherical coordinates
        position_sph = [position_sph[0].value, position_sph[1].to_value(u.degree), position_sph[2].to_value(u.degree)]
        return position_sph
    elif format=='cartesian':
        return transformed_pos  #(FORMAT: [X, Y, Z] in desired frame)
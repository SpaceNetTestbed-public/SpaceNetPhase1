#### IMPORTANT NOTES:
#    - Go to _interpolate.py in Scipy site-package (Python3.10) and change "fill_value" argument in __init__ from np.nan --> 'extrapolate' (This fixes the interpolation problems when dealing with high time frame propagation with low rtol) 
#    - build_ephem_interpolant() is important to make interp1d type to propagate one celestial body around another! Moon is not incorporated in frame mapping by default in poliastro,
#       so to add Moon, go to util.py for poliastro and Moon: {Planes.EARTH_EQUATOR: MoonICRS} in _FRAME_MAPPING! import Moon from poliastro.bodies and import MoonICRS from poliastro.frames.equatorial
#    - Installing pykep subtleties --> source make in Windows or WSL and normla pip install in linux
#################################### PACKAGE IMPORTS ##############################################
from astropy import units as u
import numpy as np

from poliastro.core.propagation import func_twobody
from poliastro.core.perturbations import J2_perturbation, third_body
from poliastro.bodies import Earth, Sun, Moon
from poliastro.twobody import Orbit
from poliastro.twobody.propagation import CowellPropagator
from poliastro.frames import Planes

from poliastro.twobody.sampling import EpochsArray, TrueAnomalyBounds, EpochBounds
from poliastro.util import time_range, Time

from jplephem.pck import PCK
pck_moon = PCK.open('data/moon/moon_pa_de421_1900-2050.bpc')

from astropy.coordinates import solar_system_ephemeris, get_body, cartesian_to_spherical, spherical_to_cartesian
solar_system_ephemeris.set("jpl")
from astropy.time import Time, TimeDelta
from poliastro.ephem import Ephem, build_ephem_interpolant
import time

from poliastro.plotting import OrbitPlotter2D, OrbitPlotter3D
import plotly.io as pio  
import matplotlib.pyplot as plt
import datetime 

import pyshtools as pysh
hlm = pysh.SHCoeffs.from_file("data/moon/gggrx_1200a_sha.txt", format='shtools',header=True)  #header=True mandatory for files with first line as header
lunar_harmonics_arr = hlm.to_array()
zonals_only = True
if zonals_only == True:
    for i, C_i in enumerate(lunar_harmonics_arr[0]):
        for j, orders_j in enumerate(C_i):
            if j != 0:
                lunar_harmonics_arr[0][i][j] = 0.e+00  # All Cm,n = 0 for n>0
    
    lunar_harmonics_arr[1] = np.zeros((len(lunar_harmonics_arr[1]),len(lunar_harmonics_arr[1][0])))

max_degrees_considered = 50

epoch = Time("2460676.5",format="jd",scale="tdb")  #1st Jan 2025
####################################################################################################

# r = [-6045, -3490, 2500] << u.km
# v = [-3.457, 6.618, 2.533] << u.km/u.s
# orb = Orbit.from_vectors(Earth, r, v)
num_of_planes = 1
raan_bias = 5   # in degrees
################# Data for lunar satellites in constellation
##### OE from Sergey.et.al
alt = 261  #kms
a = Moon.R_mean.to_value(u.km) + alt << u.km
ecc = 0.01 << u.one
inc = 84 << u.deg
argp = -90 << u.deg
nu = 0.0 << u.deg

sats = []
import lunar_dyn_utils as lunar_dyn
lunar_dyn.load_pck('data/moon/moon_pa_de421_1900-2050.bpc')
for i in range(num_of_planes):
    raan = (raan_bias*u.deg)+(i/num_of_planes)*180*u.deg
    orb = Orbit.from_classical(Moon, a, ecc, inc, raan, argp, nu, epoch)
    rr,vv = orb.rv()
    del orb
    rr_fix = lunar_dyn.frame_conversions(epoch, rr, 'PA', 'ICRS')
    vv_fix = lunar_dyn.frame_conversions(epoch, vv, 'PA', 'ICRS')
    orb_fix = Orbit.from_vectors(Moon, rr_fix, vv_fix, epoch)
    sats.append(orb_fix)
##################
perturber = "Earth"

prop_time = 9400*(2.5/6.677955)*sats[0].period.value*u.second #in seconds (Quantity type)   15*300 1yr propagation
print(prop_time.to_value(u.day))
body_r = build_ephem_interpolant(
    Earth,
    100*365 * u.day,  #Inverse Time period btw two intervals
    (epoch.value*u.day, epoch.value*u.day + prop_time.to(u.day)),  #interpolated motion till this end time
    rtol=1e-8,  #relative tolerance
    attractor=Moon
)    ###### format: Interp1d format --> input: Float type (Epoch seconds)   Output: 1D list (XYZ coordinates in kms)

tofs = TimeDelta(np.linspace(0, prop_time, num=20000*10**(1)))   # jd format
######################### Visualizing third body motion around central body ##############
# fig = plt.figure()
# plt.xlabel("X")
# plt.ylabel("Y")
# ax = plt.axes(projection="3d")
# ax.set_zlabel("Z")
# for p in range(len(tofs)):
#     pos = body_r(tofs[p].to_value('s')) #numpy array
#     ax.scatter(pos[0],pos[1],pos[2])
# plt.show()
##########################################################################################


##########################################################################################
def plot_sphere(ax, radius):

    # Make data
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)
    x = radius * np.outer(np.cos(u), np.sin(v))
    y = radius * np.outer(np.sin(u), np.sin(v))
    z = radius * np.outer(np.ones(np.size(u)), np.cos(v))

    # Plot the surface
    ax.plot_surface(x, y, z)



def return_OEs(ephem, main_body, times, flag, ecc_vec_plot=0):
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
    smallest_dist_from_moon = 100000000
    for r, v in zip(rr.to_value(u.km), vv.to_value(u.km / u.s)):
        SMA.append(rv2coe(k, r, v)[0]/(1-rv2coe(k, r, v)[1]**2))  #km
        ecc.append(rv2coe(k, r, v)[1])
        incl.append(rv2coe(k, r, v)[2]*180/np.pi)  #deg
        raans.append(rv2coe(k, r, v)[3]*180/np.pi)  #deg
        aops.append(rv2coe(k, r, v)[4]*180/np.pi)  #deg
        true_anoms.append(rv2coe(k, r, v)[5]*180/np.pi)  #deg
        if (np.linalg.norm(r)-main_body.R_mean.to_value(u.km)) < smallest_dist_from_moon:
            smallest_dist_from_moon = (np.linalg.norm(r)-main_body.R_mean.to_value(u.km))
            ecc_at_smallest_dist = rv2coe(k, r, v)[1]
    
    SMA_min = min(SMA)
    collision_ecc = 1 - main_body.R_mean.to_value(u.km)/SMA_min
    
    print('closest distance to Moon: ', smallest_dist_from_moon)
    print('Osculating Ecc at closest distance to Moon: ', ecc_at_smallest_dist)
    print('Minimum SMA value: ', SMA_min)
    print('Collision Ecc value: ', collision_ecc)

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
        ax.set_xlim([-0.3,0.3])
        ax.set_ylim([-0.3,0.3])
        ax.set(xlabel='ecos(AOP)', ylabel='esin(AOP)')
        plt.plot(ECC_X, ECC_Y)
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
    if perturber=="Earth":
        central_body_to_Earth = body_r(t)
    r_to_sun = -state[:3] + central_body_to_Earth + sun_GCRS_vector  #km
    
    sun_moonICRS_vector = central_body_to_Earth + sun_GCRS_vector

    if np.dot(state[:3],-sun_moonICRS_vector)<=0:
        shadow_factor = 1
    else:
        theta = np.arccos(np.dot(state[:3],-sun_moonICRS_vector)/(np.linalg.norm(state[:3])*np.linalg.norm(-sun_moonICRS_vector)))
        normal_vec = state[:3]*np.sin(theta)  #kms
        if np.linalg.norm(normal_vec) > Moon.R_mean.to_value(u.km):
            shadow_factor = 1
        else:
            shadow_factor = 0

    sc_factor = -Cr*p_sr*Area*AU**2/np.linalg.norm(1000*r_to_sun)**3
    SRP_acc = shadow_factor*(sc_factor/mass)*(1000*r_to_sun)  #m/s^2

    du_srp = np.array([0, 0, 0, SRP_acc[0], SRP_acc[1], SRP_acc[2]])*10**(-3)  #km/s^2
    return du_srp




def f(t0, state, k):
    global lunar_harmonics_arr

    du_kep = func_twobody(t0, state, k)   #km/s^2
    du_harmonics = lunar_harmonics_propagator(t0,state,Moon,epoch)
    data = [1.3,23,14000] #[coeff of reflectivity, Surface Area m^2, Mass kg] # Lunar Gateway properties
    du_srp = solar_radiation_dynamics(t0, state, body_r, epoch, data, perturber)
    ax, ay, az = third_body(
        t0,
        state,
        k,
        k_third= Earth.k.to(u.km**3 / u.s**2).value,
        perturbation_body=body_r,
    )   # Accelerations are in km/s^2
    du_ad = np.array([0, 0, 0, ax, ay, az])

    return du_kep + du_harmonics + du_ad









EPHEM = []
EPHEM_nopert = []
for sat in sats:
    ephems = sat.to_ephem(
    EpochsArray(sat.epoch + tofs, method=CowellPropagator(rtol=2e-14, f=f)),
    )
    ephems_nopert = sat.to_ephem(
    EpochsArray(sat.epoch + tofs, method=CowellPropagator(rtol=2e-14, f=func_twobody)),
    )
    EPHEM.append(ephems)
    EPHEM_nopert.append(ephems_nopert)

OEs = return_OEs(ephems, Moon, tofs, 1, 1)
#print(OEs[-2])
#plt.plot(tofs.value, OEs[-2])

fig = plt.figure()
plt.xlabel("X")
plt.ylabel("Y")
ax = plt.axes(projection="3d")
ax.set_zlabel("Z")
plot_sphere(ax, Moon.R_mean.to_value(u.km))
for itr, ephems in enumerate(EPHEM):
    pos = []
    X = []
    Y = []
    Z = []
    vel = []
    X_nopert = []
    Y_nopert = []
    Z_nopert = []
    for i in range(len(ephems.epochs)):
        data = ephems.rv(ephems.epochs[i])
        position_tuple = data[0]
        velocity_tuple = data[1]
        X.append(position_tuple[0].value)
        Y.append(position_tuple[1].value)
        Z.append(position_tuple[2].value)
        data_nopert = EPHEM_nopert[itr].rv(EPHEM_nopert[itr].epochs[i])
        position_tuple_nopert = data_nopert[0]
        velocity_tuple_nopert = data_nopert[1]
        X_nopert.append(position_tuple_nopert[0].value)
        Y_nopert.append(position_tuple_nopert[1].value)
        Z_nopert.append(position_tuple_nopert[2].value)
        pos.append(list(position_tuple.value))  #km
        vel.append(list(velocity_tuple.value))  #km/s
    ax.plot3D(X,Y,Z,'green')                     # Perturbation plot
    #ax.plot3D(X_nopert,Y_nopert,Z_nopert,'red')  # No perturbation plot
    ax.set_aspect('equal')

#print([np.linalg.norm(vel_ele) for vel_ele in vel])


############### Plotting
plt.show()





# import sys
# sys.path.insert(0, "/home/suryaryan/pykep/build_pykep")
# import pykep as pk
# r_m, mu_m, c, s, max_degree, max_order = pk.util.load_gravity_model("glgm3_150.txt")

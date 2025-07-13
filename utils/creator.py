''''

SpaceNet: Generate Fake Two-Line Element Sets (TLE) for custom constellation (default: ELFO type)

AUTHOR:         Suryansh Aryan, 2025
                Virginia Tech


This Python script supplies a function to generate Two-Line Element-like file for the Sets of
virtual (fake) satellites with a Non-Earth central body.
'''

# =================================================================================== #
# ------------------------------- IMPORT PACKAGES ----------------------------------- #
# =================================================================================== #

import os
import time
import datetime
import numpy as np
import calendar
from scipy.optimize import fsolve
import generate_fake_TLE as gft
from astropy import units as u
from astropy.time import Time
import poliastro.bodies as pbodies

# =================================================================================== #
# -------------------------------- MAIN FUNCTION ------------------------------------ #
# =================================================================================== #

def EA_func(E, MA, ecc):
    return E - ecc*np.sin(E) - MA

def basic_generate_fake_TLE(
                                operator_name : str,
                                main_body     : str,
                                epoch         : np.ndarray,
                                oe            : np.ndarray,
                                ma_range      : np.ndarray,
                                inc_range     : np.ndarray,
                                raan_range    : np.ndarray,
                                ipp_angle     : float
                           )   -> list[str]:
    """
    Generates virtual satellite TLEs using the user-input epoch, orbital elements (oe), and
    range of values of True Anomaly (MA) and Inclincation. This is a basic function that
    considers only variation in TA and Inclination. In addition, spacecraft collisions are 
    NOT considered. 

    Args:
        epoch (np.ndarray):     Last two digits of epoch year and day of year (with fraction)  -> e.g., [97, 210.21344211]
        oe (np.ndarray):        Keplerian orbital elements (a-km, e, i-deg, aop-deg, raan-deg, ta-deg)
        ma_range (np.ndarray):  Range of Mean Anomaly values to generate virtual satellite TLEs, in deg
        inc_range (np.ndarray): Range of Inclination values to generate virtual satellite TLEs, in deg
        raan_range (np.ndarray): Range of Right Ascension of the Ascending Node values to generate virtual satellite TLEs, in deg
        ipp_angle (float): Inter Plane Phase Angle; in-track spacing angle between the first satellites in adjacent planes
    
    Returns:
        list [str]:             List of TLE formats written as String instances
    """

    # Initialize oe_sweep
    oe_sweep    = np.array([[0., ]*6, ] * len(raan_range) * len(ma_range))

    # Iterate through TA and RAAN sets
    sat_indx = 0
    for RAAN_indx, RAAN in enumerate(raan_range):
        for MA in ma_range:
            MA_x = MA + RAAN_indx*ipp_angle
            EA_func = lambda E : E - oe[1]*np.sin(E) - MA_x*np.pi/180
            EA = fsolve(EA_func, 0.5)
            TA_x = 2*np.arctan(np.sqrt((1+oe[1])/(1-oe[1]))*np.tan(EA/2))
            oe_sweep[sat_indx, :] = oe
            oe_sweep[sat_indx, 4] = RAAN
            oe_sweep[sat_indx, 5] = TA_x*180/np.pi
            sat_indx += 1
    
    # Initialize TLE list
    TLE_output  = [""] * len(oe_sweep)

    # Perform sweep of TLEs
    for indx, oe_line in enumerate(oe_sweep):

        # Generate fake TLE
        TLE_output[indx] = gft.generate_virtual_TLE(operator_name, main_body, epoch=epoch, oe=oe_line, iter_num=indx)

    
    return TLE_output



# =================================================================================== #
# ----------------------------------- RUN SIM --------------------------------------- #
# =================================================================================== #

if __name__ == "__main__":


    ########################### ONLY MAKE CHANGES WITHIN THIS ZONE (FOR NON-DEVELOPERS) #########################################

    ELFO = True  # If false then Walker delta (IMPORTANT: if ELFO is true then the input inc is not considered, inc is calculated based on input eccentricity using ELFO equation)
    Earth = {'tot_num_sats':25, 'num_orbits':5, 'perigee':4500, 'ecc':0.2, 'inc':53, 'argp':90, 'raan_bias':0.0, 'ipp_increment':1}  # If a circular constellation then (perigee:altitude and ecc:0.0)
    Moon = {'tot_num_sats':80, 'num_orbits':10, 'perigee':200, 'ecc':0.130, 'inc':84, 'argp':90, 'raan_bias':0.0, 'ipp_increment':1}

    body = Moon   # Choose the central body you want your TLE for!

    date_time    = [2029, 9, 27, 22, 15, 6] # year, month, day, hour, minute, second
    
    #############################################################################################################################


    if body == Earth:  # Earth TLEs
        R = 6378
        operator_name = "STARLINK"
        main_body = "Earth"
    elif body == Moon:  # Lunar TLEs
        R = pbodies.Moon.R_mean.to_value(u.km)
        operator_name = "LUNAR"
        main_body = "Moon"
    
    perigee             = body['perigee']  #kms
    ecc                 = body['ecc']
    SMA                 = (R+perigee)/(1-ecc)   #Earth
    if ELFO:
        inc = np.arccos((0.6*(1-ecc**2))**(1/2))*180/np.pi  #deg
        if abs(body['argp'])==90 or abs(body['argp'])==270:
            argp = body['argp']
        else:
            argp = np.random.choice([90, 270])
    else:
        inc = body['inc']
        argp = body['argp']
    raan_bias           = body['raan_bias']  #deg
    tot_num_sats        = body['tot_num_sats']
    num_orbits          = body['num_orbits']
    num_sat_per_orbit   = int(tot_num_sats/num_orbits)

    # Date/Time input
    datetime_s  = calendar.timegm(date_time)
    datetime_obj = datetime.datetime(date_time[0],date_time[1],date_time[2],date_time[3],date_time[4],date_time[5])
    datetime_gd = Time(datetime_obj,format='datetime',scale='tdb')
    year_start  = (date_time[0], 1, 1, 0, 0, 0)
    
    # File name
    filename      = operator_name.lower() + "_" + str(datetime_s)

    # Epoch (Julian date)
    datetime_gd.format = 'jd'
    #epoch         = datetime_gd
    epoch         = [(date_time[0]-2000), (datetime_s - calendar.timegm(year_start))/86400]

    # Orbital Elements (a-km, e, i-deg, w-deg, raan-deg, MA-deg)
    oe          = [SMA, ecc, inc, argp, raan_bias, 0.0000000]

    # Range of MA
    ma_range    = np.linspace(0, 360*(1-1/num_sat_per_orbit), num_sat_per_orbit)

    # Range of Inc
    inc_range     = None

    # Range of RAAN
    # raan_range = np.array(range(0, 360, 8)) # number of orbits = 360/range_stepsize
    raan_range  = np.linspace(raan_bias, (raan_bias + 360)*(1-1/num_orbits), num_orbits)

    # Inter Plane Phase Increment/Angle
    ipp_increment = body['ipp_increment'] # set to zero for no IPP Angle, otherwise set to a positive integer
    ipp_angle = ipp_increment*360/(tot_num_sats)
        # Inter Plane Phase Increment pulled from Walker constellation pattern notation, I:T/P/F
            # I: orbital inclination
            # T: Total number of satellites (must be divisible by F)
            # P: Number of equally spaced orbital planes
            # F: Inter Plane Phase Increment
                # Inter Plane Phase Angle = F*360/T



    # -------------------------------------------------------------------------------
    # Generate TLE sweep

    TLE_sweep     = basic_generate_fake_TLE(operator_name, main_body, epoch=epoch, oe=oe, ma_range=ma_range, inc_range=inc_range, raan_range=raan_range, ipp_angle=ipp_angle)

    # -------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------
    # Write to a file

    with open("../utils/"+operator_name.lower()+"_tles/"+filename, 'a') as file:
        for TLE in TLE_sweep:
            file.write(TLE + '\n')
    print("TLE Generated. Count =", len(TLE_sweep), ". No. orbits=", len(raan_range), ". No. sats per orbit=", len(ma_range))
    # -------------------------------------------------------------------------------
    



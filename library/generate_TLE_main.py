''''

SpaceNet: Generate Fake Two-Line Element Sets (TLE)

AUTHOR:         Bruce L. Barbour, 2023
                Virginia Tech
EDITED BY:      Jack Downs, 2024
                Virginia Tech


This Python script supplies a function to generate Two-Line Element Sets of
virtual (fake) satellites and writes them into a file.
'''

# =================================================================================== #
# ------------------------------- IMPORT PACKAGES ----------------------------------- #
# =================================================================================== #

import os
import time
import numpy as np
import calendar
import generate_fake_TLE as gft

# =================================================================================== #
# -------------------------------- MAIN FUNCTION ------------------------------------ #
# =================================================================================== #

def basic_generate_fake_TLE(
                                epoch       : np.ndarray,
                                oe          : np.ndarray,
                                ta_range    : np.ndarray,
                                inc_range   : np.ndarray,
                                raan_range  : np.ndarray,
                                ipp_angle   : float
                           )   -> list[str]:
    """
    Generates virtual satellite TLEs using the user-input epoch, orbital elements (oe), and
    range of values of True Anomaly (TA) and Inclincation. This is a basic function that
    considers only variation in TA and Inclination. In addition, spacecraft collisions are 
    NOT considered. 

    **HAVE NOT ADDED IN INCLINATION RANGE YET. - BRUCE

    Args:
        epoch (np.ndarray):     Last two digits of epoch year and day of year (with fraction)  -> e.g., [97, 210.21344211]
        oe (np.ndarray):        Keplerian orbital elements (a-km, e, i-deg, aop-deg, raan-deg, ta-deg)
        ta_range (np.ndarray):  Range of True Anomaly values to generate virtual satellite TLEs, in deg
        inc_range (np.ndarray): Range of Inclination values to generate virtual satellite TLEs, in deg
        raan_range (np.ndarray): Range of Right Ascension of the Ascending Node values to generate virtual satellite TLEs, in deg
        ipp_angle (float): Inter Plane Phase Angle; in-track spacing angle between the first satellites in adjacent planes
    
    Returns:
        list [str]:             List of TLE formats written as String instances
    """

    # Initialize oe_sweep
    oe_sweep    = np.array([[0., ]*6, ] * len(raan_range) * len(ta_range))

    # Iterate through TA and RAAN sets
    sat_indx = 0
    for RAAN_indx, RAAN in enumerate(raan_range):
        for TA in ta_range:
            TA_x = TA + RAAN_indx*ipp_angle
            oe_sweep[sat_indx, :] = oe
            oe_sweep[sat_indx, 4] = RAAN
            oe_sweep[sat_indx, 5] = TA_x
            sat_indx += 1


    # # Set the Inter Plane Phase Angle to zero if set to None
    # if ipp_angle is None:
    #     ipp_angle = 0

    # # Check if a TA range and RAAN range are provided:
    # # If case where a range of multiple values is given for both TA and RAAN
    # if ta_range is not None and raan_range is not None and len(raan_range) > 1 and len(ta_range) > 1:
    #     oe_sweep = np.resize(oe_sweep, (len(ta_range)*len(raan_range), 6))
    #     oe_sweep_split = np.array_split(oe_sweep, len(raan_range))
    #     for count in range(len(raan_range)):
    #         oe_sweep_split[count][:, 4] = raan_range[count]
    #         for indx, oe_line in enumerate(oe_sweep_split[count]):
    #             oe_line[5] = ta_range[indx] + count*ipp_angle
    #         if count == 0:
    #             oe_sweep = oe_sweep_split[count]
    #         else:
    #             oe_sweep = np.vstack((oe_sweep, oe_sweep_split[count]))
    # # Else-if case where a range of multiple values is given for RAAN but only one value for TA
    # elif ta_range is not None and raan_range is not None and len(raan_range) > 1 and len(ta_range) == 1:
    #     oe_sweep[5] = ta_range[0]
    #     oe_sweep = np.resize(oe_sweep, (len(raan_range), 6))
    #     for indx, oe_line in enumerate(oe_sweep):
    #         oe_line[4] = raan_range[indx]
    # # Else-if case where a range of multiple values is given for TA but only one value for RAAN
    # elif ta_range is not None and raan_range is not None and len(ta_range) > 1 and len(raan_range) == 1:
    #     oe_sweep[4] = raan_range[0]
    #     oe_sweep = np.resize(oe_sweep, (len(ta_range), 6))
    #     for indx, oe_line in enumerate(oe_sweep):
    #         oe_line[5] = ta_range[indx]
    # # Else-if case where a single value is give for RAAN, but TA is set as None
    # elif raan_range is not None and ta_range is None and len(raan_range) == 1:
    #     oe_sweep[4] = raan_range[0]
    # # Else-if case where a single value is give for TA, but RAAN is set as None
    # elif ta_range is not None and raan_range is None and len(ta_range) == 1:
    #     oe_sweep[5] = ta_range[0]


    # Initialize TLE list
    TLE_output  = [""] * len(oe_sweep)

    # Perform sweep of TLEs
    for indx, oe_line in enumerate(oe_sweep):

        # Generate fake TLE
        TLE_output[indx] = gft.generate_virtual_TLE(epoch=epoch, oe=oe_line, iter_num=indx)

    
    return TLE_output



# =================================================================================== #
# ----------------------------------- RUN SIM --------------------------------------- #
# =================================================================================== #

if __name__ == "__main__":


    # Settings
    altitude            = 500.
    inclination         = 50.
    tot_num_sats        = 880
    num_orbits          = 44
    num_sat_per_orbit   = int(tot_num_sats/num_orbits)

    # Date/Time input
    datetime    = (2024, 7, 21, 9, 0, 0) # year, month, day, hour, minute, second
    datetime_s  = calendar.timegm(datetime)
    year_start  = (datetime[0], 1, 1, 0, 0, 0)
    
    # File name
    filename      = "starlink_" + str(datetime_s)

    # Epoch Year / Fractional Day
    epoch         = [(datetime[0]-2000), (datetime_s - calendar.timegm(year_start))/86400]

    # Orbital Elements (a-km, e, i-deg, w-deg, raan-deg, TA-deg)
    oe          = [6378+altitude, 0.00000, inclination, 0.000000, 0.000000, 0.0000000]

    # Range of TA
    # ta_range    = np.array(range(0, 360, 18)) # number of satellites per orbit = 360/range_stepsize
    ta_range    = np.linspace(0, 360*(1-1/num_sat_per_orbit), num_sat_per_orbit)

    # Range of Inc
    inc_range     = None

    # Range of RAAN
    # raan_range = np.array(range(0, 360, 8)) # number of orbits = 360/range_stepsize
    raan_range  = np.linspace(0, 360*(1-1/num_orbits), num_orbits)

    # Inter Plane Phase Increment/Angle
    ipp_increment = 1 # set to zero for no IPP Angle, otherwise set to a positive integer
    ipp_angle = ipp_increment*360/(tot_num_sats)
        # Inter Plane Phase Increment pulled from Walker constellation pattern notation, I:T/P/F
            # I: orbital inclination
            # T: Total number of satellites (must be divisible by F)
            # P: Number of equally spaced orbital planes
            # F: Inter Plane Phase Increment
                # Inter Plane Phase Angle = F*360/T



    # -------------------------------------------------------------------------------
    # Generate TLE sweep

    TLE_sweep     = basic_generate_fake_TLE(epoch=epoch, oe=oe, ta_range=ta_range, inc_range=inc_range, raan_range=raan_range, ipp_angle=ipp_angle)

    # -------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------
    # Write to a file

    with open("../utils/starlink_tles/"+filename, 'a') as file:
        for TLE in TLE_sweep:
            file.write(TLE + '\n')
    print("TLE Generated. Count =", len(TLE_sweep), ". No. orbits=", len(raan_range), ". No. sats per orbit=", len(ta_range))
    # -------------------------------------------------------------------------------
    



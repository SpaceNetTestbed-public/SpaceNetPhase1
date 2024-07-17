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

    # Initialize OE sweep array
    oe_sweep    = oe

    # Set the Inter Plane Phase Angle to zero if set to None
    if ipp_angle is None:
        ipp_angle = 0

    # Check if a TA range and RAAN range are provided:
    # If case where a range of multiple values is given for both TA and RAAN
    if ta_range is not None and raan_range is not None and len(raan_range) > 1 and len(ta_range) > 1:
        oe_sweep = np.resize(oe_sweep, (len(ta_range)*len(raan_range), 6))
        oe_sweep_split = np.array_split(oe_sweep, len(raan_range))
        for count in range(len(raan_range)):
            oe_sweep_split[count][:, 4] = raan_range[count]
            for indx, oe_line in enumerate(oe_sweep_split[count]):
                oe_line[5] = ta_range[indx] + count*ipp_angle
            if count == 0:
                oe_sweep = oe_sweep_split[count]
            else:
                oe_sweep = np.vstack((oe_sweep, oe_sweep_split[count]))
    # Else-if case where a range of multiple values is given for RAAN but only one value for TA
    elif ta_range is not None and raan_range is not None and len(raan_range) > 1 and len(ta_range) == 1:
        oe_sweep[5] = ta_range[0]
        oe_sweep = np.resize(oe_sweep, (len(raan_range), 6))
        for indx, oe_line in enumerate(oe_sweep):
            oe_line[4] = raan_range[indx]
    # Else-if case where a range of multiple values is given for TA but only one value for RAAN
    elif ta_range is not None and raan_range is not None and len(ta_range) > 1 and len(raan_range) == 1:
        oe_sweep[4] = raan_range[0]
        oe_sweep = np.resize(oe_sweep, (len(ta_range), 6))
        for indx, oe_line in enumerate(oe_sweep):
            oe_line[5] = ta_range[indx]
    # Else-if case where a single value is give for RAAN, but TA is set as None
    elif raan_range is not None and ta_range is None and len(raan_range) == 1:
        oe_sweep[4] = raan_range[0]
    # Else-if case where a single value is give for TA, but RAAN is set as None
    elif ta_range is not None and raan_range is None and len(ta_range) == 1:
        oe_sweep[5] = ta_range[0]


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


    # File name
    filename    = "starlink_" + str(int(time.time()))

    # Epoch Year / Fractional Day
    epoch       = [24, 100]

    # Orbital Elements (a-km, e, i-deg, w-deg, raan-deg, TA-deg)
    oe          = [6878, 0.00000, 50.000000, 0.000000, 0.000000, 0.0000000]

    # Range of TA
    ta_range    = np.array(range(0, 360, 18)) # number of satellites per orbit = 360/range_stepsize

    # Range of Inc
    inc_range   = None

    # Range of RAAN
    raan_range = np.array(range(0, 360, 8)) # number of orbits = 360/range_stepsize

    # Total Number of Satellites
    num_sats = len(ta_range)*len(raan_range)

    # Inter Plane Phase Increment/Angle
    ipp_increment = 1 # set to zero for no IPP Angle, otherwise set to a positive integer
    ipp_angle = ipp_increment*360/num_sats
        # Inter Plane Phase Increment pulled from Walker constellation pattern notation, I:T/P/F
            # I: orbital inclination
            # T: Total number of satellites (must be divisible by F)
            # P: Number of equally spaced orbital planes
            # F: Inter Plane Phase Increment
                # Inter Plane Phase Angle = F*360/T



    # -------------------------------------------------------------------------------
    # Generate TLE sweep

    TLE_sweep = basic_generate_fake_TLE(epoch=epoch, oe=oe, ta_range=ta_range, inc_range=inc_range, raan_range=raan_range, ipp_angle=ipp_angle)

    # -------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------
    # Write to a file

    with open("/mnt/c/Users/BluBoy/Desktop/Professional/Git/Repositories/dynamic-topology-generator/utils/starlink_tles/"+filename, 'a') as file:
        for TLE in TLE_sweep:
            file.write(TLE + '\n')
    print("TLE Generated. Count =", len(TLE_sweep), ". No. orbits=", len(raan_range), ". No. sats per orbit=", len(ta_range))
    # -------------------------------------------------------------------------------
    



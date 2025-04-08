''''

SpaceNet: Generate Fake Two-Line Element Sets (TLE)

AUTHOR:         Bruce L. Barbour, 2023
                Virginia Tech


This Python script supplies utility functions for generating Two-Line Element Sets of
virtual (fake) satellites.
'''

# =================================================================================== #
# ------------------------------- IMPORT PACKAGES ----------------------------------- #
# =================================================================================== #

import time
import numpy as np

# =================================================================================== #
# ------------------------------ UTILITY FUNCTIONS ---------------------------------- #
# =================================================================================== #

# ----------------------------------------------------- #
# FUNCTION 1 : CONVERT TO TLE FORMAT                    #
# ----------------------------------------------------- #
def tle_format_field(
                        num         : float,
                        num_base    : int,
                        num_dec     : int,
                        fieldtype   : int
                    )   -> str:
    """
    Formats the content into the strict formatting guidelines for TLEs.

    Args:
        num (float):        TLE content
        num_base (int):     Number of base numbers in format
        num_dec (int):      Number of decimals in format
        fieldtype (int):    The type of content (e.g., 0-standard, 1-epoch_year, 2-ecc)
    
    Returns:
        str:                Format corrected content 
    """

    # Checks if field type is empty, then use standard
    if not fieldtype:
        fieldtype = 0

    # Determines if it's a special field
    if fieldtype == 1:
        return str(num).zfill(num_base)  # Epoch year case
    elif fieldtype == 2:
        return format(num, f".{num_dec}f").split(".")[1]  # Eccentricity case
    
    # Corrects the base number
    base_part = str(num).split(".")[0].rjust(num_base)

    # Correct the decimal number
    decimal_part = format(num, f".{num_dec}f").split(".")[1]

    # Fix in case where inputted float is significantly close to rounding up
        # Issue primarily existed for mean anomaly inputs (e.g. a mean anomaly would 
        # calculate as 89.99999999, but would be stored in the TLE as 89.0000, rather 
        # than 90.0000)
    if (num-float(base_part)+0.0001) >= 1:
        base_part = float(base_part) + 1
        base_part = str(base_part).split(".")[0].rjust(num_base)


    return f"{base_part}.{decimal_part}"



# ----------------------------------------------------- #
# FUNCTION 2 : GENERATE TLE                             #
# ----------------------------------------------------- #
def generate_virtual_TLE(
                            epoch       : np.ndarray,
                            oe          : np.ndarray,
                            iter_num    : int
                        )   -> str:
    """
    Generates fake TLE(s) using provided input parameters of the satellite(s). It assumes
    that the same epoch will be assigned to each TLE. In addition, the only required inputs 
    of the TLE are:

    (1) Epoch Year
    (2) Epoch Fractional Day
    (3) Semi-major Axis
    (4) Eccentricity
    (3) Inclination
    (4) Argument of Perigee
    (5) Right Ascension of Ascending Node
    (6) True Anomaly

    The rest of the TLE content are randomly generated as they aren't important for the SpaceNet 
    Mininet simulator.

    Args:
        epoch (np.ndarray):     Last two digits of epoch year and day of year (with fraction)  -> e.g., [97, 210.21344211]
        oe (np.ndarray):        Keplerian orbital elements (a-km, e, i-deg, aop-deg, raan-deg, ta-deg)
        iter_num (int):         Iteration number, if running through an iterative process

    Returns:
        str:                    Spacecraft orbit properties in a TLE format
    """

    # Format the epoch into String
    epoch_str       = tle_format_field(epoch[0], 2, 0, 1) + tle_format_field(epoch[1], 3, 8, None)
    
    # Compute mean motion in rev/day
    earth_mu        = 398600.435507     # Earth's gravitational parameter, km2s-2
    nbar            = np.sqrt(earth_mu/oe[0]**3) * (86400./(2*np.pi))    # rev/day

    # Compute mean anomaly in deg
    TA              = np.deg2rad(oe[5])
    MA              = np.rad2deg(np.arctan2(-np.sqrt(1 - oe[1]**2)*np.sin(TA), -oe[1] - np.cos(TA)) + np.pi - oe[1]*(np.sqrt(1 - oe[1]**2) * np.sin(TA))/(1 + oe[1]*np.cos(TA)))

    # Format orbital elements into String
    nbar_str        = tle_format_field(nbar, 2, 8, 0)
    ecc_str         = tle_format_field(oe[1], 0, 7, 2)
    inc_str         = tle_format_field(oe[2], 3, 4, 0)
    aop_str         = tle_format_field(oe[3], 3, 4, 0)
    raan_str        = tle_format_field(oe[4], 3, 4, 0)
    ma_str          = tle_format_field(MA, 3, 4, 0)

    # Add other content for Lines 1 & 2 using random/fixed generation
    sat_title_str   = str("STARLINK-") + str(1000 + iter_num)
    sat_cat_str     = str(int(np.random.uniform(10000, 99999)))
    sat_class_str   = "U"
    intl_desg_str   = "98065A  "
    d1nbar_str      = "-.00022620"
    d2nbar_str      = " 00000-0"
    bstr_str        = " 13156-2"
    eph_type        = "0"
    elem_num_str    = " 100"
    chksum1_str     = str(int(np.random.uniform(0, 9)))
    chksum2_str     = str(int(np.random.uniform(0, 9)))
    rev_num_str     = str(int(np.random.uniform(10000, 99999)))

    # Consolidate into a single TLE set
    TLE_output_L0   = sat_title_str
    TLE_output_L1   = "1 " + (sat_cat_str+sat_class_str+" ") + (intl_desg_str+" ") + (epoch_str+" ") + (d1nbar_str+" ") + (d2nbar_str+" ") + (bstr_str+" ") + (eph_type+" ") + (elem_num_str+chksum1_str)
    TLE_output_L2   = "2 " + (sat_cat_str+" ") + (inc_str+" ") + (raan_str+" ") + (ecc_str+" ") + (aop_str+" ") + (ma_str+" ") + (nbar_str+rev_num_str+chksum2_str)
    TLE_output      = TLE_output_L0 + "\n" + TLE_output_L1 + "\n" + TLE_output_L2


    return TLE_output

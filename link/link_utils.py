'''
SpaceNet Link Utils
AUTHOR:         Mohamed M. Kassem, Ph.D.
                University of Surrey
EDITOR:         Bruce Barbour
                Virginia Tech
DESCRIPTION:    This Python script supplies the primary utility functions for computing link characteristics.
CONTENTS:       LINK UTILITY FUNCTIONS
                    get_weather_info(lat, lon, init_timestamp)
                    calc_gsl_snr(satellite, ground_station, t, distance, direction)
                    calc_gsl_snr_given_distance(gsl_distance)
                MAPPING FUNCTIONS
                    map_environment(environment_str)
                    map_building_type(building_type_str)
                    map_path_condition(condition_str)
                    map_link_direction(direction_str)
'''

# =================================================================================== #
# ------------------------------- IMPORT PACKAGES ----------------------------------- #
# =================================================================================== #

import math
import itur
import requests
import sys
import logging
sys.path.append("../")

# Import NTN Channel Model and related classes
from link.ntn_channel_model import NTNChannelModel, LinkBudget
from link.enum1 import Environment, BuildingType, PathCondition, LinkDirection

# Set up logging
logger = logging.getLogger(__name__)

# =================================================================================== #
# ---------------------------- BUILT-IN ASSUMPTIONS --------------------------------- #
# =================================================================================== #
#api_key                                 = "d06b0a02f8377dff811a2a6d0882a2d6"
api_key                                 = "cab2710f043a0aeedb61b28b3a316146" #Contains weather history subscription of OpenWeather API (One Call API 3.0)
channelFreq_isls                        = 37.0      # GHz
channelFreq_sat_to_ground               = 12.7      # GHz
channelFreq_ground_to_sat               = 14.5      # GHz
channnel_bandwidth_downlink             = 240       # MHz
channnel_bandwidth_uplink               = 60        # MHz
polarization_loss                       = 3         # dBi
misalignment_attenuation_losses         = 0.5       # dB
starlink_merit_figure                   = 9.2       # dB/K
satellite_eirp                          = 80.9
satellite_eirp_dbW                      = 50.9
ground_station_tx_power                 = 36.08526  # dBm -- https://apps.fcc.gov/els/GetAtt.html?id=259301
ground_station_receive_attenna_gain     = 33.2      # dBi -- https://apps.fcc.gov/els/GetAtt.html?id=259301
ground_station_transmit_attenna_gain    = 34.6      # dBi -- https://apps.fcc.gov/els/GetAtt.html?id=259301

# =================================================================================== #
# ---------------------------- MAPPING FUNCTIONS ------------------------------------ #
# =================================================================================== #

def map_environment(environment_str):
    """
    Convert environment string to Environment enum

    Args:
        environment_str (str): String representation of environment

    Returns:
        Environment: The corresponding Environment enum value
    """
    if environment_str is None:
        logger.warning("Environment string is None, defaulting to SUBURBAN_RURAL")
        return Environment.SUBURBAN_RURAL

    env_map = {
        "DENSE_URBAN": Environment.DENSE_URBAN,
        "URBAN": Environment.URBAN,
        "SUBURBAN_RURAL": Environment.SUBURBAN_RURAL
    }

    try:
        return env_map[environment_str.upper()]
    except KeyError:
        logger.warning(f"Unknown environment: {environment_str}, defaulting to SUBURBAN_RURAL")
        return Environment.SUBURBAN_RURAL


def map_building_type(building_type_str):
    """
    Convert building type string to BuildingType enum

    Args:
        building_type_str (str): String representation of building type

    Returns:
        BuildingType: The corresponding BuildingType enum value or None if indoor=False
    """
    if building_type_str is None:
        return None

    bldg_map = {
        "TRADITIONAL": BuildingType.TRADITIONAL,
        "THERMALLY_EFFICIENT": BuildingType.THERMALLY_EFFICIENT
    }

    try:
        return bldg_map[building_type_str.upper()]
    except KeyError:
        logger.warning(f"Unknown building type: {building_type_str}, defaulting to None")
        return None


def map_path_condition(condition_str):
    """
    Convert path condition string to PathCondition enum

    Args:
        condition_str (str): String representation of path condition

    Returns:
        PathCondition: The corresponding PathCondition enum value or None for probabilistic
    """
    if condition_str is None:
        return None

    cond_map = {
        "LOS": PathCondition.LOS,
        "NLOS": PathCondition.NLOS
    }

    try:
        return cond_map[condition_str.upper()]
    except KeyError:
        logger.warning(f"Unknown path condition: {condition_str}, defaulting to None (probabilistic)")
        return None


def map_link_direction(direction_str):
    """
    Convert link direction string to LinkDirection enum

    Args:
        direction_str (str): String representation of link direction

    Returns:
        LinkDirection: The corresponding LinkDirection enum value
    """
    if direction_str is None:
        logger.warning("Link direction string is None, defaulting to DOWNLINK")
        return LinkDirection.DOWNLINK

    dir_map = {
        "UPLINK": LinkDirection.UPLINK,
        "DOWNLINK": LinkDirection.DOWNLINK
    }

    try:
        return dir_map[direction_str.upper()]
    except KeyError:
        logger.warning(f"Unknown link direction: {direction_str}, defaulting to DOWNLINK")
        return LinkDirection.DOWNLINK

# =================================================================================== #
# ---------------------------- LINK UTILITY FUNCTIONS ------------------------------- #
# =================================================================================== #
def get_weather_info(
                        lat            : float,
                        lon            : float,
                        init_timestamp : int,
                        wait_on_rate_limit=False
                    ) -> dict:
    """
    Retrieve weather information using OpenWeatherMap API based on latitude and longitude.

    Args:
        lat (float)         :    Latitude of the location
        lon (float)         :    Longitude of the location
        init_timestamp (int):    Initial simulation timestamp in unix time
        wait_on_rate_limit (bool): Whether to wait if rate limit is exceeded

    Returns:
        dict:           Dictionary containing temperature, humidity, pressure, and weather description
    """

    # Construct the OpenWeatherMap API URL using latitude, longitude, and API key
    url = "https://api.openweathermap.org/data/3.0/onecall/timemachine?lat=%s&lon=%s&dt=%s&appid=%s" % (str(lat), str(lon), str(init_timestamp), api_key)

    # Send a GET request to the API
    try:
        response = requests.get(url)
    except requests.exceptions.RequestException as e:
        logger.warning(f"Error: Unable to connect to the weather API - skipping... {e}")
        return ""

    if response.status_code == 429: # Too many requests
        if not wait_on_rate_limit:
            return ""
        retry_after = int(response.headers.get("Retry-After", 10))
        logger.info(f"Rate limit exceeded. Waiting for {retry_after} seconds...")
        import time
        time.sleep(retry_after)
        return get_weather_info(lat, lon, init_timestamp, wait_on_rate_limit=True)

    # Parse the JSON response
    try:
        data = response.json()
    except ValueError:
        logger.warning("Error: Invalid JSON response from weather API")
        return ""

    # Initialize variables to store weather information
    temp = None
    humidity = None
    pressure = None
    description = None

    # Check if the response contains weather information
    if data != "":
        try:
            # Extract weather description
            w_data = data['data'][0]
            da = w_data["weather"]
            description = da[0]["description"]

            # Extract general weather data
            temp = w_data["temp"]
            humidity = w_data["humidity"]
            pressure = w_data["pressure"]
        except KeyError as e:
            logger.warning(f"Error: Unable to extract weather data - skipping... {e}")
            return ""

    # Return a dictionary containing weather information
    return {
        "temp": temp,
        "humidity": humidity,
        "pressure": pressure,
        "description": description
    }


def calc_gsl_snr(
                    satellite           : dict,
                    ground_station      : dict,
                    t                   : float,
                    distance            : float,
                    direction           : str,
                    use_ntn_model       : bool = False,
                    satellite_altitude  : float = None,
                    frequency_ghz       : float = None,
                    environment         : str = "SUBURBAN_RURAL",
                    is_indoor           : bool = False,
                    building_type       : str = None,
                    bandwidth_hz        : float = None,
                    add_interference    : bool = False,
                    forced_condition    : str = None
                ) -> float:
    """
    Calculate Signal-to-Noise Ratio (SNR) for a ground station and satellite communication link

    Args:
        satellite (dict):       Dictionary containing satellite parameters
        ground_station (dict):  Dictionary containing ground station parameters
        t (float):              Time of communication
        distance (float):       Distance between ground station and satellite in meters
        direction (str):        Communication direction, "downlink" or "uplink"
        use_ntn_model (bool):   Whether to use the NTN Channel Model for calculation
        satellite_altitude (float): Altitude of the satellite in meters (only used with NTN model)
        frequency_ghz (float):   Carrier frequency in GHz (only used with NTN model)
        environment (str):       UE environment type: "DENSE_URBAN", "URBAN", or "SUBURBAN_RURAL"
        is_indoor (bool):        Whether the UE is indoor
        building_type (str):     Building type for indoor UEs: "TRADITIONAL" or "THERMALLY_EFFICIENT"
        bandwidth_hz (float):    Signal bandwidth in Hz
        add_interference (bool): Whether to add interference in the calculation
        forced_condition (str):  Force a specific path condition: "LOS" or "NLOS"

    Returns:
        float:                  SNR in dB
    """
    # Initialize variables
    gsl_distance = distance

    # Get ground station latitude and longitude
    lat_gs = float(ground_station["latitude_degrees_str"])
    lon_gs = float(ground_station["longitude_degrees_str"])

    # Use traditional method if NTN model is not requested
    if not use_ntn_model:
        # Calculate SNR using the traditional method

        # Frequency for downlink and uplink
        f_dl = channelFreq_sat_to_ground * itur.u.GHz
        f_ul = channelFreq_ground_to_sat * itur.u.GHz

        # Free space path loss calculation
        if direction.lower() == "downlink":
            fspl = 20 * math.log10(gsl_distance/1000) + 20 * math.log10(channelFreq_sat_to_ground) + 92.45
            freq = f_dl
            bandwidth = channnel_bandwidth_downlink
        else:  # uplink
            fspl = 20 * math.log10(gsl_distance/1000) + 20 * math.log10(channelFreq_ground_to_sat) + 92.45
            freq = f_ul
            bandwidth = channnel_bandwidth_uplink

        # Get weather information for the ground station
        if 'weather_data' in ground_station:  # data already queried
            weather_data = ground_station['weather_data']
        else:
            # Need init_timestamp from satellite or t
            init_timestamp = int(t) if isinstance(t, (int, float)) else int(satellite.get("init_timestamp", t))
            weather_data = get_weather_info(lat_gs, lon_gs, init_timestamp)

        # Weather attenuation calculation
        weather_attenuation = 0
        if weather_data and weather_data != "":
            # Antenna size and elevation angle
            D = 0.58 * itur.u.m   # Size of the receiver antenna (starlink dish v.1 diameter)
            el = 70                # Elevation angle constant of 70 degrees

            # Percentage of time that attenuation values are exceeded
            p = 0.01

            # Set rain rate based on weather description
            if "drizzle" in str(weather_data["description"]):
                r001 = 0.25
            elif "light rain" in str(weather_data["description"]):
                r001 = 2.5
            elif "moderate rain" in str(weather_data["description"]):
                r001 = 12.5
            elif str(weather_data["description"]) == "heavy rain":
                r001 = 25
            elif str(weather_data["description"]) == "very heavy rain" or str(weather_data["description"]) == "extreme rain":
                r001 = 50
            elif str(weather_data["description"]) == "heavy intensity shower rain" or str(weather_data["description"]) == "shower rain":
                r001 = 100
            elif str(weather_data["description"]) == "ragged shower rain":
                r001 = 150
            else:
                r001 = None

            # Atmospheric attenuation calculation based on weather conditions
            temp = float(weather_data["temp"])
            humidity = float(weather_data["humidity"])
            pressure = float(weather_data["pressure"])

            try:
                # Calculate attenuation with weather data
                weather_attenuation = itur.atmospheric_attenuation_slant_path(
                    lat_gs, lon_gs, freq, el, p, D, R001=r001, T=temp, H=humidity, P=pressure
                )
                weather_attenuation = weather_attenuation.value
            except Exception as e:
                logger.warning(f"Error calculating weather attenuation: {e}")
                # Calculate attenuation without weather data
                weather_attenuation = itur.atmospheric_attenuation_slant_path(
                    lat_gs, lon_gs, freq, el, p, D, return_contributions=True
                )
                if isinstance(weather_attenuation, tuple):
                    weather_attenuation = weather_attenuation[4]  # 4th index is the total attenuation
                weather_attenuation = weather_attenuation.value
        else:
            # Calculate attenuation without weather data
            # Antenna size and elevation angle
            D = 0.58 * itur.u.m   # Size of the receiver antenna (starlink dish v.1 diameter)
            el = 70                # Elevation angle constant of 70 degrees

            # Percentage of time that attenuation values are exceeded
            p = 0.01

            try:
                weather_attenuation = itur.atmospheric_attenuation_slant_path(
                    lat_gs, lon_gs, freq, el, p, D, return_contributions=True
                )
                if isinstance(weather_attenuation, tuple):
                    weather_attenuation = weather_attenuation[4]  # 4th index is the total attenuation
                weather_attenuation = weather_attenuation.value
            except Exception as e:
                logger.warning(f"Error calculating atmospheric attenuation: {e}")
                weather_attenuation = 0

        # Calculate elevation angle from satellite (simplified model)
        satellite_alt_km = satellite_altitude / 1000 if satellite_altitude is not None else 550  # Default to 550 km if not provided
        earth_radius_km = 6371  # Earth radius in km
        slant_range_km = gsl_distance / 1000  # Convert to km

        # Approximate elevation angle calculation
        try:
            elevation_angle = math.degrees(math.atan2(satellite_alt_km,
                                                     math.sqrt(max(0, slant_range_km**2 - satellite_alt_km**2))))
        except (ValueError, ZeroDivisionError):
            # Fallback to a default angle if calculation fails
            elevation_angle = 45
            logger.warning(f"Elevation angle calculation failed, using default of {elevation_angle} degrees")

        # Convert environment string to enum or use default
        environment_enum = Environment.SUBURBAN_RURAL
        building_type_enum = None

        # Create NTN channel model
        channel_model = NTNChannelModel(
            satellite_altitude=satellite_altitude if satellite_altitude is not None else 550000,  # Default to 550 km
            frequency_ghz=channelFreq_sat_to_ground if direction.lower() == "downlink" else channelFreq_ground_to_sat,
            environment=environment_enum,
            is_indoor=False,  # Assuming ground stations are not indoor
            building_type=building_type_enum
        )

        # Create link budget calculator based on direction
        if direction.lower() == "downlink":
            # Downlink parameters (satellite to ground station)
            tx_power_dbm = 43  # Satellite transmit power (20 Watts = 43 dBm)
            tx_gain_dbi = 30  # Satellite antenna gain
            rx_gain_dbi = ground_station_receive_attenna_gain  # Ground station receive antenna gain
            noise_figure_db = 7  # Ground station noise figure

            # Create link budget calculator for downlink
            link_budget = LinkBudget(
                channel_model=channel_model,
                tx_power_dbm=tx_power_dbm,
                tx_gain_dbi=tx_gain_dbi,
                rx_gain_dbi=rx_gain_dbi,
                noise_figure_db=noise_figure_db,
                bandwidth_hz=channnel_bandwidth_downlink * 1e6,  # Convert MHz to Hz
                link_direction=LinkDirection.DOWNLINK,
                add_interference=False
            )

        else:  # uplink
            # Uplink parameters (ground station to satellite)
            tx_power_dbm = ground_station_tx_power  # Ground station transmit power
            tx_gain_dbi = ground_station_transmit_attenna_gain  # Ground station transmit antenna gain
            rx_gain_dbi = 30  # Satellite receive antenna gain
            noise_figure_db = 2  # Satellite noise figure

            # Create link budget calculator for uplink
            link_budget = LinkBudget(
                channel_model=channel_model,
                tx_power_dbm=tx_power_dbm,
                tx_gain_dbi=tx_gain_dbi,
                rx_gain_dbi=rx_gain_dbi,
                noise_figure_db=noise_figure_db,
                bandwidth_hz=channnel_bandwidth_uplink * 1e6,  # Convert MHz to Hz
                link_direction=LinkDirection.UPLINK,
                add_interference=False
            )

        # Calculate SNR using the NTN channel model
        snr, path_loss_results = link_budget.calculate_snr(elevation_angle, None)

        # Apply weather attenuation if available (as additional loss)
        if weather_attenuation > 0:
            snr -= weather_attenuation

        return snr

    else:
        # Use the NTN Channel Model
        if satellite_altitude is None:
            raise ValueError("satellite_altitude must be provided when using NTN model")
        if frequency_ghz is None:
            # Use default frequencies if not provided
            frequency_ghz = channelFreq_sat_to_ground if direction.lower() == "downlink" else channelFreq_ground_to_sat

        # Set bandwidth if not provided
        if bandwidth_hz is None:
            bandwidth_hz = channnel_bandwidth_downlink * pow(10, 6) if direction.lower() == "downlink" else channnel_bandwidth_uplink * pow(10, 6)

        # Convert string parameters to enums
        environment_enum = map_environment(environment)
        building_type_enum = map_building_type(building_type) if is_indoor else None
        forced_condition_enum = map_path_condition(forced_condition)

        # Calculate elevation angle from satellite
        # This should be calculated based on satellite position relative to ground station
        # For demonstration, we'll use a simplistic calculation
        # In a real implementation, this would use a more accurate model
        satellite_alt_km = satellite_altitude / 1000  # Convert to km
        earth_radius_km = 6371  # Earth radius in km

        # Approximate elevation angle calculation
        # This is a simplified calculation assuming flat Earth and satellite directly above
        # A more accurate calculation would use the actual satellite position
        slant_range_km = gsl_distance / 1000  # Convert to km
        elevation_angle = math.degrees(math.atan2(satellite_alt_km,
                                                 math.sqrt(slant_range_km**2 - satellite_alt_km**2)))

        # Create NTN channel model
        channel_model = NTNChannelModel(
            satellite_altitude=satellite_altitude,
            frequency_ghz=frequency_ghz,
            environment=environment_enum,
            is_indoor=is_indoor,
            building_type=building_type_enum
        )

        # Create link budget calculator based on direction
        if direction.lower() == "downlink":
            # Downlink parameters (satellite to UE)
            tx_power_dbm = 43  # Satellite transmit power (20 Watts = 43 dBm)
            tx_gain_dbi = 30  # Satellite antenna gain
            rx_gain_dbi = ground_station_receive_attenna_gain  # Ground station antenna gain
            noise_figure_db = 7  # Ground station noise figure

            # Create link budget calculator
            link_budget = LinkBudget(
                channel_model=channel_model,
                tx_power_dbm=tx_power_dbm,
                tx_gain_dbi=tx_gain_dbi,
                rx_gain_dbi=rx_gain_dbi,
                noise_figure_db=noise_figure_db,
                bandwidth_hz=bandwidth_hz,
                link_direction=LinkDirection.DOWNLINK,
                add_interference=add_interference
            )

        elif direction.lower() == "uplink":
            # Uplink parameters (UE to satellite)
            tx_power_dbm = ground_station_tx_power  # Ground station transmit power
            tx_gain_dbi = ground_station_transmit_attenna_gain  # Ground station antenna gain
            rx_gain_dbi = 30  # Satellite antenna gain
            noise_figure_db = 2  # Satellite noise figure

            # Create link budget calculator
            link_budget = LinkBudget(
                channel_model=channel_model,
                tx_power_dbm=tx_power_dbm,
                tx_gain_dbi=tx_gain_dbi,
                rx_gain_dbi=rx_gain_dbi,
                noise_figure_db=noise_figure_db,
                bandwidth_hz=bandwidth_hz,
                link_direction=LinkDirection.UPLINK,
                add_interference=add_interference
            )

        else:
            raise ValueError("Invalid link direction. Must be 'uplink' or 'downlink'.")

        # Calculate SNR using the NTN channel model
        snr, _ = link_budget.calculate_snr(elevation_angle, forced_condition_enum)

        return snr


def calc_gsl_snr_given_distance(
                                    gsl_distance    : float
                               ) -> float:
    """
    Calculate Signal-to-Noise Ratio (SNR) for a ground station and satellite communication link given the distance between them.
    Does not use the Weather API.

    Args:
        gsl_distance (float):   Distance between ground station and satellite in meters

    Returns:
        float:                  SNR in dB
    """
    # Free space path loss calculation
    fspl = 20 * math.log10(gsl_distance/1000) + 20 * math.log10(channelFreq_sat_to_ground) + 92.45

    # Receive Signal Strength (RSS) calculation
    rss_dBm = satellite_eirp - 2 + ground_station_receive_attenna_gain - fspl - polarization_loss - misalignment_attenuation_losses - 1.0
    rss_watt = pow(10, ((rss_dBm - 30) / 10))

    # Noise power calculation
    noise_watt = 290 * 1.38064852 * pow(10, -23) * channnel_bandwidth_downlink * pow(10, 6)   # ktB

    # Calculate Signal-to-Noise Ratio (SNR)
    snr_linear = rss_watt / noise_watt
    snr_db = 10 * math.log10(snr_linear)

    # Return the calculated SNR
    return snr_db
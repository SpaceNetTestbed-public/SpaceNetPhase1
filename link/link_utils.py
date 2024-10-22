''''

SpaceNet Link Utils/

AUTHOR:         Mohamed M. Kassem, Ph.D.
                University of Surrey

EDITOR:         Bruce Barbour
                Virginia Tech

DESCRIPTION:    This Python script supplies the primary utility functions for computing link characteristics.


CONTENTS:       LINK UTILITY FUNCTIONS/ (STARTS AT 55)
                    get_weather_info(lat, lon)
                    calc_gsl_snr(satellite, ground_station, t, distance, direction)
                    calc_gsl_snr_given_distance(gsl_distance)
'''

# =================================================================================== #
# ------------------------------- IMPORT PACKAGES ----------------------------------- #
# =================================================================================== #

import math
import itur
import requests
import sys
sys.path.append("../")
from mobility.mobility_utils import *

# =================================================================================== #
# ---------------------------- BUILT-IN ASSUMPTIONS --------------------------------- #
# =================================================================================== #
api_key                                 = "d06b0a02f8377dff811a2a6d0882a2d6"
channelFreq_isls                        = 37.0      # GHz
channelFreq_sat_to_ground               = 12.7      # GHz
channelFreq_ground_to_sat               = 14.5      # GH
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
# ---------------------------- LINK UTILITY FUNCTIONS ------------------------------- #
# =================================================================================== #
def get_weather_info(
                        lat : float, 
                        lon : float,
                        wait_on_rate_limit=False
                    ) -> dict:
    """
    Retrieve weather information using OpenWeatherMap API based on latitude and longitude.

    Args:
        lat (float):    Latitude of the location
        lon (float):    Longitude of the location

    Returns:
        dict:           Dictionary containing temperature, humidity, pressure, and weather description
    """

    # Construct the OpenWeatherMap API URL using latitude, longitude, and API key
    url = "https://api.openweathermap.org/data/2.5/weather?lat=%s&lon=%s&appid=%s&units=standard" % (str(lat), str(lon), api_key)
    
    # Send a GET request to the API
    try:
        response = requests.get(url)
    except requests.exceptions.RequestException as e:
        print("Error: Unable to connect to the weather API - skipping...")
        return ""
    if response.status_code == 429: # Too many requests
        if not wait_on_rate_limit:
            return ""
        retry_after = int(response.headers.get("Retry-After", 10))
        print(f"Rate limit exceeded. Waiting for {retry_after} seconds...")
        import time
        time.sleep(retry_after)
        return get_weather_info(lat, lon, wait_on_rate_limit=True)
    # Parse the JSON response
    data = response.json()

    # Initialize variables to store weather information
    temp = None
    humidity = None
    pressure = None
    description = None

    # Check if the response contains weather information
    if data != "":

        try:
            # Extract weather description
            da = data["weather"]
            description =  da[0]["description"]

            # Extract general weather data
            general = data["main"]
            temp  = general["temp"]
            humidity = general["humidity"]
            pressure = general["pressure"]
        except KeyError:
            #print("Error: Unable to extract weather data - skipping...")
            return ""


    # Return a dictionary containing weather information
    return {"temp": temp,
            "humidity": humidity,
            "pressure": pressure,
            "description": description
            }


def calc_gsl_snr(
                    satellite           : dict, 
                    ground_station      : dict, 
                    t                   : float, 
                    distance            : float, 
                    direction           : str
                ) -> float:
    """
    Calculate Signal-to-Noise Ratio (SNR) for a ground station and satellite communication link

    Args:
        satellite (dict):       Dictionary containing satellite parameters (unused)
        ground_station (dict):  Dictionary containing ground station parameters
        t (float):              Time of communication (unused)
        distance (float):       Distance between ground station and satellite (if provided, then ignore 'satellite' and 't' arguments)
        direction (str):        Communication direction, "downlink" or "uplink"

    Returns:
        float:                  SNR
    """

    # Initialize variable
    gsl_distance = distance

    # *Computes distance already in the mobility_utils
    # distance_between_ground_station_satellite(ground_station, satellite, t)

    # Free space path loss calculation
    fspl = 20 * math.log10(gsl_distance/1000) + 20 * math.log10(channelFreq_sat_to_ground) + 92.45

    # Ground station latitude and longitude
    lat_gs  = float(ground_station["latitude_degrees_str"])
    lon_gs  = float(ground_station["longitude_degrees_str"])

    # Frequency for downlink and uplink
    f_dl    = channelFreq_sat_to_ground * itur.u.GHz 
    f_ul    = channelFreq_ground_to_sat * itur.u.GHz

    # Antenna size and elevation angle
    D       = 0.58 * itur.u.m   # Size of the receiver antenna (this is the diamter of starlink dish v.1)
    el      = 70                # Elevation angle constant of 60 degrees
    
    # Percentage of time that attenuation values are exceeded
    p       = 0.01

    # Get weather information for the ground station
    if 'weather_data' in ground_station: # data already queried
        weather_data = ground_station['weather_data']
    else:
        weather_data = get_weather_info(lat_gs, lon_gs)

    # Check if weather data isn't empty
    if weather_data != "":

        # Drizzle 
        if "drizzle" in str(weather_data["description"]):
            r001 = 0.25

        # Light rain
        elif "light rain" in str(weather_data["description"]):
            r001 = 2.5

        # Moderate rain
        elif "moderate rain" in str(weather_data["description"]):
            r001 = 12.5

        # Heavy rain
        elif str(weather_data["description"]) == "heavy rain":
            r001 = 25

        # Very heavy rain
        elif str(weather_data["description"]) == "very heavy rain" or str(weather_data["description"]) == "extreme rain":
            r001 = 50

        # Extreme rain
        elif str(weather_data["description"]) == "heavy intensity shower rain" or str(weather_data["description"]) == "shower rain":
            r001 = 100

        # Ragged shower rain
        elif str(weather_data["description"]) == "ragged shower rain":
            r001 = 150

        # Clear conditions
        else:
            r001 = None

        # Atmospheric attenuation calculation based on weather conditions
        temp = float(weather_data["temp"])	#temp in Kelvin
        humidity = float(weather_data["humidity"])
        pressure = float(weather_data["pressure"])
        # temp = temp + 273.15
        # print(temp, humidity, pressure)

        # Downlink attenuation with weather 
        weather_attenuation_dl = itur.atmospheric_attenuation_slant_path(lat_gs, lon_gs, f_dl, el, p, D, R001=r001, T=temp, H=humidity, P=pressure)
        weather_attenuation_dl = weather_attenuation_dl.value

        # Uplink attenuation with weather
        weather_attenuation_ul = itur.atmospheric_attenuation_slant_path(lat_gs, lon_gs, f_ul, el, p, D, R001=r001, T=temp, H=humidity, P=pressure)
        weather_attenuation_ul = weather_attenuation_ul.value
        # print(weather_attenuation)

    else:

        #print("no weather data -- ")

        # Downlink attenuation without weather 
        weather_attenuation_dl = itur.atmospheric_attenuation_slant_path(lat_gs, lon_gs, f_dl, el, p, D, return_contributions=True)
        if type(weather_attenuation_dl) == tuple:
            weather_attenuation_dl = weather_attenuation_dl[4] # 4th index is the total attenuation
        weather_attenuation_dl = weather_attenuation_dl.value

        # Uplink attenuation without weather
        weather_attenuation_ul = itur.atmospheric_attenuation_slant_path(lat_gs, lon_gs, f_ul, el, p, D, return_contributions=True)
        if type(weather_attenuation_ul) == tuple:
            weather_attenuation_ul = weather_attenuation_ul[4] # 4th index is the total attenuation
        weather_attenuation_ul = weather_attenuation_ul.value


    # Calculate SNR for downlink (sat -> gs)
    if direction == "downlink":
        snr_db = satellite_eirp_dbW - (10*(math.log(channnel_bandwidth_downlink*pow(10,6))/math.log(10))) - fspl - polarization_loss - misalignment_attenuation_losses - weather_attenuation_dl - 3 + starlink_merit_figure+228.6
        snr = pow(10,(snr_db/10))
        # rss_dBm = satellite_eirp - 2 + ground_station_receive_attenna_gain - fspl - polarization_loss - misalignment_attenuation_losses - weather_attenuation - 1.0;
        # rss_watt = pow(10,((rss_dBm - 30)/10));
        # print(fspl, snr_db, weather_attenuation)
        return snr

    # Calculate SNR for uplink (gs -> sat)
    if direction == "uplink":
        rss_dBm = ground_station_tx_power + ground_station_transmit_attenna_gain + 10 - 2 - fspl - polarization_loss - misalignment_attenuation_losses - weather_attenuation_ul - 1.0
        rss_watt = pow(10,((rss_dBm - 30)/10))
        # snr = pow(10,(snr_db/10))
        # return snr

    # Noise power calculation
    noise_watt = 200 * 1.38064852 * pow(10, -23) * 250*pow(10, 6); # ktB channnel_bandwidth_downlink
    # print(fspl, rss_dBm, rss_watt, noise_watt)

    # Calculate SNR
    snr = rss_watt/noise_watt

    # Return the calculated SNR
    return snr


def calc_gsl_snr_given_distance(
                                    gsl_distance    : float
                               ) -> float:
    """
    Calculate Signal-to-Noise Ratio (SNR) for a ground station and satellite communication link given the distance between them.
    Does not use the Weather API.

    Args:
        gsl_distance (float):   Distance between ground station and satellite

    Returns:
        float:                  SNR
    """

    # Free space path loss calculation
    fspl = 20 * math.log10(gsl_distance/1000) + 20 * math.log10(channelFreq_sat_to_ground) + 92.45

    # Receive Signal Strength (RSS) calculation
    rss_dBm = satellite_eirp - 2 + ground_station_receive_attenna_gain - fspl - polarization_loss - misalignment_attenuation_losses - 1.0
    rss_watt = pow(10,((rss_dBm - 30)/10))

    # Noise power calculation
    noise_watt = 200 * 1.38064852 * pow(10, -23) * channnel_bandwidth_downlink*pow(10, 6)   # ktB

    # Calculate Signal-to-Noise Ratio (SNR)
    snr = rss_watt/noise_watt

    # Return the calculated SNR
    return snr
from skyfield.api import load, wgs84
from datetime import datetime


def calculate_satellite_elevation(tle_lines, ground_station_coords, time=None):
    """
    Calculate satellite elevation angle using Skyfield

    This function calculates the actual elevation angle and altitude of a satellite
    from a ground station using TLE data.

    Args:
        tle_lines (list): TLE data (2 lines)
        ground_station_coords (tuple): (latitude, longitude, altitude in km)
        time (datetime): Time for calculation (default: current time)

    Returns:
        tuple: (elevation angle in degrees, satellite altitude in meters)
    """
    # Load time scale
    ts = load.timescale()

    # Parse TLE lines to create a satellite
    from skyfield.sgp4lib import EarthSatellite

    line1 = tle_lines[0]
    line2 = tle_lines[1]
    satellite = EarthSatellite(line1, line2, 'SAT', ts)

    # Define the ground station
    ground_lat, ground_lon, ground_alt = ground_station_coords
    ground_station = wgs84.latlon(ground_lat, ground_lon, ground_alt * 1000)  # Skyfield uses meters

    # Set time for calculation
    if time is None:
        t = ts.now()
    else:
        t = ts.from_datetime(time)

    # Calculate the position of the satellite relative to the ground station
    difference = satellite - ground_station
    topocentric = difference.at(t)

    # Extract altitude, azimuth and distance
    alt, az, distance = topocentric.altaz()

    # Get elevation angle in degrees
    elevation_angle = alt.degrees

    # Get the satellite's distance from Earth's center
    satellite_altitude = satellite.at(t).distance().km * 1000  # in meters

    return elevation_angle, satellite_altitude


import numpy as np
from scipy.stats import norm
from link.enum1 import Environment, FrequencyBand, PathCondition, BuildingType, LinkDirection
import math

class NTNChannelModel:
    """NTN Channel Model based on 3GPP 38.811 standard"""

    # Earth radius in meters
    EARTH_RADIUS = 6371e3

    # LOS probability tables
    LOS_PROBABILITY = {
        Environment.DENSE_URBAN: {
            10: 0.282, 20: 0.331, 30: 0.398, 40: 0.468,
            50: 0.537, 60: 0.612, 70: 0.738, 80: 0.820, 90: 0.981
        },
        Environment.URBAN: {
            10: 0.246, 20: 0.386, 30: 0.493, 40: 0.613,
            50: 0.726, 60: 0.805, 70: 0.919, 80: 0.968, 90: 0.992
        },
        Environment.SUBURBAN_RURAL: {
            10: 0.782, 20: 0.869, 30: 0.919, 40: 0.929,
            50: 0.935, 60: 0.940, 70: 0.949, 80: 0.952, 90: 0.998
        }
    }

    # Shadow fading and clutter loss table for LOS and NLOS conditions
    # Format: {environment: {band: {condition: {elevation: (shadow_fading, clutter_loss)}}}}
    SF_CL_TABLES = {
        Environment.DENSE_URBAN: {
            FrequencyBand.S_BAND: {
                PathCondition.LOS: {
                    10: (3.5, 0), 20: (3.4, 0), 30: (2.9, 0), 40: (3.0, 0),
                    50: (3.1, 0), 60: (2.7, 0), 70: (2.5, 0), 80: (2.3, 0), 90: (1.2, 0)
                },
                PathCondition.NLOS: {
                    10: (15.5, 34.3), 20: (13.9, 30.9), 30: (12.4, 29.0), 40: (11.7, 27.7),
                    50: (10.6, 26.8), 60: (10.5, 26.2), 70: (10.1, 25.8), 80: (9.2, 25.5), 90: (9.2, 25.5)
                }
            },
            FrequencyBand.KA_BAND: {
                PathCondition.LOS: {
                    10: (2.9, 0), 20: (2.4, 0), 30: (2.7, 0), 40: (2.4, 0),
                    50: (2.4, 0), 60: (2.7, 0), 70: (2.6, 0), 80: (2.8, 0), 90: (0.6, 0)
                },
                PathCondition.NLOS: {
                    10: (17.1, 44.3), 20: (17.1, 39.9), 30: (15.6, 37.5), 40: (14.6, 35.8),
                    50: (14.2, 34.6), 60: (12.6, 33.8), 70: (12.1, 33.3), 80: (12.3, 33.0), 90: (12.3, 32.9)
                }
            }
        },
        Environment.URBAN: {
            FrequencyBand.S_BAND: {
                PathCondition.LOS: {
                    10: (4, 0), 20: (4, 0), 30: (4, 0), 40: (4, 0),
                    50: (4, 0), 60: (4, 0), 70: (4, 0), 80: (4, 0), 90: (4, 0)
                },
                PathCondition.NLOS: {
                    10: (6, 34.3), 20: (6, 30.9), 30: (6, 29.0), 40: (6, 27.7),
                    50: (6, 26.8), 60: (6, 26.2), 70: (6, 25.8), 80: (6, 25.5), 90: (6, 25.5)
                }
            },
            FrequencyBand.KA_BAND: {
                PathCondition.LOS: {
                    10: (4, 0), 20: (4, 0), 30: (4, 0), 40: (4, 0),
                    50: (4, 0), 60: (4, 0), 70: (4, 0), 80: (4, 0), 90: (4, 0)
                },
                PathCondition.NLOS: {
                    10: (6, 44.3), 20: (6, 39.9), 30: (6, 37.5), 40: (6, 35.8),
                    50: (6, 34.6), 60: (6, 33.8), 70: (6, 33.3), 80: (6, 33.0), 90: (6, 32.9)
                }
            }
        },
        Environment.SUBURBAN_RURAL: {
            FrequencyBand.S_BAND: {
                PathCondition.LOS: {
                    10: (1.79, 0), 20: (1.14, 0), 30: (1.14, 0), 40: (0.92, 0),
                    50: (1.42, 0), 60: (1.56, 0), 70: (0.85, 0), 80: (0.72, 0), 90: (0.72, 0)
                },
                PathCondition.NLOS: {
                    10: (8.93, 19.52), 20: (9.08, 18.17), 30: (8.78, 18.42), 40: (10.25, 18.28),
                    50: (10.56, 18.63), 60: (10.74, 17.68), 70: (10.17, 16.50), 80: (11.52, 16.30), 90: (11.52, 16.30)
                }
            },
            FrequencyBand.KA_BAND: {
                PathCondition.LOS: {
                    10: (1.9, 0), 20: (1.6, 0), 30: (1.9, 0), 40: (2.3, 0),
                    50: (2.7, 0), 60: (3.1, 0), 70: (3.0, 0), 80: (3.6, 0), 90: (0.4, 0)
                },
                PathCondition.NLOS: {
                    10: (10.7, 29.5), 20: (10.0, 24.6), 30: (11.2, 21.9), 40: (11.6, 20.0),
                    50: (11.8, 18.7), 60: (10.8, 17.8), 70: (10.8, 17.2), 80: (10.8, 16.9), 90: (10.8, 16.8)
                }
            }
        }
    }

    # Coefficients for O2I penetration loss
    O2I_COEFFICIENTS = {
        BuildingType.TRADITIONAL: {
            'r': 12.64, 's': 3.72, 't': 0.96, 'u': 9.6, 'v': 2.0,
            'w': 9.1, 'x': -3.0, 'y': 4.5, 'z': -2.0
        },
        BuildingType.THERMALLY_EFFICIENT: {
            'r': 28.19, 's': -3.00, 't': 8.48, 'u': 13.5, 'v': 3.8,
            'w': 27.8, 'x': -2.9, 'y': 9.4, 'z': -2.1
        }
    }

    # Atmospheric gaseous attenuation at zenith (from ITU-R P.676)
    # Approximated values for sea level and standard atmosphere
    ATMOSPHERIC_ATTENUATION_ZENITH = {
        1: 0.0075,  # 1 GHz
        2: 0.012,  # 2 GHz (S-band)
        6: 0.03,  # 6 GHz
        10: 0.1,  # 10 GHz
        20: 0.4,  # 20 GHz (Ka-band)
        30: 0.22  # 30 GHz (Ka-band)
    }

    # Tropospheric scintillation attenuation (dB) with 99% probability at 20 GHz
    TROPOSPHERIC_SCINTILLATION = {
        10: 1.08, 20: 0.48, 30: 0.30, 40: 0.22,
        50: 0.17, 60: 0.13, 70: 0.12, 80: 0.12, 90: 0.12
    }

    def __init__(self, satellite_altitude, frequency_ghz, environment, is_indoor=False, building_type=None):
        """
        Initialize the NTN channel model

        Args:
            satellite_altitude (float): Altitude of the satellite in meters
            frequency_ghz (float): Carrier frequency in GHz
            environment (Environment): UE environment type
            is_indoor (bool): Whether the UE is indoor
            building_type (BuildingType): Building type for indoor UEs
        """
        self.satellite_altitude = satellite_altitude
        self.frequency_ghz = frequency_ghz
        self.environment = environment
        self.is_indoor = is_indoor
        self.building_type = building_type

        # Set frequency band based on frequency
        if frequency_ghz < 6:
            self.frequency_band = FrequencyBand.S_BAND
        else:
            self.frequency_band = FrequencyBand.KA_BAND

    def get_nearest_elevation_angle(self, elevation_angle):
        """
        Get the nearest reference elevation angle

        This function maps a calculated elevation angle to the nearest reference angle
        (10°, 20°, 30°, etc.) used in the 3GPP tables. It doesn't calculate the actual
        elevation angle of a satellite.

        Args:
            elevation_angle (float): Elevation angle in degrees

        Returns:
            int: Nearest reference elevation angle
        """
        reference_angles = [10, 20, 30, 40, 50, 60, 70, 80, 90]
        return min(reference_angles, key=lambda x: abs(x - elevation_angle))

    def determine_los_condition(self, elevation_angle, forced_condition=None):
        """
        Determine path condition (LOS or NLOS)

        According to 3GPP 38.811 section 6.6.1, LOS probability depends on UE environment
        and elevation angle. This function can be used in three ways:
        1. Forced LOS: Always return LOS condition
        2. Forced NLOS: Always return NLOS condition
        3. Probabilistic: Return LOS or NLOS based on probability in the specification

        Args:
            elevation_angle (float): Elevation angle in degrees
            forced_condition (PathCondition): Force a specific condition (LOS or NLOS),
                                              if None, use probabilistic determination

        Returns:
            PathCondition: LOS or NLOS
        """
        if forced_condition is not None:
            return forced_condition

        nearest_angle = self.get_nearest_elevation_angle(elevation_angle)
        los_probability = self.LOS_PROBABILITY[self.environment][nearest_angle]

        # Probabilistic determination based on the LOS probability table
        if np.random.random() < los_probability:
            return PathCondition.LOS
        else:
            return PathCondition.NLOS

    def calculate_slant_range(self, elevation_angle):
        """
        Calculate the slant range between satellite and ground terminal

        Args:
            elevation_angle (float): Elevation angle in degrees

        Returns:
            float: Slant range in meters
        """
        alpha_rad = np.radians(elevation_angle)
        distance = np.sqrt(self.EARTH_RADIUS ** 2 * np.sin(alpha_rad) ** 2 +
                           2 * self.EARTH_RADIUS * self.satellite_altitude +
                           self.satellite_altitude ** 2) - self.EARTH_RADIUS * np.sin(alpha_rad)
        return distance

    def calculate_fspl(self, distance):
        """
        Calculate the free space path loss

        Args:
            distance (float): Distance in meters

        Returns:
            float: Free space path loss in dB
        """
        # FSPL = 32.45 + 20*log10(f) + 20*log10(d)
        # where f is in MHz and d is in km
        f_mhz = self.frequency_ghz * 1000
        d_km = distance / 1000
        return 32.45 + 20 * np.log10(f_mhz) + 20 * np.log10(d_km)

    def get_shadow_fading_and_clutter_loss(self, elevation_angle, path_condition):
        """
        Get shadow fading and clutter loss values

        Args:
            elevation_angle (float): Elevation angle in degrees
            path_condition (PathCondition): LOS or NLOS condition

        Returns:
            tuple: (shadow_fading_std, clutter_loss)
        """
        nearest_angle = self.get_nearest_elevation_angle(elevation_angle)
        sf_std, cl = self.SF_CL_TABLES[self.environment][self.frequency_band][path_condition][nearest_angle]
        return sf_std, cl

    def calculate_o2i_penetration_loss(self, elevation_angle, probability=0.5):
        """
        Calculate building entry loss for indoor UEs using ITU-R P.2109

        Args:
            elevation_angle (float): Elevation angle in degrees
            probability (float): Probability that the loss is not exceeded (0 < p < 1)

        Returns:
            float: Building entry loss in dB
        """
        if not self.is_indoor or self.building_type is None:
            return 0.0

        # Get coefficients
        coef = self.O2I_COEFFICIENTS[self.building_type]

        # Calculate median loss for horizontal paths (equation 6.6-6)
        Lh = coef['r'] + coef['s'] * np.log10(self.frequency_ghz) + coef['t'] * coef['u']

        # Calculate correction for elevation angle (equation 6.6-7)
        Le = coef['v'] * (1 - np.exp(-elevation_angle / coef['w'])) + \
             coef['x'] * (1 - np.exp(-self.frequency_ghz / coef['y'])) + \
             coef['z']

        # Calculate combined loss (equation 6.6-5)
        F_inv = norm.ppf(probability)

        mu1 = Lh + Le
        mu2 = coef['u']
        sigma1 = coef['v']
        sigma2 = coef['y']

        p1 = 1 / (1 + (probability / (1 - probability)) * np.exp(
            F_inv * np.sqrt(np.log(10) ** 2 / 100 * (sigma2 ** 2 - sigma1 ** 2))))
        loss = mu1 + F_inv * sigma1 * np.sqrt(
            np.log(10) ** 2 / 100) if probability <= p1 else mu2 + F_inv * sigma2 * np.sqrt(np.log(10) ** 2 / 100)

        return loss

    def calculate_atmospheric_absorption(self, elevation_angle):
        """
        Calculate attenuation due to atmospheric gases according to ITU-R P.676

        Args:
            elevation_angle (float): Elevation angle in degrees

        Returns:
            float: Atmospheric absorption in dB
        """
        # Get closest frequency for zenith attenuation
        frequencies = list(self.ATMOSPHERIC_ATTENUATION_ZENITH.keys())
        closest_freq = min(frequencies, key=lambda x: abs(x - self.frequency_ghz))

        # Get zenith attenuation and apply elevation angle correction
        zenith_attenuation = self.ATMOSPHERIC_ATTENUATION_ZENITH[closest_freq]

        # Apply frequency scaling if not exact match
        if closest_freq != self.frequency_ghz:
            # Simple linear interpolation
            if self.frequency_ghz > closest_freq and frequencies.index(closest_freq) < len(frequencies) - 1:
                next_freq = frequencies[frequencies.index(closest_freq) + 1]
                next_atten = self.ATMOSPHERIC_ATTENUATION_ZENITH[next_freq]
                zenith_attenuation = zenith_attenuation + (next_atten - zenith_attenuation) * \
                                     (self.frequency_ghz - closest_freq) / (next_freq - closest_freq)
            elif self.frequency_ghz < closest_freq and frequencies.index(closest_freq) > 0:
                prev_freq = frequencies[frequencies.index(closest_freq) - 1]
                prev_atten = self.ATMOSPHERIC_ATTENUATION_ZENITH[prev_freq]
                zenith_attenuation = prev_atten + (zenith_attenuation - prev_atten) * \
                                     (self.frequency_ghz - prev_freq) / (closest_freq - prev_freq)

        # Apply elevation angle correction (equation 6.6-8)
        attenuation = zenith_attenuation / np.sin(np.radians(elevation_angle))

        # Limit to reasonable values
        return min(attenuation, 10.0)  # Cap at 10 dB to avoid extreme values at very low elevation angles

    def calculate_scintillation_loss(self, elevation_angle):
        """
        Calculate attenuation due to scintillation (ionospheric or tropospheric)

        Args:
            elevation_angle (float): Elevation angle in degrees

        Returns:
            float: Scintillation loss in dB
        """
        nearest_angle = self.get_nearest_elevation_angle(elevation_angle)

        # Ionospheric scintillation (for frequencies below 6 GHz)
        if self.frequency_band == FrequencyBand.S_BAND:
            # Simplified model based on section 6.6.6.1.4
            scintillation_4ghz = 0.76  # At 99% of the P3 curve for 4 GHz
            # Scale to actual frequency (equation 6.6-11)
            scintillation = scintillation_4ghz * (self.frequency_ghz / 4.0) ** (-1.5)
            return scintillation

        # Tropospheric scintillation (for frequencies above 6 GHz)
        elif self.frequency_band == FrequencyBand.KA_BAND:
            scintillation_20ghz = self.TROPOSPHERIC_SCINTILLATION[nearest_angle]
            # Scale to actual frequency (tropospheric scintillation increases with frequency)
            scintillation = scintillation_20ghz * (self.frequency_ghz / 20.0) ** (7 / 12)
            return scintillation

        return 0.0

    def calculate_path_loss(self, elevation_angle, forced_condition=None):
        """
        Calculate the total path loss according to 3GPP 38.811 section 6.6.2

        Args:
            elevation_angle (float): Elevation angle in degrees
            forced_condition (PathCondition): Force a specific condition (LOS or NLOS),
                                             if None, use probabilistic determination

        Returns:
            dict: Path loss components and total path loss
        """
        # Determine path condition (LOS or NLOS)
        path_condition = self.determine_los_condition(elevation_angle, forced_condition)

        # Calculate slant range (equation 6.6-3)
        distance = self.calculate_slant_range(elevation_angle)

        # Calculate free space path loss (equation 6.6-2)
        fspl = self.calculate_fspl(distance)

        # Get shadow fading and clutter loss values from tables 6.6.2-1/2/3
        sf_std, cl = self.get_shadow_fading_and_clutter_loss(elevation_angle, path_condition)

        # Generate shadow fading (log-normal distribution with σ from tables)
        sf = np.random.normal(0, sf_std)

        # Calculate basic path loss (equation 6.6-4)
        # When the UE is in LOS condition, clutter loss is negligible and should be set to 0 dB
        if path_condition == PathCondition.LOS:
            cl = 0
        basic_pl = fspl + cl + sf

        # Calculate O2I penetration loss (section 6.6.3)
        o2i_loss = self.calculate_o2i_penetration_loss(elevation_angle)

        # Calculate atmospheric absorption (section 6.6.4)
        atmos_loss = self.calculate_atmospheric_absorption(elevation_angle)

        # Calculate scintillation loss (section 6.6.6)
        scint_loss = self.calculate_scintillation_loss(elevation_angle)

        # Calculate total path loss (equation 6.6-1)
        total_pl = basic_pl + atmos_loss + scint_loss + o2i_loss

        # Return components and total
        # return total_pl
        return {
            'path_condition': path_condition,
            'distance': distance,
            'free_space_path_loss': fspl,
            'clutter_loss': cl,
            'shadow_fading': sf,
            'basic_path_loss': basic_pl,
            'o2i_penetration_loss': o2i_loss,
            'atmospheric_loss': atmos_loss,
            'scintillation_loss': scint_loss,
            'total_path_loss': total_pl
        }

class LinkBudget:
    """Link budget calculator for NTN communications"""

    def __init__(self, channel_model, tx_power_dbm, tx_gain_dbi, rx_gain_dbi, noise_figure_db,
               bandwidth_hz, link_direction=LinkDirection.DOWNLINK, add_interference=False):
        """
        Initialize the link budget calculator

        Args:
            channel_model (NTNChannelModel): NTN channel model
            tx_power_dbm (float): Transmit power in dBm
            tx_gain_dbi (float): Transmit antenna gain in dBi
            rx_gain_dbi (float): Receive antenna gain in dBi
            noise_figure_db (float): Receiver noise figure in dB
            bandwidth_hz (float): Signal bandwidth in Hz
            link_direction (LinkDirection): UPLINK or DOWNLINK
            add_interference (bool): Whether to add interference
        """
        self.channel_model = channel_model
        self.tx_power_dbm = tx_power_dbm
        self.tx_gain_dbi = tx_gain_dbi
        self.rx_gain_dbi = rx_gain_dbi
        self.noise_figure_db = noise_figure_db
        self.bandwidth_hz = bandwidth_hz
        self.link_direction = link_direction
        self.add_interference = add_interference

    def calculate_thermal_noise_dbm(self):
        """
        Calculate thermal noise power

        Returns:
            float: Thermal noise power in dBm
        """
        # N = k * T * B
        k = 1.38e-23  # Boltzmann constant
        T = 290  # Temperature in Kelvin
        noise_power_watts = k * T * self.bandwidth_hz
        noise_power_dbm = 10 * np.log10(noise_power_watts * 1000)
        return noise_power_dbm

    def calculate_interference_dbm(self, elevation_angle):
        """
        Calculate interference power (simplified model)

        Args:
            elevation_angle (float): Elevation angle in degrees

        Returns:
            float: Interference power in dBm
        """
        if not self.add_interference:
            return -float('inf')  # No interference

        # Simplified model: interference decreases with elevation angle
        # and is proportional to received signal power
        path_loss_results = self.channel_model.calculate_path_loss(elevation_angle)
        path_loss = path_loss_results['total_path_loss']
        received_power = self.tx_power_dbm + self.tx_gain_dbi + self.rx_gain_dbi - path_loss

        # Interference is modeled as 10-20 dB below signal power, stronger at low elevation angles
        interference_offset = -10 - 10 * (elevation_angle / 90.0)
        return received_power + interference_offset

    def calculate_snr(self, elevation_angle, forced_condition=None):
        """
        Calculate Signal-to-Noise Ratio (SNR)

        Args:
            elevation_angle (float): Elevation angle in degrees
            forced_condition (PathCondition): Force a specific path condition

        Returns:
            tuple: (SNR in dB, path loss components)
        """
        # Calculate path loss
        path_loss_results = self.channel_model.calculate_path_loss(elevation_angle, forced_condition)
        path_loss = path_loss_results['total_path_loss']

        # Calculate received power
        received_power_dbm = self.tx_power_dbm + self.tx_gain_dbi + self.rx_gain_dbi - path_loss

        # Calculate noise power
        noise_power_dbm = self.calculate_thermal_noise_dbm() + self.noise_figure_db

        # Calculate SNR
        snr_db = received_power_dbm - noise_power_dbm

        return snr_db, path_loss_results

    def calculate_sinr(self, elevation_angle, forced_condition=None):
        """
        Calculate Signal-to-Interference-plus-Noise Ratio (SINR)

        Args:
            elevation_angle (float): Elevation angle in degrees
            forced_condition (PathCondition): Force a specific path condition

        Returns:
            tuple: (SINR in dB, path loss components)
        """
        # Calculate SNR
        snr_db, path_loss_results = self.calculate_snr(elevation_angle, forced_condition)

        # Calculate interference power
        interference_dbm = self.calculate_interference_dbm(elevation_angle)

        # Calculate noise power
        noise_power_dbm = self.calculate_thermal_noise_dbm() + self.noise_figure_db

        # Calculate total interference plus noise power
        if interference_dbm == -float('inf'):
            # No interference, SINR = SNR
            sinr_db = snr_db
        else:
            # Convert to linear, add, then convert back to dB
            noise_linear = 10 ** (noise_power_dbm / 10)
            interference_linear = 10 ** (interference_dbm / 10)
            total_noise_interference_linear = noise_linear + interference_linear
            total_noise_interference_dbm = 10 * np.log10(total_noise_interference_linear)

            # Calculate received power
            path_loss = path_loss_results['total_path_loss']
            received_power_dbm = self.tx_power_dbm + self.tx_gain_dbi + self.rx_gain_dbi - path_loss

            # Calculate SINR
            sinr_db = received_power_dbm - total_noise_interference_dbm

        return sinr_db, path_loss_results


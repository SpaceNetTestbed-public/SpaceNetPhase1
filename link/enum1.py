from enum import Enum
# import sys
# print(sys.version)  # Check your Python version
#
# try:
#    import enum
#    print("Enum module found at:", enum.__file__)
# except ImportError:
#    print("Enum module not found")

class Environment(Enum):
    """UE environment types"""
    DENSE_URBAN = 1
    URBAN = 2
    SUBURBAN_RURAL = 3

class FrequencyBand(Enum):
    """Frequency bands"""
    S_BAND = 1  # ~2 GHz
    KA_BAND = 2  # ~20-30 GHz

class LinkDirection(Enum):
    """Link direction"""
    UPLINK = 1
    DOWNLINK = 2

class PathCondition(Enum):
    """Path condition"""
    LOS = 1
    NLOS = 2

class BuildingType(Enum):
    """Building types for O2I penetration loss"""
    TRADITIONAL = 1
    THERMALLY_EFFICIENT = 2


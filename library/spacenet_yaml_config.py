import yaml
import datetime

def load_sim_and_constellation_config_file(sim_config_path, sim_config_file_name, constellation_config_subdir):
    sim_config = read_config(sim_config_path + sim_config_file_name)
    if sim_config is None:
        print(f"Error: Unable to load simulation configuration file '{sim_config_path + sim_config_file_name}'")
        return None, None
    constellationName = sim_config['ConstellationName']
    constellation_config = read_config(sim_config_path + constellation_config_subdir + constellationName + '.yaml')
    if constellation_config is None:
        print(f"Error: Unable to load constellation configuration file '{sim_config_path + constellation_config_subdir + constellationName + '.yaml'}'")
        return None, None
    constellation_config["EpochStartDateTime"] = datetime.datetime(int(constellation_config["EpochStartYear"]), int(constellation_config["EpochStartMonth"]), int(constellation_config["EpochStartDay"]), int(constellation_config["EpochStartHour"]), int(constellation_config["EpochStartMinute"]), int(constellation_config["EpochStartSecond"]))
    return sim_config, constellation_config

def read_config(file_path):
    try:
        with open(file_path, 'r') as stream:
            config = yaml.safe_load(stream)
    except Exception as exc:
        print(exc)
        config = None
    print(f"Loaded configuration file '{file_path}'")
    return config
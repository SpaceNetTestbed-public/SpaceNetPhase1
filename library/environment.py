import re
import warnings
from datetime import datetime
from astropy.time import Time
from library.spacenet_yaml_config import load_sim_and_constellation_config_file

class simulation_variables:

    def __init__(self, config_file_path:str=None, config_file_name:str=None, sat_config_sub_path:str=None):
        '''
        TODO: Instead of main_config and sat_config inputs take the paths as inputs and call spacenet_yaml_config here! 
        Need to change logic in main for this (MINOR CHANGE)
        '''
        class SimTime:

            def __init__(time_self, year:int=1970, month:int=1, day:int=1, hour:int=0, minute:int=0, second:int=0):
                time_self.year = year
                time_self.month = month
                time_self.day = day
                time_self.hour = hour
                time_self.minute = minute
                time_self.second = second

                time_self.t_datetime = datetime(time_self.year, time_self.month, time_self.day, time_self.hour, time_self.minute, time_self.second)
                time_self.t_Time = Time(time_self.t_datetime, format='datetime', scale='utc')
                time_self.t_jd = time_self.t_Time.jd
                time_self.t_unix = time_self.t_Time.unix

        self.main_config_path = config_file_path if config_file_path else "config_files/"
        self.main_filename = config_file_name if config_file_name else "main_config.yaml"
        self.sat_config_path = self.main_config_path + sat_config_sub_path if sat_config_sub_path else self.main_config_path + "sat_config_files/"
        m_config, s_config, PATHS = load_sim_and_constellation_config_file(self.main_config_path, self.main_filename, self.sat_config_path)
        self.ALL_CONFIG_PATHS = PATHS
        self.M_CONFIG = m_config
        self.S_CONFIG = s_config

        self.boolmap = lambda x: {"false":0, "true":1}[x.lower()]
        self.simtime = SimTime(int(s_config['Sim_Date_Time']['StartYear']), int(s_config['Sim_Date_Time']['StartMonth']), int(s_config['Sim_Date_Time']['StartDay']), int(s_config['Sim_Date_Time']['StartHour']), int(s_config['Sim_Date_Time']['StartMinute']), int(s_config['Sim_Date_Time']['StartSecond']))
        self.dt = s_config['Sim_Length']['TimeStepDuration']
        self.N = s_config['Sim_Length']['TimeStepCount']
        self.generate_TLE = s_config['generate_TLE']
        self.N_shells = len(s_config['shells'].keys())
        self.store_shell_data(s_config['shells']) # Initializes self.shell_data that also contains self.bodymodel, self.perturbermodel types

        self.get_operator()
        self.tle_path = s_config['TLEFilePath'] if 'TLEFilePath' in s_config.keys() else "utils/" + self.operator_name + "_tles"
        self.isl_type = m_config['ISL_Type']
        self.output_path = m_config['OutputFilePath']
        self.monitor_resource = bool(m_config['MonitorResource'])
        self.source_node = int(m_config['SourceNode'])
        self.dest_node = int(m_config['DestNode'])
        self.routing_metric = m_config['RouteWeight']
        self.gs_data_loc = m_config['GroundStationFile']
        self.min_elev = float(m_config['min_elevation_angle'])
        self.gsl_criterion = m_config['AssociationCritGSL']
        self.use_weather = bool(m_config['UseWeatherData'])
        self.topo_criterion = int(m_config['TopoCrit'])
        self.use_azure = bool(m_config['Azure']['t2t_use_azure'])
        self.use_wonderproxy = bool(m_config['WonderProxy']['t2t_use_wonderproxy'])
        self.debug = bool(m_config['Debug'])

        # self.sanity_checks()


    def get_operator(self):
        '''
        Parses operator name from main config file or satellite_config file
        '''
        self.operator_name = re.match(r'[a-zA-Z]+', self.M_CONFIG["ConstellationName"]).group(0)
        if 'operator_name' in self.S_CONFIG:
            self.operator_name = re.match(r'[a-zA-Z]+', self.S_CONFIG["operator_name"]).group(0)

    
    def store_shell_data(self, s_conf:dict) -> None:
        '''
        Method to store shell data directly from the config files. Calls SimSatShells class to store data.
        '''
        class SimSatShells:
        
            def __init__(self, s_config:dict):

                class CelesModel:
                    def __init__(inner_self, info:dict):
                            inner_self.name = info['name']
                            inner_self.model = info['model']
                            inner_self.pck = info['pck']
    
                    def update_pck(inner_self, new_path):
                        '''
                        Updates PCK
                        '''
                        inner_self.pck = new_path
    
                    def update_model(inner_self, new_path):
                        '''
                        Updates model
                        '''
                        inner_self.pck = new_path
    
                class BodyModel(CelesModel):
                    def __init__(inner_self, info:dict):
                        super().__init__(info)
    
                class PertModel(CelesModel):
                    def __init__(inner_self, info:dict):
                        pert_info = info
                        pert_info['degree'] = None
                        pert_info['order'] = None
                        pert_info['zonal_only'] = False
                        super().__init__(pert_info)

                self.name = s_config['name']
                self.orbits = s_config['orbits']
                self.sat_per_orbits = s_config['sat_per_orbit']
                self.altitude = s_config['altitude']
                self.inc = s_config['inclination']
                self.pattern = s_config['pattern']
                self.ipp = s_config['ipp_increment']
                self.bodymodel = BodyModel(s_config['body_model']) if 'body_model' in s_config.keys() else s_config['body']
                self.perturbermodel = PertModel(s_config['perturber']) if 'perturber' in s_config.keys() else None

        shells = []
        for sh in s_conf.values():
            shells.append(SimSatShells(sh))
        self.shell_data = shells


    def sanity_checks(self):
        '''
        The method responsible to perform sanity checks (incompatibilities within provided configuration)
        '''
        pass
        


    def load_environment(self):
        '''
        Loads the simulators configuration into the system environment!
        '''
        pass
# Dynamic Topology Generator (SUBSET OF BBARBOUR_JP BRANCH)



## Getting started

To run a basic experiment on the default constellation (Starlink) follow these steps:
- [Make sure your files exists](#checking-your-files)
- Run main.py file 

To design your own custom experiment on an arbitrary starlink TLE:
- Run `sh get_tles.sh` to extract an actual TLE or generate a custom TLE by [setting up the TLE generator](#tle-generator). Store this TLE file at your desired location or the default location (utils/starlink_tles/) 
- [Make sure your files exists](#checking-your-files)
- Run main.py file

## Checking your files
A successful experiment runs when all the supporting files are located correctly. Look for the following files:
- Your constellation configuration file (constel_config) at `config_files/sat_config_files/`
- Your constel_config mentioned in the main_config.yaml at `config_files/`
- Your TLE file with correct filename (unix timestamp) corresponding to your experiment's datetime at `utils/starlink_tles/` or TLE location at your constel_config.
- Mention the location of your ground station file in the main_config.yaml or keep it unchanged to the default location.

## TLE generator
The testbed can generate custom constellation and can even emulate network topology over it. To make your own built, its as simple 
as defining your shell specs in your constel_config file and enable `generate_TLE = true`. You can also add multiple shells in your custom constellation.

## T2T Links
For t2t links to work, you need to have the following configurations in the main_mn_config.yaml file:
Use_t2t: set to 'True'
t2t_gateway_kmz_type: set to 'local' if using a KMZ/KML file that contains all the data; use 'link' if KMZ data has URL to pull data from the Internet
t2t_gateway_kmz_path: the location of the KMZ file for the script to retrieve the Gateway data (download from Unofficial Starlink Gateway map)
t2t_dict_output_file: location to save the Gateway dictionary after the data has been pulled from the KMZ file

t2t_use_azure: set to 'True' to use Azure data center data for Endpoints
t2t_azure_endpoint_location_file: location of csv file that provides location data for each data center
t2t_azure_endpoint_latency_url: URL to use to scrape latency data
t2t_azure_dict_output_file: location to save the Endpoint dictionary after the data has been scraped and compiled

t2t_use_wonderproxy: set to 'True' to use WonderProxy server data for Endpoints
t2t_wonderproxy_endpoint_location_file: location of csv file that provides location data for WonderProxy servers
t2t_wonderproxy_endpoint_latency_file: location of csv file that provides latency data
t2t_wonderproxy_dict_output_file: location to save the Endpoint dictionary after the data has been scraped and compiled

## Project status
If your system runs out of energy or time for your project, put a note at the top of the README saying that development has slowed down or stopped completely. Someone may choose to fork your project or volunteer to step in as a maintainer or owner, allowing your project to keep going. You can also make an explicit request for maintainers.

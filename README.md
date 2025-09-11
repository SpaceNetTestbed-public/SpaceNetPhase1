# SpaceNet Phase1 (Simulation Phase)

Welcome to the public codebase for the Phase1 block of Virginia Tech's space network emulator *SpaceNet*. This is the simulation phase of the testbed that designs the network topology for the 
SpaceNet can be used to design your own experiments on space constellation designs and Mininet-based realistic virtual network performance trends on communication use-case scenarios between any nodes in the system. By default, a low-fidelity experiment between source and destination stations only requires their location in geograhic coordinate system (lat, long) added in the Ground Station files, and the testbed would instantiate such a ground station with the default communication module specifications provided in `DATASHEET.txt`. However in the upcoming release, any user can also provide their own set of wireless parameters for their devices to be added in the simulations. The Phase1's intention is to provide a connected graph-type network topology over which the TCP/IP protocol can be administered. Therefore, in the current release, the testbed demands a well-connected constellation with essentially no isolated graphs at any time instance. This limits the usage of the testbed to only near-circular evenly-spread constellations with sufficient satellites (Starlink, Kuiper, Iridium). However from the upcoming release onwards the testbed would support deep-space high-eccentric unevenly-spaced constellations that can function with a store & forward mechanism ([Bundle protocol](https://github.com/nasa/bp)) for the packet transfers that's not available due to the current IP support. 

For the realism, the upcoming version would also introduce the availibility of traffic emulation features. A basic burst-aggregation traffic model would be added for a congestion-aware routing algorithm to produce a static routing table in a centralised congesiton routing manner. Moreover, the Phase2 (Emulation Phase) would also have its own low-level implementation of dynamic routing on individual satellite nodes if the static routing table is not sufficient enough to provide data on actual traffic. Please refer Phase2 codebase [here](https://github.com/SpaceNetTestbed-public/SpaceNetPhase2). 

## Getting started

To run a basic experiment on the default constellation (Starlink) follow these steps:
- [Make sure your files exists](#checking-your-files)
- [Tune the parameters for your experiment](#table-of-parameters)
- Run main.py file 

To design your own custom experiment on an arbitrary starlink TLE:
- Run `sh get_tles.sh` to extract an actual TLE or generate a custom TLE by [setting up the TLE generator](#tle-generator). Store this TLE file at your desired location or the default location (utils/starlink_tles/) 
- [Make sure your files exists](#checking-your-files)
- [Tune the parameters for your experiment](#table-of-parameters)
- Run main.py file

## Table of parameters
Important parameters that the users can toggle specific to their testbed experiments:

| Parameters  | Options | Definition |
| ------------- | ------------- | ------------- |
| ConstellationName  | string value  | constel_config filename  |
| Debug  | 0/1  | Verbose mode  |
| SourceNode, DestNode  | string value  | ground station names in str  |
| RouteWeight  | latency/distance/capacity/congestion  | Routing strategy  |
| AssociationCritGSL  | BASED_ON_DISTANCE_ONLY_MININET  | Ground Station Link conneciton strategy  |
| UseWeatherData  | 0/1  | Using weather information  |
| TopoCrit  | 0,1 or 2  | Enabling/disabling ISTN  |
| generate_TLE  | true/false  | generates custom TLE  |

## Checking your files
A successful experiment runs when all the supporting files are located correctly. Look for the following files:
- Your constellation configuration file (constel_config) at `config_files/sat_config_files/`
- Your constel_config mentioned in the main_config.yaml at `config_files/`
- Your TLE file with correct filename (unix timestamp) corresponding to your experiment's datetime at `utils/starlink_tles/` or TLE location at your constel_config. If it doesn't exist then you can also [generate TLEs](#tle-generator) from scratch.
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

## Ground device specifications 
Multiple device type support with their respective wireless design values.

|  | VSAT | Starlink | Handheld UEs | Gateways | IOT (class 1) | IOT (class 2) | IOT (class 3) |
| ------------- | ------------- | ------------- | ------------- | ------------- | ------------- | ------------- | ------------- |
| Transmission power (dBm)  | 33 | 47  | 23  | 23  | 14  | 20  | 23  |
| Antenna Type  | 60cm aperture dia  | phased array  | omnidirectional antenna  | omnidirectional antenna  | omnidirectional antenna  | omnidirectional antenna  | omnidirectional antenna  |
| TX Gain (dBi)  | 43.2  | 33.5  | 0  | 0  | 0  | 0  | 0  |
| RX Gain (dBi)  | 39.7  | 33.0  | 0  | 0  | 0  | 0  | 0  |
| Noise figure | 1.2  | 2.5  | 9  | 9  | 9  | 9  | 9  |
| RX Cable Loss  | 3  | 1.5  | 0  | 3  | 0  | 0  |0  |
| Polarization  | 0  | 0  | 3  | 0  | 3  | 3  | 3  |

## Project status
If your system runs out of energy or time for your project, put a note at the top of the README saying that development has slowed down or stopped completely. Someone may choose to fork your project or volunteer to step in as a maintainer or owner, allowing your project to keep going. You can also make an explicit request for maintainers.




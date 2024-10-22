from utils.file_utils import *

# Azure variables
azure_latency_json_filename = 'gs_files/azure_latency_data.json'
azure_csv_filename = 'gs_files/AzureDataCenterLocations.csv'
azure_dict_filename = 'gs_files/azure_server_dict.json'
url = 'https://learn.microsoft.com/en-us/azure/networking/azure-network-latency'
sectionIdList = ['tabpanel_2_WestUS_Americas', 'tabpanel_2_EastUS_Americas', 'tabpanel_2_CentralUS_Americas', 'tabpanel_2_Canada_Americas', 'tabpanel_2_Australia_APAC', 'tabpanel_2_Japan_APAC', 'tabpanel_2_WesternEurope_Europe', 'tabpanel_2_CentralEurope_Europe', 'tabpanel_2_NorwaySweden_Europe', 'tabpanel_2_UKNorthEurope_Europe', 'tabpanel_2_Korea_APAC', 'tabpanel_2_India_APAC', 'tabpanel_2_Asia_APAC', 'tabpanel_2_israel-qatar-uae_MiddleEast', 'tabpanel_2_southafrica_MiddleEast']
# Wonderproxy variables
# WonderProxy files downloaded from: https://wonderproxy.com/blog/a-day-in-the-life-of-the-internet/
wonderproxy_server_csv_filename = 'gs_files/wonderproxy_servers-2020-07-19.csv'
wonderproxy_latency_csv_filename = 'gs_files/wonderproxy_pings-2020-07-19-2020-07-20.csv'
wonderproxy_latency_json_filename = 'gs_files/wonderproxy_latency_data.json'
wonderproxy_dict_filename = 'gs_files/wonderproxy_server_dict.json'
# Unofficial Starlink Global Gateways and PoPs KMZ file variables
kmz_file = "gs_files/UnofficialStarlinkGlobalGatewaysNPoPs_noLinks.kmz"
usable_bands = ['Ka', 'E']

adjacency_threshold = 2000 # Threshold distance in kilometers for considering two points as adjacent (be within 5ms latency (10ms round-trip))
t2t_dict_filename = 'gs_files/t2t_dict.json'
topology_filename = 'gs_files/t2t_topology.txt'

terrestrial_link_bandwidth = 10000 # Bandwidth in Mbps for terrestrial links (10 Gbps)

# >>>>>>>>> WonderProxy Functions <<<<<<<<<<<<
def load_wonderproxy_server_dict(wonderproxy_server_csv_filename, wonderproxy_latency_csv_filename, wonderproxy_dict_filename = None):
    if wonderproxy_dict_filename is not None:
        wonderproxy_server_dict = read_json_file(wonderproxy_dict_filename)
        if wonderproxy_server_dict is not None:
            return wonderproxy_server_dict
    wonderproxy_server_dict = load_wonderproxy_server_location_coordinates(wonderproxy_server_csv_filename)
    wonderproxy_latency_dict = load_wonderproxy_latency_dict(None, wonderproxy_latency_csv_filename, wonderproxy_server_dict)
    unavail_server_list = []
    for server in wonderproxy_server_dict:
        if server not in wonderproxy_latency_dict:
            if server not in unavail_server_list:
                unavail_server_list.append(server)
            continue
        wonderproxy_server_dict[server].update(wonderproxy_latency_dict[server]) # Add latency data to server dictionary
    print(f"{len(unavail_server_list)} servers not found in latency dictionary: {unavail_server_list}")
    # Save the data to a json file for future use
    if wonderproxy_dict_filename is not None:
        write_json_file(wonderproxy_dict_filename, wonderproxy_server_dict)
    return wonderproxy_server_dict

def load_wonderproxy_server_location_coordinates(csv_filename):
    import csv
    if csv_filename is None:
        raise ValueError("csv_filename must be provided.")
    encoding = detect_text_file_encoding(csv_filename)
    if encoding is None:
        raise ValueError("Failed to detect the encoding of the csv file.")
    wonderproxy_server_dict = {}
    with open(csv_filename, 'r', encoding=encoding) as file:
        reader = csv.reader(file)
        next(reader) # Skip the header row
        for row in reader:
            locationCountry = row[5]
            if locationCountry == 'United States':
                locationName = row[3] + ', ' + row[4]
            else:
                locationName = row[3]
            latitude = float(row[-2])
            longitude = float(row[-1])
            index = int(row[0])
            wonderproxy_server_dict[locationName] = {'coordinates': (latitude, longitude), 'wp_index': index, 'data_source': 'WonderProxy'}
    return wonderproxy_server_dict

def gen_wonderproxy_server_latency_json_from_csv(csv_filename, json_filename, wonderproxy_server_dict):
    if csv_filename is None or wonderproxy_server_dict is None:
        raise ValueError(f"(gen_wonderproxy_server_latency_json_from_csv) csv_filename and server dict must be provided. Received csv_filename: {csv_filename}, server dict: {wonderproxy_server_dict}")
    # Build wonderproxy_server_dict lookup table by index
    wonderproxy_server_lookup_dict = {}
    for locationName, data in wonderproxy_server_dict.items():
        index = data['wp_index']
        wonderproxy_server_lookup_dict[index] = locationName
    wonderproxy_latency_data = read_csv_file(csv_filename)
    if wonderproxy_latency_data is None:
        raise ValueError(f"Failed to read data from {csv_filename}")
    wonderproxy_latency_dict = {}
    unavail_index_list = []
    for row in wonderproxy_latency_data[1:]:
        source_index = int(row[0])
        if source_index not in wonderproxy_server_lookup_dict:
            if source_index not in unavail_index_list:
                unavail_index_list.append(source_index)
            continue
        dest_index = int(row[1])
        if dest_index not in wonderproxy_server_lookup_dict:
            if dest_index not in unavail_index_list:
                unavail_index_list.append(dest_index)
            continue
        latency = float(row[4]) # Avg latency in ms
        source_name = wonderproxy_server_lookup_dict[source_index]
        dest_name = wonderproxy_server_lookup_dict[dest_index]
        if source_name not in wonderproxy_latency_dict:
            wonderproxy_latency_dict[source_name] = {}
        wonderproxy_latency_dict[source_name][dest_name] = latency
        if dest_name not in wonderproxy_latency_dict:
            wonderproxy_latency_dict[dest_name] = {}
        wonderproxy_latency_dict[dest_name][source_name] = latency
    print(f"{len(unavail_index_list)} indices not found in server dictionary: {unavail_index_list}")
    if json_filename is not None:
        write_json_file(json_filename, wonderproxy_latency_dict)
        print(f"Data saved to {json_filename}")
    return wonderproxy_latency_dict

def load_wonderproxy_latency_dict(wonderproxy_latency_dict_filename = None, csv_filename = None, wonderproxy_server_dict = None):
    wonderproxy_latency_dict = None
    if wonderproxy_latency_dict_filename is None and csv_filename is None:
        raise ValueError("(load_wonderproxy_latency_dict) Either csv_filename or json_filename must be provided.")
    if wonderproxy_latency_dict_filename is not None:
        wonderproxy_latency_dict = read_json_file(wonderproxy_latency_dict_filename)
    if wonderproxy_latency_dict is None and csv_filename is None:
        raise ValueError("(load_wonderproxy_latency_dict) Cached data not found and csv_filename is not provided.")
    elif wonderproxy_latency_dict is None:
        wonderproxy_latency_dict = gen_wonderproxy_server_latency_json_from_csv(csv_filename, wonderproxy_latency_dict_filename, wonderproxy_server_dict)
    return wonderproxy_latency_dict

# >>>>>>>>> MS Azure Functions <<<<<<<<<<<<
# Function that reads AzureDataCenterLocation csv and adds location name and coordinates to the dictionary
def add_azure_location_coordinates(azure_data_centers, csv_filename):
    import csv
    if csv_filename is None:
        raise ValueError("csv_filename must be provided.")
    encoding = detect_text_file_encoding(csv_filename)
    if encoding is None:
        raise ValueError("Failed to detect the encoding of the csv file.")
    unavail_server_list = []
    with open(csv_filename, 'r', encoding=encoding) as file:
        reader = csv.reader(file)
        next(reader) # Skip the header row
        for row in reader:
            centerName = row[0]
            locationName = row[1]
            latitude = float(row[2])
            longitude = float(row[3])
            if centerName in azure_data_centers:
                azure_data_centers[centerName]['coordinates'] = (latitude, longitude)
                azure_data_centers[centerName]['friendly_name'] = locationName
            else:
                if centerName not in unavail_server_list:
                    unavail_server_list.append(centerName)
    print(f"{len(unavail_server_list)} Azure data centers not found in azure_data_centers dictionary: {unavail_server_list}")
    return azure_data_centers

def scrape_azure_latency_data(url, section_id='tabpanel_2_WestUS_Americas'):
    import requests
    from bs4 import BeautifulSoup
    print(f"Scraping data for section {section_id}")
    response = requests.get(url)
    response.raise_for_status() # Raise an exception for 4xx/5xx status codes

    soup = BeautifulSoup(response.text, 'html.parser')
    latency_data = []

    # Find the section that contains the table for latency data
    section = soup.find('section', {'id': section_id})
    if not section:
        print(f"No section found for section {section_id}")
        return []
    table = section.find('table')
    if not table:
        print(f"No table found for section {section_id}")
        return []
    
    # Extract table headers
    headers = [header.text for header in table.find_all('th')]
    # Extract table rows
    rows = table.find_all('tr')
    for row in rows[1:]: # Skip the first row (headers)
        cells = row.find_all('td')
        if len(cells) == len(headers):
            data = {headers[i]: cells[i].get_text(strip=True) for i in range(len(headers))}
            latency_data.append(data)
    return latency_data

def load_azure_data_center_latency(azure_latency_json_filename = None, url = None):
    if azure_latency_json_filename is None and url is None:
        raise ValueError("Either csv_filename or url must be provided.")
    azure_data_centers = None
    if azure_latency_json_filename is not None:
        azure_data_centers = read_json_file(azure_latency_json_filename)
    if azure_data_centers is None and url is None:
        raise ValueError("Cached data not found and url is not provided.")
    elif azure_data_centers is None:
        print("No cached data found. Scraping data from Azure website.")
        azure_data_centers = {}
        latency_data = []
        for section_id in sectionIdList:
            latency_data += scrape_azure_latency_data(url, section_id)
        for data in latency_data:
            source = data['Source']
            del data['Source']
            if source not in azure_data_centers:
                azure_data_centers[source] = data
            else:
                azure_data_centers[source].update(data)
        # Now that we have the data, convert the latencies to floats; any non-conforming entries will be marked for removal
        # Also use this opportunity to label the source of the data
        remove_list = []
        for sourceName, data in azure_data_centers.items():
            for destName, latencyStr in data.items():
                try:
                    data[destName] = float(latencyStr)
                except ValueError:
                    print(f"Error converting {latencyStr} to float for {sourceName} --> {destName}. Marking entry for removal.")
                    remove_list.append((sourceName, destName))
            data['data_source'] = 'Azure' # Add data source label
        print(f"Marked {len(remove_list)} entries for removal. Deleting...: ", end="")
        for sourceName, destName in remove_list:
            del azure_data_centers[sourceName][destName]
        print("Done.")

        # Save the data to a json file if filename is provided for future use
        if azure_latency_json_filename is not None:
            write_json_file(azure_latency_json_filename, azure_data_centers)
            print(f"Data saved to {azure_latency_json_filename}")
    return azure_data_centers

def load_azure_data_center_dict(azure_latency_json_filename = None, azure_url = None, csv_filename = None, azure_data_center_dict_json_filename = None):
    if azure_data_center_dict_json_filename is not None:
        azure_data_centers = read_json_file(azure_data_center_dict_json_filename)
        if azure_data_centers is not None:
            return azure_data_centers
    azure_data_centers = load_azure_data_center_latency(azure_latency_json_filename, azure_url)
    add_azure_location_coordinates(azure_data_centers, csv_filename)
    # List any azure_data_centers that do not have coordinates
    for region, data in azure_data_centers.items():
        if 'coordinates' not in data:
            print(f"No coordinates found for region: {region}")
    # Save the data to a json file for future use
    if azure_data_center_dict_json_filename is not None:
        write_json_file(azure_data_center_dict_json_filename, azure_data_centers)
    return azure_data_centers

# >>>>>>>>> Coordinate/Distance Functions <<<<<<<<<<<<
# Function to find adjacency pairs between server locations and ground stations
def find_adjacency_pairs(server_location_dict, ground_station_dict, adjacency_pair_dict = None, adjacency_distance_threshold = 2000):
    #adjacency_pairs = []
    if adjacency_pair_dict is None:
        adjacency_pair_dict = {} # will use gs_name as key; value is tuple with region, distance
    for region, data in server_location_dict.items():
        if 'coordinates' not in data:
            continue
        server_coords = data['coordinates']
        for gs_name, gs_data in ground_station_dict.items():
            gs_coords = gs_data['coordinates']
            distance = calc_coord_distance(server_coords, gs_coords)
            if distance is None:
                print(f"Error calculating distance between server {data} and ground station {gs_data}")
                exit()
            if distance <= adjacency_distance_threshold:
                if gs_name in adjacency_pair_dict:
                    if distance < adjacency_pair_dict[gs_name][1]:
                        adjacency_pair_dict[gs_name] = (region, distance)
                else:
                    adjacency_pair_dict[gs_name] = (region, distance)
                #adjacency_pairs.append((region, gs_name, distance))
    return adjacency_pair_dict #adjacency_pairs

# Function to calculate the distance between two points on Earth's surface using the Haversine formula
def haversine(lat1, lon1, lat2, lon2):
    import math
    # Radius of the Earth in kilometers
    R = 6371.0
    
    # Convert latitude and longitude from degrees to radians
    lat1 = math.radians(lat1)
    lon1 = math.radians(lon1)
    lat2 = math.radians(lat2)
    lon2 = math.radians(lon2)
    
    # Differences in coordinates
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    
    # Haversine formula
    a = math.sin(dlat / 2)**2 + math.cos(lat1) * math.cos(lat2) * math.sin(dlon / 2)**2
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))
    
    # Distance in kilometers
    distance = R * c
    
    return distance

# Function to calculate the distance between two points on Earth's surface using the Haversine formula
def calc_coord_distance(point1, point2):
    lat1, lon1 = point1
    lat2, lon2 = point2
    return haversine(lat1, lon1, lat2, lon2)
    #from geopy.distance import great_circle
    #try:
    #    distance = great_circle(point1, point2).kilometers
    #except ValueError as e:
    #    print(f"Error calculating distance between {point1} and {point2}: {e}")
    #    distance = None
    #return distance

# >>>>>>>>> KML Functions <<<<<<<<<<<<
def parse_local_kml_file(file_path, print_content=True):
    import zipfile
    import xml.etree.ElementTree as ET
    """Parse a local KMZ file and print the contents of the KML file inside."""
    with zipfile.ZipFile(file_path, 'r') as kmz:
        # List all the files in the KMZ archive
        file_list = kmz.namelist()
        print(f"(parse_local_kml_file) Files in KMZ archive: {file_list}")
        
        # Extract the KML file (assuming it's the first file in the archive)
        kml_filename = file_list[0]
        print(f"(parse_local_kml_file) Extracting and parsing KML file: {kml_filename}")
        
        kml_content = kmz.read(kml_filename)

        # Parse the KML content
        root = ET.fromstring(kml_content)
        
        if print_content:
            # Pretty print the KML content
            pretty_print_xml(kml_content.decode('utf-8'))
        
        return root
    
def pretty_print_xml(xml_string):
    """Pretty print XML string for easier reading."""
    from xml.dom import minidom
    parsed_xml = minidom.parseString(xml_string)
    pretty_xml = parsed_xml.toprettyxml(indent="  ")
    print(pretty_xml)

def extract_placemarks(root):
    """Extract and print relevant placemarks from a KMZ file."""
    ns = {'kml': 'http://www.opengis.net/kml/2.2'}

    placemark_dict = {}

    # Parse through folders and placemarks
    for folder in root.findall('.//kml:Folder', ns):
        folder_name = folder.find('kml:name', ns).text
        print(f"Folder: {folder_name}")

        for placemark in folder.findall('kml:Placemark', ns):
            placemark_name = placemark.find('kml:name', ns).text
            coordinates = placemark.find('.//kml:coordinates', ns).text.strip()
            description = placemark.find('kml:description', ns).text if placemark.find('kml:description', ns) is not None else 'No description'

            if folder_name not in placemark_dict:
                placemark_dict[folder_name] = []
            placemark_dict[folder_name].append({placemark_name: {'coordinates': coordinates, 'description': description}})
            #print(f"  Placemark: {placemark_name}")
            #print(f"    Coordinates: {coordinates}")
            #print(f"    Description: {description}")
    return placemark_dict

def get_ground_stations_from_placemarks(placemark_dict):
    """Extract ground station data from the placemark dictionary."""
    ground_station_dict = {}
    for folder_name, placemarks in placemark_dict.items():
        if folder_name == "PoPs & Backbone":
            continue
        for placemark in placemarks:
            for name, data in placemark.items():
                if 'coordinates' not in data:
                    print(f"Skipping ground station {name} without coordinates.")
                    continue
                if name == "Unnamed" or name == "Unnamed Placemark" or name == "":
                    print(f"Skipping ground station without name.")
                    continue
                operational = False
                for band in usable_bands:
                    if f"{band} Operational: TRUE" in data['description']:
                        operational = True
                        break
                if not operational:
                    continue
                coordinateStr = data['coordinates']
                #coord_split = coordinateStr.split(',')
                # Remove leading/trailing whitespace and convert to float
                coord_list = coordinateStr.split(',')
                if len(coord_list) == 1:
                    continue
                #x_coordStr, y_coordStr = coordinateStr.split(',') # in KMZ, the coordinates are in x, y
                x_coordStr, y_coordStr = coord_list[0], coord_list[1]
                lat = float(y_coordStr.strip())
                lon = float(x_coordStr.strip())
                #coordinates = [float(coord.strip()) for coord in coordinateStr.split(',')]
                coordinates = [lat, lon]
                ground_station_dict[name] = {
                    'locationName': name,
                    'coordinates': (coordinates[0], coordinates[1]),
                    'description': data['description'],
                    'data_source': 'Unofficial Starlink Global Gateways and PoPs KMZ'
                }
    return ground_station_dict

def load_gateways_from_local_kml(kmz_file):
    kml_root = parse_local_kml_file(kmz_file, print_content=False)
    placemark_dict = extract_placemarks(kml_root)
    ground_station_dict = get_ground_stations_from_placemarks(placemark_dict)
    return ground_station_dict

def genT2tDict(gateway_dict, endpoint_dict_with_latency_list, t2t_dict_output_filename = None):
    """
    Generates dictionary of nodes with the following format:
    source_id:
        {'type': 'endpoint' or 'gateway',
        'name': source_endpoint,
        'coordinates': (latitude, longitude),
        'friendly_name': friendly_name, (optional)
        'data_source': data_source,
        'wp_index': wp_index, (optional)
        dest_id: latency,
        dest_id2: latency2,
        ...}
    Later, the dictionary can be augmented with the following keys:
    'conn_mat_index_to_t2t_index': {conn_mat_index: t2t_index, ...}
    't2t_gw_to_ep_link_list': [(t2t_index, t2t_index), ...]
    
    """
    t2t_dict = {} # build dictionary as 2D array of latencies between each pair of endpoints
    aggr_endpoint_dict = {} # build dictionary of all endpoints from all dictionaries
    for endpoint_dict in endpoint_dict_with_latency_list:
        aggr_endpoint_dict.update(endpoint_dict) # --> TO DO: Consider case when endpoint names are not unique

    # build list of all endpoints from all dictionaries
    endpoint_name_list = aggr_endpoint_dict.keys()
    
    # Ensure latency values have reverse mappings to sacrifice memory for lookup performance -- WHAT TO DO REGARDING LACK OF COORDINATES?
    endpoint_non_latency_keys = ['coordinates', 'friendly_name', 'data_source', 'wp_index']
    #for sourceEndpoint in aggr_endpoint_dict.keys():
    #    for key, value in aggr_endpoint_dict[sourceEndpoint].items():
    #        if key in endpoint_non_latency_keys:
    #            continue
    #        if value == "":
    #            continue
    #        if key not in aggr_endpoint_dict:
    #            aggr_endpoint_dict[key] = {'data_source': aggr_endpoint_dict[sourceEndpoint]['data_source'] + '_derived'}
    #        aggr_endpoint_dict[key][sourceEndpoint] = value

    # Assign unique numbers to endpoints and build mappings of endpoint-to-number and number-to-endpoint
    endpoint_to_id_dict = {endpoint: id_num for id_num, endpoint in enumerate(endpoint_name_list)}
    id_to_endpoint_dict = {id_num: endpoint for id_num, endpoint in enumerate(endpoint_name_list)}
    # Loop through endpoints to build t2t_dict latency values using source_ids as keys
    for source_endpoint in endpoint_name_list:
        source_id = int(endpoint_to_id_dict[source_endpoint])
        if source_id not in t2t_dict:
            t2t_dict[source_id] = {'type': 'endpoint', 'name': source_endpoint}
            for non_latency_key in endpoint_non_latency_keys:
                if non_latency_key in aggr_endpoint_dict[source_endpoint]:
                    t2t_dict[source_id][non_latency_key] = aggr_endpoint_dict[source_endpoint][non_latency_key]
        for dest_endpoint in endpoint_name_list: # Now get latency values
            if source_endpoint == dest_endpoint:
                continue
            dest_id = int(endpoint_to_id_dict[dest_endpoint])
            if dest_endpoint in aggr_endpoint_dict[source_endpoint]:
                latency = aggr_endpoint_dict[source_endpoint][dest_endpoint]
                t2t_dict[source_id][dest_id] = latency
                #if dest_id not in t2t_dict:
                #    t2t_dict[dest_id] = {'type': 'endpoint', 'name': dest_endpoint}
                #t2t_dict[dest_id][source_id] = latency # sacrifice memory for lookup performance by adding in reverse direction
    # Find endpoints that are within adjacency threshold of each other and add connection with latency of 10ms
    print("Adding links between endpoints within adjacency threshold (to form links between seperate endpoint sources)...")
    added_links_list = []
    for source_id, source_data in t2t_dict.items():
        for dest_id, dest_data in t2t_dict.items():
            if (source_id == dest_id) or ((dest_id in source_data) and (source_id in dest_data)): # Skip if source and dest are the same or if the latency is already present
                continue
            if calc_coord_distance(source_data['coordinates'], dest_data['coordinates']) <= adjacency_threshold // 4: # tightening the restriction to be confident that the endpoints are adjacent
                if dest_id not in source_data:
                    source_data[dest_id] = 10
                if source_id not in dest_data:
                    dest_data[source_id] = 10
                added_links_list.append((id_to_endpoint_dict[source_id], id_to_endpoint_dict[dest_id]))
    print(f"Added {len(added_links_list)} links with 10ms latency.")
    #print(added_links_list)
    # Now add in ground station to endpoint latencies
    gs_name_list = gateway_dict.keys()
    # Give ground stations IDs and add to mappings
    next_id = int(max(id_to_endpoint_dict.keys())) + 1
    for gs_name in gs_name_list:
        gs_id = next_id
        next_id += 1
        endpoint_to_id_dict[gs_name] = gs_id
        id_to_endpoint_dict[gs_id] = gs_name
    
    adjacency_pair_dict = find_adjacency_pairs(aggr_endpoint_dict, gateway_dict, None, adjacency_threshold) # Find adjacency pairs between endpoints and ground stations
    gs_non_latency_keys = ['coordinates', 'description', 'data_source']
    for gs_name, data in adjacency_pair_dict.items():
        gs_id = int(endpoint_to_id_dict[gs_name])
        adj_endpoint, _ = data # Get the ground station's adjacent endpoint (don't need distance since we're adding a flat 10ms latency)
        adj_endpoint_id = int(endpoint_to_id_dict[adj_endpoint])
        if gs_id not in t2t_dict:
            t2t_dict[gs_id] = {'type': 'gateway', 'name': gs_name} # Add ground station to t2t_dict with base entries
            for non_latency_key in gs_non_latency_keys:
                if non_latency_key in gateway_dict[gs_name]: # If any other non-latency keys are present in the ground station data, add them to the t2t_dict
                    t2t_dict[gs_id][non_latency_key] = gateway_dict[gs_name][non_latency_key]
        for key, value in t2t_dict[adj_endpoint_id].items(): # For each latency entry for the given adjacent endpoint
            if key in endpoint_non_latency_keys or key in t2t_dict[gs_id]: # Skip non-latency keys and entries already present in the ground station's entry
                continue
            # assume all that's left are latency values
            dest_id = key
            dest_latency = value
            if type(dest_latency) is str and dest_latency.isnumeric(): # check for latencies that are strings (shouldn't be, but just in case)
                dest_latency = float(dest_latency) # Convert to float
            gs_latency = float(dest_latency) + 10 # assume 10ms latency from ground station to endpoint
            t2t_dict[gs_id][dest_id] = gs_latency
            #t2t_dict[dest_endpoint][gs_name] = gs_latency # sacrifice memory for lookup performance by adding in reverse direction

    # Write the t2t_dict to file
    if t2t_dict_output_filename is not None:
        write_json_file(t2t_dict_output_filename, t2t_dict)
    return t2t_dict

def genT2tTopology(t2t_dict):
    t2t_topology_dict = {}
    admin_key_list = ['type', 'name', 'coordinates', 'friendly_name', 'data_source', 'wp_index', 'description', 'locationName']
    for source_id, data_dict in t2t_dict.items(): # Loop through all entries in the t2t_dict and generate dictionary of latencies
        for key, value in data_dict.items():
            if key not in admin_key_list: # Skip administrative keys
                dest_id = key
                latency = value
                if source_id not in t2t_topology_dict:
                    t2t_topology_dict[source_id] = {}
                t2t_topology_dict[source_id][dest_id] = latency
                if dest_id not in t2t_topology_dict:
                    t2t_topology_dict[dest_id] = {}
                t2t_topology_dict[dest_id][source_id] = latency

    # Write topology to disk
    with open(topology_filename, 'w') as file:
        for source_endpoint, latency_dict in t2t_topology_dict.items():
            for dest_endpoint, latency in latency_dict.items():
                file.write(f"{source_endpoint},{dest_endpoint},{latency},{terrestrial_link_bandwidth}\n")
    print(f"Topology file written to {topology_filename}")
    
    return t2t_topology_dict

def add_gateway_gs(ground_station_dict_list, t2t_dict):
    highest_gid = 0
    for ground_station_dict in ground_station_dict_list:
        gid = ground_station_dict['gid']
        if gid > highest_gid:
            highest_gid = gid
    #next_gid = highest_gid + 1
    #for id in range(len(t2t_dict.keys())):
    next_gid = highest_gid + 1
    for id in t2t_dict.keys():
        if t2t_dict[id]['type'] == 'gateway':
            gateway_gid = next_gid
            next_gid += 1
            ground_station_entry_dict = {
                "gid": gateway_gid,
                "name": t2t_dict[id]['name'],
                "latitude_degrees_str": str(t2t_dict[id]['coordinates'][0]),
                "longitude_degrees_str": str(t2t_dict[id]['coordinates'][1]),
                "elevation_m_float": 0.0, # Does not appear to be used at this time
                "cartesian_x": 0.0, # cartesian coordinates do not appear to be used at this time
                "cartesian_y": 0.0,
                "cartesian_z": 0.0,
                "type": 9, # Using value of 9 to indicate gateway
                "next_update": "",
                "sat_re_LAC": -1
            }
            ground_station_dict_list.append(ground_station_entry_dict)
            # Now add GID to t2t_dict entry for this Gateway
            t2t_dict[id]['gid'] = gateway_gid
    # After assigning GID's to all Gateways, add GID's to all Endpoints
    for id in t2t_dict.keys():
        if 'type' in t2t_dict[id] and t2t_dict[id]['type'] == 'endpoint':
            t2t_dict[id]['gid'] = next_gid 
            next_gid += 1
    return ground_station_dict_list, t2t_dict

def load_t2t_dict(t2t_settings):
    import os
    t2t_dict = None
    # Check if t2t_dict file exists, if so load it
    if t2t_settings['t2t_dict_output_file'] != None and os.path.exists(t2t_settings['t2t_dict_output_file']): 
        t2t_dict = read_json_file(t2t_settings['t2t_dict_output_file'])
        if t2t_dict != None:
            print("..........(load_t2t_dict) Loaded t2t_dict from file.")
            return t2t_dict
    # No pre-existing dictionary; Generate t2t dictionary
    print("..........(load_t2t_dict) Generating t2t_dict...")
    endpoint_dict_list = []
    # Generate endpoint dictionary
    if t2t_settings['use_azure']:
        print("..........(load_t2t_dict) Loading Azure data center dictionary...")
        azure_latency_url = t2t_settings['azure_endpoint_latency_url']
        azure_location_file = t2t_settings['azure_endpoint_location_file']
        azure_dict_output_file = t2t_settings['azure_dict_output_file']
        if azure_dict_output_file == "": azure_dict_output_file = None
        azure_data_center_dict = load_azure_data_center_dict(None, azure_latency_url, azure_location_file, azure_dict_output_file)
        endpoint_dict_list.append(azure_data_center_dict)
    # Generate endpoint dictionary
    if t2t_settings['use_wonderproxy']:
        print("..........(load_t2t_dict) Loading WonderProxy server dictionary...")
        wonderproxy_endpoint_location_file = t2t_settings['wonderproxy_endpoint_location_file']
        wonderproxy_endpoint_latency_file = t2t_settings['wonderproxy_endpoint_latency_file']
        wonderproxy_dict_output_file = t2t_settings['wonderproxy_dict_output_file']
        wonderproxy_dict = load_wonderproxy_server_dict(wonderproxy_server_csv_filename=wonderproxy_endpoint_location_file, wonderproxy_latency_csv_filename=wonderproxy_endpoint_latency_file, wonderproxy_dict_filename=wonderproxy_dict_output_file)
        endpoint_dict_list.append(wonderproxy_dict)
    # Generate gateway dictionary
    print("..........(load_t2t_dict) Loading gateway dictionary...")
    gateway_dict = load_gateways_from_local_kml(t2t_settings['gateway_kmz_path'])
    # Generate t2t dictionary from endpoint and gateway dictionaries
    print("..........(load_t2t_dict) Generating t2t_dict from endpoint and gateway dictionaries...")
    t2t_dict = genT2tDict(gateway_dict, endpoint_dict_list, t2t_settings['t2t_dict_output_file'])

    return t2t_dict

def get_t2t_settings(main_config, output_filepath):
    #t2t_settings = {'use_t2t': True}
    t2t_settings = {}
    t2t_output_path = output_filepath+"/t2t/"
    t2t_settings['t2t_output_path'] = t2t_output_path # Save output path seperately in case it is needed for other functions
    t2t_settings['kmz_type'] = main_config["t2t_gateway_kmz_type"] # local or link
    t2t_settings['gateway_kmz_path'] = main_config["t2t_gateway_kmz_path"] # either a file path or a url
    if "t2t_dict_output_file" in main_config and main_config["t2t_dict_output_file"] != "":
        t2t_settings['t2t_dict_output_file'] = main_config["t2t_dict_output_file"] 
    else:
        t2t_settings['t2t_dict_output_file'] = None
    if "t2t_use_azure" in main_config:
        t2t_settings['use_azure'] = main_config["t2t_use_azure"]
        if t2t_settings['use_azure']:
            t2t_settings['azure_endpoint_location_file'] = main_config["t2t_azure_endpoint_location_file"]
            t2t_settings['azure_endpoint_latency_url'] = main_config["t2t_azure_endpoint_latency_url"]
            if "t2t_azure_dict_output_file" in main_config and main_config["t2t_azure_dict_output_file"] != "":
                t2t_settings['azure_dict_output_file'] = main_config["t2t_azure_dict_output_file"]
        else:
            t2t_settings['azure_dict_output_file'] = None
    else:
        t2t_settings['use_azure'] = False
    if "t2t_use_wonderproxy" in main_config:
        t2t_settings['use_wonderproxy'] = main_config["t2t_use_wonderproxy"]
        if t2t_settings['use_wonderproxy']:
            t2t_settings['wonderproxy_endpoint_location_file'] = main_config["t2t_wonderproxy_endpoint_location_file"]
            t2t_settings['wonderproxy_endpoint_latency_file'] = main_config["t2t_wonderproxy_endpoint_latency_file"]
            if "t2t_wonderproxy_dict_output_file" in main_config and main_config["t2t_wonderproxy_dict_output_file"] != "":
                t2t_settings['wonderproxy_dict_output_file'] = main_config["t2t_wonderproxy_dict_output_file"]
            else:
                t2t_settings['wonderproxy_dict_output_file'] = None
    else:
        t2t_settings['use_wonderproxy'] = False
    return t2t_settings

def add_t2t_links_to_connectivity_matrix(connectivity_matrix, links_characteristics, satellites_by_index, ground_station_dict_list, t2t_dict):
    # Expand the size of the connectivity matrix to include Endpoint nodes (gateways are already included)
    # Since links between gateways and endpoints are not time dependent, calculate only once
    # Find the number of entries in t2t_dict that have a 'type' of 'endpoint'
    if 'num_endpoints' in t2t_dict:
        num_endpoints = t2t_dict['num_endpoints']
    else:
        num_endpoints = len([key for key, value in t2t_dict.items() if value['type'] == 'endpoint'])
        t2t_dict['num_endpoints'] = num_endpoints
    current_conn_mat_size = len(connectivity_matrix)
    new_conn_mat_size = current_conn_mat_size + num_endpoints

    latency_matrix = links_characteristics['latency_matrix']
    throughput_matrix = links_characteristics['throughput_matrix']
    # Expand the existing rows of the matracies, initializing the new rows with 0s (or 0.0)
    for row in connectivity_matrix:
        row.extend([0] * num_endpoints)
    for row in latency_matrix:
        row.extend([0.0] * num_endpoints)
    for row in throughput_matrix:
        row.extend([0.0] * num_endpoints)
    # Add new rows to the matrices, initializing them with 0s (or 0.0)
    for _ in range(num_endpoints):
        connectivity_matrix.append([0] * new_conn_mat_size)
        latency_matrix.append([0.0] * new_conn_mat_size)
        throughput_matrix.append([0.0] * new_conn_mat_size)

    if 't2t_gw_to_ep_link_list' in t2t_dict: # Values already calculated, so just reference and update matrices
        for tuple in t2t_dict['t2t_gw_to_ep_link_list']:
            x, y = tuple
            connectivity_matrix[x][y] = 1 # Add bi-directional link to connectivity matrix
            connectivity_matrix[y][x] = 1
            source_id = t2t_dict['conn_mat_index_to_t2t_index'][x]
            dest_id = t2t_dict['conn_mat_index_to_t2t_index'][y]
            latency_matrix[x][y] = t2t_dict[source_id][dest_id] # Add latency value to latency matrix
            throughput_matrix[x][y] = terrestrial_link_bandwidth # Add throughput value to throughput matrix
    else: # Calculate values and update matrices
        t2t_dict['t2t_gw_to_ep_link_list'] = []
        t2t_dict['conn_mat_index_to_t2t_index'] = {}
        num_satellites = len(satellites_by_index)
        # Find next available GID along with number of gateways
        num_gateways = 0
        highest_gid = 0
        #gid_list = [] # TESTING:  Ensure assigned GID values are unique and incrementing sequentially
        for gs in ground_station_dict_list:
            #gid_list.append(int(gs['gid']))
            if gs['gid'] > highest_gid:
                    highest_gid = gs['gid']
            if gs['type'] == 9:
                num_gateways += 1
        next_gid = highest_gid + 1
        #gid_list.sort()
        #missing_gid = False
        #for i in range(len(gid_list)):
            #if gid_list[i] != i:
                #print(f"Missing GID value: {i}")
                #missing_gid = True
        #if not missing_gid:
            #print(f"Assigned GID values are unique and incrementing sequentially with max value of {highest_gid}.")
        # Loop through all Gateways and identify links to Endpoints; add these links to the connectivity matrix
        #for source_id, value in t2t_dict.items():
        for key in t2t_dict.keys():
            if type(t2t_dict[key]) is not dict or 'type' not in t2t_dict[key]: # Skip non-node entries
                continue
            source_id = key
            node_dict = t2t_dict[source_id]
            if node_dict['type'] == 'gateway':
                gateway_gid = node_dict['gid']
                for dest_id, _ in node_dict.items():
                    if dest_id in t2t_dict and t2t_dict[dest_id]['type'] == 'endpoint': # Verify destination entry is a node and that node is an endpoint
                        # Check if Endpoint has a GID; assign one if not
                        if 'gid' in t2t_dict[dest_id]:
                            endpoint_gid = t2t_dict[dest_id]['gid']
                        else:
                            t2t_dict[dest_id]['gid'] = next_gid
                            endpoint_gid = next_gid
                            next_gid += 1
                        x = num_satellites + gateway_gid
                        y = num_satellites + endpoint_gid
                        # Update matrices
                        if x >= new_conn_mat_size or y >= new_conn_mat_size:
                            print(f"Error: x or y value exceeds size of matrices: {x}, {y}")
                            print(f"  Old matrix size: {current_conn_mat_size} x {current_conn_mat_size}")
                            print(f"  New matrix size: {new_conn_mat_size} x {new_conn_mat_size}")
                            print(f"  Number of satellites: {len(satellites_by_index)}, number of existing ground stations: {len(ground_station_dict_list)}, number of Endpoints being added: {num_endpoints}")
                            print(f"  Highest GID value: {highest_gid}, GID value trying to be used: {endpoint_gid}")
                            print(f"  Number of Endpoints added: {len(t2t_dict['t2t_gw_to_ep_link_list'])}")
                            exit()
                        connectivity_matrix[x][y] = 1 # Add bi-directional link to connectivity matrix
                        connectivity_matrix[y][x] = 1
                        latency_matrix[x][y] = t2t_dict[source_id][dest_id] # Add latency value to latency matrix
                        latency_matrix[y][x] = t2t_dict[source_id][dest_id]
                        throughput_matrix[x][y] = terrestrial_link_bandwidth # Add throughput value to throughput matrix (currently a fixed value)
                        throughput_matrix[y][x] = terrestrial_link_bandwidth

                        # Record link in t2t_dict for reference in future time increments
                        t2t_dict['t2t_gw_to_ep_link_list'].append((x, y))
                        # Record index mappings for future reference
                        if x not in t2t_dict['conn_mat_index_to_t2t_index']:
                            t2t_dict['conn_mat_index_to_t2t_index'][x] = source_id
                        if y not in t2t_dict['conn_mat_index_to_t2t_index']:
                            t2t_dict['conn_mat_index_to_t2t_index'][y] = dest_id
    links_characteristics = {'latency_matrix': latency_matrix, 'throughput_matrix': throughput_matrix, 'distance_matrix': links_characteristics['throughput_matrix']}
    return connectivity_matrix, links_characteristics, t2t_dict
    
def plot_nodes_links_ground_stations(satellites_by_index, ground_stations, connectivity_matrix, operator_name, time_stamp, topology = None, routes = None, t2t_dict = None):
    if topology is not None and routes is not None:
        print("Both topology and routes provided. Please provide only one or none. Exiting.")
        return
    import matplotlib.pyplot as plt
    from mpl_toolkits.basemap import Basemap
    plt.figure(figsize=(10, 5))
    # Plot basemap
    m = Basemap(projection='cyl', llcrnrlat=-60, urcrnrlat=60, llcrnrlon=-180, urcrnrlon=180, resolution='c')
    m.drawcoastlines()
    m.drawcountries()
    m.fillcontinents(color='lightgray')
    gs_types = {0 : 'Customer Terminal', 9: 'Gateway', 1: 'Internet Endpoint'}
    gs_colors = {'Customer Terminal': 'blue', 'Gateway': 'green', 'Internet Endpoint': 'orange'}
    gs_marker = {'Customer Terminal': 'D', 'Gateway': '^', 'Internet Endpoint': 's'}
    gs_marker_size = {'Customer Terminal': 5, 'Gateway': 4, 'Internet Endpoint': 3}
    gs_alias_prefix = {'Customer Terminal': 'CT-', 'Gateway': 'GW-', 'Internet Endpoint': 'IE-'}
    plot_edges = False
    for ground_station in ground_stations:
        x, y = m(float(ground_station["longitude_degrees_str"]), float(ground_station["latitude_degrees_str"]))
        station_type = gs_types[ground_station["type"]]
        m.plot(x, y, marker=gs_marker[station_type], 
               markersize=gs_marker_size[station_type], 
               color = gs_colors[station_type], 
               latlon=True, alpha=0.5)#, label=station_type)#markersize=10)
        plt.annotate(gs_alias_prefix[station_type]+str(ground_station['gid']), 
                     xy=(x, y), 
                     textcoords="offset points", 
                     fontsize=5,
                     alpha=0.5,
                     xytext=(5,5), 
                     ha='center')
    
    for key in t2t_dict.keys():
        #if type(t2t_dict[key]) is not dict or "gid" not in t2t_dict[key]:
        #    continue
        if "type" in t2t_dict[key] and t2t_dict[key]["type"] == "endpoint":
            coordinates = t2t_dict[key]["coordinates"]
            #print(f"Endpoint {t2t_dict[key]['name']} at coordinates {coordinates}")
            x, y = m(coordinates[1], coordinates[0])
            station_type = 'Internet Endpoint'
            m.plot(x, y, marker=gs_marker[station_type], 
                   markersize=gs_marker_size[station_type], 
                   color = gs_colors[station_type], 
                   latlon=True, 
                   alpha=0.5)#, label=station_type)
            plt.annotate(gs_alias_prefix[station_type]+str(t2t_dict[key]['gid']), 
                         xy=(x, y), 
                         textcoords="offset points",
                         fontsize=5,
                         alpha=0.5, 
                         xytext=(5,5), 
                         ha='center')
    if plot_edges:
        pass
    # Plot information
    if time_stamp is None:
        title_string = f"{operator_name.capitalize()} Ground Station Topology"
    else:
        title_string = f"{operator_name.capitalize()} Ground Station Topology at {time_stamp}"
    plt.title(title_string)
    plt.xlabel("Longitude")
    plt.ylabel("Latitude")
    # Plot legend
    for gs_type in gs_types:
        station_type = gs_types[gs_type]
        plt.plot([], [], marker=gs_marker[station_type], markersize=gs_marker_size[station_type], color = gs_colors[station_type], label=station_type, linewidth=0)
    plt.legend(loc='lower left', title='Ground Station Types')
    plt.show()

if __name__ == '__main__':
    azure_data_center_dict_with_latencies = load_azure_data_center_dict(azure_latency_json_filename = azure_latency_json_filename, azure_url = url, csv_filename = azure_csv_filename)
    #azure_data_center_latency_dict = load_azure_data_center_latency(azure_latency_json_filename)
    #wonderproxy_server_dict = load_wonderproxy_server_location_coordinates(wonderproxy_server_csv_filename)
    #wonderproxy_server_latency_dict = load_wonderproxy_latency_dict(wonderproxy_json_filename = wonderproxy_latency_json_filename, csv_filename = wonderproxy_latency_csv_filename, wonderproxy_server_dict = wonderproxy_server_dict)
    wonderproxy_server_dict_with_latencies = load_wonderproxy_server_dict(wonderproxy_server_csv_filename, wonderproxy_latency_json_filename, wonderproxy_latency_csv_filename)
    ground_station_dict = load_gateways_from_local_kml(kmz_file)
    t2t_dict = genT2tDict(ground_station_dict, [azure_data_center_dict_with_latencies, wonderproxy_server_dict_with_latencies], topology_filename)
    print(f"{len(t2t_dict.keys())} endpoints in t2t topology loaded.")
    t2t_topology_dict = genT2tTopology(t2t_dict)
    print(f"{len(t2t_topology_dict.keys())} endpoints in t2t topology written to file.")
    
    exit()
    print(f"{len(ground_station_dict)} ground stations loaded.")
    #print("Ground station locations:")
    #pprint.pprint(ground_station_dict)
    azure_adjacency_pair_dict = find_adjacency_pairs(azure_data_centers, ground_station_dict, None)
    print(f"{len(azure_adjacency_pair_dict.keys())} Azure adjacency pairs:")
    #for adjacency_pair in azure_adjacency_pair_list:
    #    print(adjacency_pair)
    gs_without_adjacency = []
    for item in ground_station_dict.keys():
        if item in azure_adjacency_pair_dict:
            continue
        else:
            gs_without_adjacency.append(item)
    print(f"{len(gs_without_adjacency)} ground stations without Azure adjacency:")
    #for item in gs_without_adjacency:
    #    print(item)
    no_adjacency_gs_dict = {location: ground_station_dict[location] for location in gs_without_adjacency}
    
    print(f"{len(wonderproxy_server_dict)} WonderProxy server locations loaded.")
    #import pprint
    #print("WonderProxy server locations:")
    #pprint.pprint(wonderproxy_server_dict)
    wonderproxy_adjacency_pair_dict = find_adjacency_pairs(wonderproxy_server_dict, no_adjacency_gs_dict)
    print(f"{len(wonderproxy_adjacency_pair_dict.keys())} WonderProxy adjacency pairs:")
    #for gs_name, data in wonderproxy_adjacency_pair_dict.items():
    #    print(f"{gs_name}: {data}")
    gs_without_adjacency = []
    for item in no_adjacency_gs_dict.keys():
        if item in wonderproxy_adjacency_pair_dict.keys():
            continue
        else:
            gs_without_adjacency.append(item)
    print(f"{len(gs_without_adjacency)} ground stations without adjacency:")
    for item in gs_without_adjacency:
        print(item)
    #import pprint
    #pprint.pprint(azure_data_centers)
# Azure variables
azure_latency_json_filename = 'gs_files/azure_latency_data.json'
azure_csv_filename = 'gs_files/AzureDataCenterLocations.csv'
url = 'https://learn.microsoft.com/en-us/azure/networking/azure-network-latency'
sectionIdList = ['tabpanel_2_WestUS_Americas', 'tabpanel_2_EastUS_Americas', 'tabpanel_2_CentralUS_Americas', 'tabpanel_2_Canada_Americas', 'tabpanel_2_Australia_APAC', 'tabpanel_2_Japan_APAC', 'tabpanel_2_WesternEurope_Europe', 'tabpanel_2_CentralEurope_Europe', 'tabpanel_2_NorwaySweden_Europe', 'tabpanel_2_UKNorthEurope_Europe', 'tabpanel_2_Korea_APAC', 'tabpanel_2_India_APAC', 'tabpanel_2_Asia_APAC', 'tabpanel_2_israel-qatar-uae_MiddleEast', 'tabpanel_2_southafrica_MiddleEast']
# Wonderproxy variables
# WonderProxy files downloaded from: https://wonderproxy.com/blog/a-day-in-the-life-of-the-internet/
wonderproxy_server_csv_filename = 'gs_files/wonderproxy_servers-2020-07-19.csv'
wonderproxy_latency_csv_filename = 'gs_files/wonderproxy_pings-2020-07-19-2020-07-20.csv'
wonderproxy_latency_json_filename = 'gs_files/wonderproxy_latency_data.json'
# Unofficial Starlink Global Gateways and PoPs KMZ file variables
kmz_file = "gs_files/UnofficialStarlinkGlobalGatewaysNPoPs_noLinks.kmz"
usable_bands = ['Ka', 'E']

adjacency_threshold = 2000 # Threshold distance in kilometers for considering two points as adjacent (be within 5ms latency (10ms round-trip))
topology_filename = 'gs_files/t2t_topology.txt'
# >>>>>>>>> WonderProxy Functions <<<<<<<<<<<<
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
            wonderproxy_server_dict[locationName] = {'coordinates': (latitude, longitude), 'index': index}
    return wonderproxy_server_dict

def gen_wonderproxy_server_latency_json_from_csv(csv_filename, json_filename, wonderproxy_server_dict):
    if csv_filename is None or json_filename is None or wonderproxy_server_dict is None:
        print(f"Received: csv_filename: {True if csv_filename else False}, json_filename: {True if json_filename else False}, wonderproxy_server_dict: {True if wonderproxy_server_dict else False}")
        raise ValueError("csv_filename, json_filename, and server dict must be provided.")
    # Build wonderproxy_server_dict lookup table by index
    wonderproxy_server_lookup_dict = {}
    for locationName, data in wonderproxy_server_dict.items():
        index = data['index']
        wonderproxy_server_lookup_dict[index] = locationName
    wonderproxy_latency_data = read_csv_file(csv_filename)
    if wonderproxy_latency_data is None:
        raise ValueError(f"Failed to read data from {csv_filename}")
    wonderproxy_latency_dict = {}
    for row in wonderproxy_latency_data[1:]:
        source_index = int(row[0])
        if source_index not in wonderproxy_server_lookup_dict:
            print(f"Source index {source_index} not found in server dictionary.")
            continue
        dest_index = int(row[1])
        if dest_index not in wonderproxy_server_lookup_dict:
            print(f"Destination index {dest_index} not found in server dictionary.")
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
    write_json_file(json_filename, wonderproxy_latency_dict)
    print(f"Data saved to {json_filename}")
    return wonderproxy_latency_dict

def load_wonderproxy_latency_dict(wonderproxy_json_filename = None, csv_filename = None, wonderproxy_server_dict = None):
    wonderproxy_latency_dict = None
    if wonderproxy_json_filename is None and csv_filename is None:
        raise ValueError("Either csv_filename or json_filename must be provided.")
    if wonderproxy_json_filename is not None:
        wonderproxy_latency_dict = read_json_file(wonderproxy_json_filename)
    if wonderproxy_latency_dict is None and csv_filename is None:
        raise ValueError("Cached data not found and csv_filename is not provided.")
    elif wonderproxy_latency_dict is None:
        wonderproxy_latency_dict = gen_wonderproxy_server_latency_json_from_csv(csv_filename, wonderproxy_json_filename, wonderproxy_server_dict)
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
    with open(csv_filename, 'r', encoding=encoding) as file:
        reader = csv.reader(file)
        next(reader) # Skip the header row
        for row in reader:
            regionName = row[0]
            locationName = row[1]
            latitude = float(row[2])
            longitude = float(row[3])
            if regionName in azure_data_centers:
                azure_data_centers[regionName]['coordinates'] = (latitude, longitude)
                azure_data_centers[regionName]['locationName'] = locationName
            else:
                print(f"Region not found in Azure Data Center dictionary: {regionName}")
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
        if azure_latency_json_filename is not None:
            write_json_file(azure_latency_json_filename, azure_data_centers)
            print(f"Data saved to {azure_latency_json_filename}")
    return azure_data_centers

def load_azure_data_centers(azure_json_filename = None, azure_url = None, csv_filename = None):
    azure_data_centers = load_azure_data_center_latency(azure_json_filename, azure_url)
    add_azure_location_coordinates(azure_data_centers, csv_filename)
    # List any azure_data_centers that do not have coordinates
    for region, data in azure_data_centers.items():
        if 'coordinates' not in data:
            print(f"No coordinates found for region: {region}")
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

# >>>>>>>>> CSV/JSON Functions <<<<<<<<<<<<
# Function to read file as json and return the data as a dictionary
# If file is not present, return None
def read_json_file(file_path):
    import json
    try:
        with open(file_path, 'r') as f:
            data = json.load(f)
        return data
    except FileNotFoundError:
        return None
    
# Function to write data as json to a file
def write_json_file(file_path, data):
    import json
    with open(file_path, 'w') as f:
        json.dump(data, f)

def read_csv_file(file_path):
    import csv
    encoding = detect_text_file_encoding(file_path)
    if encoding is None:
        raise ValueError("Failed to detect the encoding of the csv file.")
    try:
        with open(file_path, 'r', encoding=encoding) as file:
            reader = csv.reader(file)
            data = [row for row in reader]
        return data
    except FileNotFoundError:
        return None

def detect_text_file_encoding(file_path):
    #import csv
    # List of encodings to try, in order of likelihood
    encodings_to_try = ['utf-8', 'cp1252', 'ISO-8859-1', 'latin1']
    for encoding in encodings_to_try:
        try:
            with open(file_path, mode='r', encoding=encoding) as file:
                # Read first line of text file
                first_line = file.readline()
                #reader = csv.reader(file)
                ## Try reading the first line to see if the encoding works
                #next(reader)
                print(f"File opened successfully with encoding: {encoding}")
                return encoding
        except UnicodeDecodeError:
            continue  # Try the next encoding if decoding failed
    else:
        print(f"Failed to open the file with any of the tried encodings ({encodings_to_try}).")
        return None

# >>>>>>>>> KML Functions <<<<<<<<<<<<
def parse_local_kml_file(file_path, print_content=True):
    import zipfile
    import xml.etree.ElementTree as ET
    """Parse a local KMZ file and print the contents of the KML file inside."""
    with zipfile.ZipFile(file_path, 'r') as kmz:
        # List all the files in the KMZ archive
        file_list = kmz.namelist()
        print(f"Files in KMZ archive: {file_list}")
        
        # Extract the KML file (assuming it's the first file in the archive)
        kml_filename = file_list[0]
        print(f"Extracting and parsing KML file: {kml_filename}")
        
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
                    'description': data['description']
                }
    return ground_station_dict

def load_groundstations_from_local_kml(kmz_file):
    kml_root = parse_local_kml_file(kmz_file, print_content=False)
    placemark_dict = extract_placemarks(kml_root)
    ground_station_dict = get_ground_stations_from_placemarks(placemark_dict)
    return ground_station_dict

def genT2tTopologyFile(gs_dict, endpoint_dict_list, endpoint_latency_dict_list, topology_filename):
    latency_matrix_dict = {} # build dictionary as 2D array of latencies between each pair of endpoints
    
    aggr_endpoint_dict = {} # build dictionary of all endpoints from all dictionaries
    for endpoint_dict in endpoint_dict_list:
        aggr_endpoint_dict.update(endpoint_dict)
    aggr_endpoint_latency_dict = {} # build dictionary of all latencies between each pair of endpoints
    for endpoint_latency_dict in endpoint_latency_dict_list:
        aggr_endpoint_latency_dict.update(endpoint_latency_dict)
    # build list of all endpoints from all dictionaries
    endpoint_name_list = aggr_endpoint_dict.keys()
    # Loop through endpoints to build latency matrix
    for source_endpoint in endpoint_name_list:
        if source_endpoint not in latency_matrix_dict:
            latency_matrix_dict[source_endpoint] = {}
        for dest_endpoint in endpoint_name_list:
            if source_endpoint == dest_endpoint:
                continue
            for endpoint_latency_dict in endpoint_latency_dict_list:
                if source_endpoint in endpoint_latency_dict:
                    if dest_endpoint in endpoint_latency_dict[source_endpoint]:
                        latency = endpoint_latency_dict[source_endpoint][dest_endpoint]
                        latency_matrix_dict[source_endpoint][dest_endpoint] = latency
                        if dest_endpoint not in latency_matrix_dict:
                            latency_matrix_dict[dest_endpoint] = {}
                        latency_matrix_dict[dest_endpoint][source_endpoint] = latency # sacrifice memory for lookup performance by adding in reverse direction
                        break
    # Now add in ground station to endpoint latencies
    adjacency_pair_dict = find_adjacency_pairs(aggr_endpoint_dict, gs_dict, None, adjacency_threshold)
    for gs_name, data in adjacency_pair_dict.items():
        adj_endpoint, _ = data
        adj_endpoint_latency_dict = latency_matrix_dict[adj_endpoint]
        if gs_name not in latency_matrix_dict:
            latency_matrix_dict[gs_name] = {}
        for dest_endpoint, dest_latency in adj_endpoint_latency_dict.items():
            gs_latency = float(dest_latency) + 10 # assume 10ms latency from ground station to endpoint
            latency_matrix_dict[gs_name][dest_endpoint] = gs_latency
            latency_matrix_dict[dest_endpoint][gs_name] = gs_latency # sacrifice memory for lookup performance by adding in reverse direction

    # Write the topology file
    with open(topology_filename, 'w') as file:
        for source_endpoint, latency_dict in latency_matrix_dict.items():
            for dest_endpoint, latency in latency_dict.items():
                file.write(f"{source_endpoint} {dest_endpoint} {latency}\n")
    print(f"Topology file written to {topology_filename}")
    return latency_matrix_dict

if __name__ == '__main__':
    azure_data_center_dict = load_azure_data_centers(azure_json_filename = azure_latency_json_filename, azure_url = url, csv_filename = azure_csv_filename)
    azure_data_center_latency_dict = load_azure_data_center_latency(azure_latency_json_filename)
    wonderproxy_server_dict = load_wonderproxy_server_location_coordinates(wonderproxy_server_csv_filename)
    wonderproxy_server_latency_dict = load_wonderproxy_latency_dict(wonderproxy_json_filename = wonderproxy_latency_json_filename, csv_filename = wonderproxy_latency_csv_filename, wonderproxy_server_dict = wonderproxy_server_dict)
    ground_station_dict = load_groundstations_from_local_kml(kmz_file)
    t2t_topology_dict = genT2tTopologyFile(ground_station_dict, [azure_data_center_dict, wonderproxy_server_dict], [azure_data_center_latency_dict, wonderproxy_server_latency_dict], topology_filename)
    import pprint
    print(f"{len(t2t_topology_dict.keys())} endpoints in t2t topology loaded.")
    pprint.pprint(t2t_topology_dict)
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
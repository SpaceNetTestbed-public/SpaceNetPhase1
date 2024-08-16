import requests
from bs4 import BeautifulSoup
from pprint import pprint
import json

json_filename = 'azure_latency_data.json'

url = 'https://learn.microsoft.com/en-us/azure/networking/azure-network-latency'
sectionIdList = ['tabpanel_2_WestUS_Americas', 'tabpanel_2_EastUS_Americas', 'tabpanel_2_CentralUS_Americas', 'tabpanel_2_Canada_Americas', 'tabpanel_2_Australia_APAC', 'tabpanel_2_Japan_APAC', 'tabpanel_2_WesternEurope_Europe', 'tabpanel_2_CentralEurope_Europe', 'tabpanel_2_NorwaySweden_Europe', 'tabpanel_2_UKNorthEurope_Europe', 'tabpanel_2_Korea_APAC', 'tabpanel_2_India_APAC', 'tabpanel_2_Asia_APAC', 'tabpanel_2_israel-qatar-uae_MiddleEast', 'tabpanel_2_southafrica_MiddleEast']

def scrape_azure_latency_data(url, section_id='tabpanel_2_WestUS_Americas'):
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

# Function to read file as json and return the data as a dictionary
# If file is not present, return None
def read_json_file(file_path):
    try:
        with open(file_path, 'r') as f:
            data = json.load(f)
        return data
    except FileNotFoundError:
        return None
    
# Function to write data as json to a file
def write_json_file(file_path, data):
    with open(file_path, 'w') as f:
        json.dump(data, f)

def main():
    latency_dict = read_json_file(json_filename)
    if latency_dict is None:
        print("No cached data found. Scraping data from Azure website.")
        latency_data = []
        for section_id in sectionIdList:
            latency_data += scrape_azure_latency_data(url, section_id)
        # Latency data is a list of dictionaries
        # Each dictionary has one key named 'Source' and two other keys named as Azure regions
        # Value of 'Source' key is the name of the source region
        # Value of the other two keys are the latency values in milliseconds
        # Consume list and build a dictionary with source region as key and a dictionary of destination regions and their latencies as value
        latency_dict = {}
        for data in latency_data:
            source = data['Source']
            del data['Source']
            if source not in latency_dict:
                latency_dict[source] = data
            else:
                latency_dict[source].update(data)
        # Save combined dictionary as a JSON file
        write_json_file(json_filename, latency_dict)
    else:
        print("Cached data found.")
    
    pprint(latency_dict)



if __name__ == '__main__':
    main()
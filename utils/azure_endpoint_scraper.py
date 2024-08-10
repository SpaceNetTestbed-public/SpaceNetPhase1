import requests
from bs4 import BeautifulSoup

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

def main():
    url = 'https://learn.microsoft.com/en-us/azure/networking/azure-network-latency'
    sectionIdList = ['tabpanel_2_WestUS_Americas', 'tabpanel_2_EastUS_Americas', 'tabpanel_2_CentralUS_Americas', 'tabpanel_2_Canada_Americas', 'tabpanel_2_Australia_APAC', 'tabpanel_2_Japan_APAC', 'tabpanel_2_WesternEurope_Europe', 'tabpanel_2_CentralEurope_Europe', 'tabpanel_2_NorwaySweden_Europe', 'tabpanel_2_UKNorthEurope_Europe', 'tabpanel_2_Korea_APAC', 'tabpanel_2_India_APAC', 'tabpanel_2_Asia_APAC', 'tabpanel_2_israel-qatar-uae_MiddleEast', 'tabpanel_2_southafrica_MiddleEast']
    latency_data = []
    for section_id in sectionIdList:
        latency_data += scrape_azure_latency_data(url, section_id)
    print(latency_data)


if __name__ == '__main__':
    main()
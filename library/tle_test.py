from mpl_toolkits.basemap import Basemap
from skyfield.api import load
import matplotlib.pyplot as plt
import numpy as np

# Load TLE file
tle_file = '/home/spacenet/simulator/dynamic-topology-generator/utils/starlink_tles/starlink_1724155200'
satellites = load.tle_file(tle_file)

# Load timescale
ts = load.timescale()
t = ts.now()

# Create a new figure
fig = plt.figure(figsize=(12, 8))

# Create a Basemap
m = Basemap(projection='ortho', lat_0=10, lon_0=-50, resolution='c')

# Draw coastlines, and the edges of the map.
m.drawcoastlines(color='white')
m.fillcontinents(color='white',lake_color=(200/255, 204/255, 208/255))
m.drawmapboundary(fill_color=(200/255, 204/255, 208/255), linewidth=0)  # Set linewidth to 0

# Draw meridians
meridians = np.arange(-180.,181.,45.)
m.drawmeridians(meridians, labels=[False,False,False,True], color='grey')

# Draw equator
parallels = np.arange(-90.,91.,90.)
m.drawparallels(parallels, labels=[True,False,False,False], color='grey')

# Loop over all satellites and plot their positions
sats_list = [sat.at(t).subpoint() for sat in satellites]
num_subarrays = len(sats_list) // 22
sats_per_orbit_list = np.array_split(sats_list, num_subarrays)
for orbit in sats_per_orbit_list:
    lon_lat = [(point.longitude.degrees, point.latitude.degrees) for point in orbit]
    lon_lat.append(lon_lat[0])
    x, y = m(*zip(*lon_lat))
    m.plot(x, y, '-*', markersize=5, color='k', lw=1)
    
plt.tight_layout()
plt.show()
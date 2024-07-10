import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np

# Initialize list for semi-major axis values
semi_major_axes = []
inclination     = []

# Open TLE file and read lines
with open('C:/Users/BluBoy/Desktop/Professional/Git/Repositories/dynamic-topology-generator/utils/starlink_tles/starlink_1719866022', 'r') as tle_file:
    Lines = tle_file.readlines()

# Process TLEs
for i in range(0, len(Lines), 3):
    tle_second_line = list([_f for _f in Lines[i+2].strip("\n").split(" ") if _f])
    tle_n = float(tle_second_line[7]) * 2 * np.pi / 86400
    tle_a = (396800. / (tle_n ** 2)) ** (1. / 3.) - 6378.135
    semi_major_axes.append((tle_a, tle_second_line[2]))

print(len(range(0, len(Lines), 3)))


filtered_data = []
for sma, inc in semi_major_axes:
    if sma > 525 and sma < 535:
        filtered_data.append(inc)

# Plot histogram of semi-major axes
ax = plt.gca()
ax.hist(filtered_data, bins=20, color='blue', alpha=0.7)
ax.xaxis.set_major_locator(MaxNLocator(nbins=6))
ax.set_xlabel('X')
ax.set_ylabel('Frequency')
plt.show()

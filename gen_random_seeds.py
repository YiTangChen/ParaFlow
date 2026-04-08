"""
Generate random seed points in the Atlantic Ocean.
Seeds are written as binary file of (x, y, z) double triples
in Cartesian coordinates on Earth's sphere (radius = 6371228 m).

Atlantic Ocean bounding box (approximate):
  Latitude:  -60 to 70 degrees
  Longitude: -80 to 20 degrees
"""

import struct
import math
import random
import os

dst_dir = "seeds"
os.makedirs(dst_dir, exist_ok=True)

# create other folders
folders = ["drawSubdomain", "png", "streamlines", "pathlines", "block_outputs"]
for folder in folders:
    os.makedirs(folder, exist_ok=True)

EARTH_RADIUS = 6371220.0  # meters, matching MPASO mesh sphere_radius
EARTH_RADIUS_1 = EARTH_RADIUS - 1.0  # meters

# Atlantic Ocean bounding box
# LAT_MIN = -60.0
# LAT_MAX =  70.0
# LON_MIN = -80.0
# LON_MAX =  20.0
LAT_MIN = -90.0
LAT_MAX =  90.0
LON_MIN = -180.0
LON_MAX =  180.0

N_SEEDS = 1000
SEED_FILE = os.path.join(dst_dir, "random_seeds_1000.bin")

random.seed(42)

def latlon_to_xyz(lat_deg, lon_deg):
    lat = math.radians(lat_deg)
    lon = math.radians(lon_deg)
    x = EARTH_RADIUS_1 * math.cos(lat) * math.cos(lon)
    y = EARTH_RADIUS_1 * math.cos(lat) * math.sin(lon)
    z = EARTH_RADIUS_1 * math.sin(lat)
    return x, y, z

with open(SEED_FILE, 'wb') as f:
    count = 0
    # Use rejection sampling on a uniform sphere within the bounding box.
    # For uniform distribution on the sphere, sample sin(lat) uniformly.
    sin_lat_min = math.sin(math.radians(LAT_MIN))
    sin_lat_max = math.sin(math.radians(LAT_MAX))

    while count < N_SEEDS:
        # Uniform sampling on sphere surface within lat/lon box
        sin_lat = random.uniform(sin_lat_min, sin_lat_max)
        lat = math.degrees(math.asin(sin_lat))
        lon = random.uniform(LON_MIN, LON_MAX)

        x, y, z = latlon_to_xyz(lat, lon)
        f.write(struct.pack('ddd', x, y, z))
        count += 1

print(f"Generated {count} seeds in '{SEED_FILE}'")
print(f"Lat range: [{LAT_MIN}, {LAT_MAX}], Lon range: [{LON_MIN}, {LON_MAX}]")
print(f"Earth radius: {EARTH_RADIUS} m")

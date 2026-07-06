"""
Generate random seed points in the ocean.
Seeds are written as binary file of (x, y, z) double triples
in Cartesian coordinates on Earth's sphere (radius = 6371228 m).

Atlantic Ocean bounding box (approximate):
  Latitude:  -60 to 70 degrees
  Longitude: -80 to 20 degrees

East Coast presets:
  eastcoast_surface  -- surface seeds near US East Coast / Gulf Stream
  eastcoast_depth    -- deep seeds same region, ~200-2000 m depth

Usage:
    python3 gen_random_seeds.py [n_seeds]
    python3 gen_random_seeds.py [n_seeds] --preset eastcoast_surface
    python3 gen_random_seeds.py [n_seeds] --preset eastcoast_depth
    python3 gen_random_seeds.py [n_seeds] --lat-min 25 --lat-max 45 --lon-min -80 --lon-max -60 --output seeds/custom.bin
"""

import argparse
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

DEFAULT_N_SEEDS = 10000

# East Coast / Gulf Stream region
EASTCOAST_LAT_MIN = 25.0
EASTCOAST_LAT_MAX = 45.0
EASTCOAST_LON_MIN = -80.0
EASTCOAST_LON_MAX = -60.0

PRESETS = {
    "eastcoast_surface": {
        "lat_min": EASTCOAST_LAT_MIN,
        "lat_max": EASTCOAST_LAT_MAX,
        "lon_min": EASTCOAST_LON_MIN,
        "lon_max": EASTCOAST_LON_MAX,
        "depth_min": 0.0,
        "depth_max": 10.0,
        "output": os.path.join(dst_dir, "eastcoast_surface_seeds.bin"),
    },
    "eastcoast_depth": {
        "lat_min": EASTCOAST_LAT_MIN,
        "lat_max": EASTCOAST_LAT_MAX,
        "lon_min": EASTCOAST_LON_MIN,
        "lon_max": EASTCOAST_LON_MAX,
        "depth_min": 200.0,
        "depth_max": 2000.0,
        "output": os.path.join(dst_dir, "eastcoast_depth_seeds.bin"),
    },
}


def latlon_to_xyz(lat_deg, lon_deg, depth=0.0):
    """Convert lat/lon/depth to Cartesian XYZ. depth is meters below surface (positive = deeper)."""
    r = EARTH_RADIUS_1 - depth
    lat = math.radians(lat_deg)
    lon = math.radians(lon_deg)
    x = r * math.cos(lat) * math.cos(lon)
    y = r * math.cos(lat) * math.sin(lon)
    z = r * math.sin(lat)
    return x, y, z


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "n_seeds",
        nargs="?",
        type=int,
        default=DEFAULT_N_SEEDS,
        help=f"number of random seeds to generate (default: {DEFAULT_N_SEEDS})",
    )
    parser.add_argument("--preset", choices=list(PRESETS.keys()), default=None,
                        help="use a named preset (overrides bounding box defaults)")
    parser.add_argument("--lat-min", type=float, default=None)
    parser.add_argument("--lat-max", type=float, default=None)
    parser.add_argument("--lon-min", type=float, default=None)
    parser.add_argument("--lon-max", type=float, default=None)
    parser.add_argument("--depth-min", type=float, default=None,
                        help="minimum depth in meters below surface (default: 0)")
    parser.add_argument("--depth-max", type=float, default=None,
                        help="maximum depth in meters below surface (default: 0)")
    parser.add_argument("--output", type=str, default=None,
                        help="output file path (default: seeds/random_seeds.bin)")
    return parser.parse_args()


def main():
    args = parse_args()
    if args.n_seeds <= 0:
        raise ValueError("n_seeds must be a positive integer")

    n_seeds = args.n_seeds

    # resolve settings: preset < explicit args
    lat_min   = LAT_MIN
    lat_max   = LAT_MAX
    lon_min   = LON_MIN
    lon_max   = LON_MAX
    depth_min = 0.0
    depth_max = 0.0
    seed_file = os.path.join(dst_dir, "random_seeds.bin")

    if args.preset is not None:
        p = PRESETS[args.preset]
        lat_min   = p["lat_min"]
        lat_max   = p["lat_max"]
        lon_min   = p["lon_min"]
        lon_max   = p["lon_max"]
        depth_min = p["depth_min"]
        depth_max = p["depth_max"]
        seed_file = p["output"]

    if args.lat_min  is not None: lat_min   = args.lat_min
    if args.lat_max  is not None: lat_max   = args.lat_max
    if args.lon_min  is not None: lon_min   = args.lon_min
    if args.lon_max  is not None: lon_max   = args.lon_max
    if args.depth_min is not None: depth_min = args.depth_min
    if args.depth_max is not None: depth_max = args.depth_max
    if args.output   is not None: seed_file = args.output

    random.seed(42)

    with open(seed_file, "wb") as f:
        count = 0
        sin_lat_min = math.sin(math.radians(lat_min))
        sin_lat_max = math.sin(math.radians(lat_max))

        while count < n_seeds:
            sin_lat = random.uniform(sin_lat_min, sin_lat_max)
            lat = math.degrees(math.asin(sin_lat))
            lon = random.uniform(lon_min, lon_max)
            depth = random.uniform(depth_min, depth_max)

            x, y, z = latlon_to_xyz(lat, lon, depth)
            f.write(struct.pack("ddd", x, y, z))
            count += 1

    print(f"Generated {count} seeds in '{seed_file}'")
    print(f"Lat range: [{lat_min}, {lat_max}], Lon range: [{lon_min}, {lon_max}]")
    print(f"Depth range: [{depth_min}, {depth_max}] m")
    print(f"Earth radius: {EARTH_RADIUS} m")


if __name__ == "__main__":
    main()

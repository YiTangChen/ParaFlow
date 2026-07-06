"""
Check that atlantic_seeds.bin are correctly placed in the Atlantic Ocean.

Reads the combined subdomain binary (all blocks) as an ocean mask,
then overlays the atlantic seeds to verify:
  1. Seeds are within lat/lon bounds
  2. Seeds land on MPASO ocean cells (not land)
  3. Visual check on map

Usage:
  python check_atlantic_seeds.py [draw_dir] [seed_file]
  python check_atlantic_seeds.py drawSubdomain atlantic_seeds.bin
"""

import struct
import math
import sys
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

LAT_RES = 600
LON_RES = 1200

def read_ocean_mask(draw_dir):
    """Read combined subdomain binary as ocean mask (value > 0 = ocean cell)."""
    path = os.path.join(draw_dir, f"{LAT_RES}_{LON_RES}.bin")
    data = np.fromfile(path, dtype=np.float64).reshape((LAT_RES, LON_RES))
    return data  # value > 0 means at least one block covers this cell

def read_seeds_bin(seed_file):
    """Read XYZ binary seed file, return list of (x, y, z)."""
    seeds = []
    with open(seed_file, 'rb') as f:
        while True:
            raw = f.read(24)
            if len(raw) < 24:
                break
            x, y, z = struct.unpack('ddd', raw)
            seeds.append((x, y, z))
    return seeds

def xyz_to_latlon_deg(x, y, z):
    lon = math.degrees(math.atan2(y, x))
    lat = math.degrees(math.atan2(z, math.hypot(x, y)))
    return lat, lon

def latlon_to_idx(lat_deg, lon_deg):
    """Convert lat/lon (degrees) to (latIdx, lonIdx) matching DrawSubdomain.cpp."""
    lat = math.radians(lat_deg)
    lon = math.radians(lon_deg)
    lon = math.fmod(lon, 2 * math.pi)
    if lon < 0:
        lon += 2 * math.pi
    lat_idx = int(round((lat + math.pi / 2) / math.pi * LAT_RES - 0.5))
    lon_idx = int(round(lon / (2 * math.pi) * LON_RES - 0.5))
    # clamp to valid range
    lat_idx = max(0, min(LAT_RES - 1, lat_idx))
    lon_idx = max(0, min(LON_RES - 1, lon_idx))
    return lat_idx, lon_idx

def main():
    draw_dir  = sys.argv[1] if len(sys.argv) > 1 else 'drawSubdomain'
    seed_file = sys.argv[2] if len(sys.argv) > 2 else 'atlantic_seeds.bin'

    print(f"Ocean mask : {draw_dir}/{LAT_RES}_{LON_RES}.bin")
    print(f"Seed file  : {seed_file}")

    ocean_mask = read_ocean_mask(draw_dir)
    seeds = read_seeds_bin(seed_file)
    print(f"Loaded {len(seeds)} seeds\n")

    # Convert seeds to lat/lon and check
    lats, lons = [], []
    on_land   = []  # seeds with ocean_mask value = 0
    out_bbox  = []  # seeds outside the nominal Atlantic bbox

    LAT_MIN, LAT_MAX = -60.0,  70.0
    LON_MIN, LON_MAX = -80.0,  20.0

    for i, (x, y, z) in enumerate(seeds):
        lat, lon = xyz_to_latlon_deg(x, y, z)
        lats.append(lat)
        lons.append(lon)

        lat_idx, lon_idx = latlon_to_idx(lat, lon)
        ocean_val = ocean_mask[lat_idx, lon_idx]

        if ocean_val == 0.0:
            on_land.append((i, lat, lon))
        if not (LAT_MIN <= lat <= LAT_MAX and LON_MIN <= lon <= LON_MAX):
            out_bbox.append((i, lat, lon))

    # Summary
    print(f"Seeds on land (ocean_mask=0) : {len(on_land)}")
    if on_land:
        for idx, lat, lon in on_land[:10]:
            print(f"  seed {idx:4d}: lat={lat:.2f}, lon={lon:.2f}")
        if len(on_land) > 10:
            print(f"  ... and {len(on_land)-10} more")

    print(f"\nSeeds outside bbox [{LAT_MIN},{LAT_MAX}]x[{LON_MIN},{LON_MAX}]: {len(out_bbox)}")
    if out_bbox:
        for idx, lat, lon in out_bbox[:10]:
            print(f"  seed {idx:4d}: lat={lat:.2f}, lon={lon:.2f}")

    # Plot
    fig, ax = plt.subplots(figsize=(14, 7))

    # Background: ocean mask (clipped to [0,1])
    # ocean_mask column 0 = lon≈0° (prime meridian).
    # Roll by LON_RES//2 so column 0 of the display = lon 180° = -180°,
    # which correctly aligns with extent=[-180, 180].
    display = np.clip(np.roll(ocean_mask, LON_RES // 2, axis=1), 0, 1)
    ax.imshow(display, cmap='ocean', origin='lower', vmin=0, vmax=1,
              extent=[-180, 180, -90, 90], aspect='auto')

    # Plot seeds: green = ok, red = on land
    land_idx_set = {i for i, _, _ in on_land}
    ok_lons  = [lons[i] for i in range(len(seeds)) if i not in land_idx_set]
    ok_lats  = [lats[i] for i in range(len(seeds)) if i not in land_idx_set]
    bad_lons = [lons[i] for i in land_idx_set]
    bad_lats = [lats[i] for i in land_idx_set]

    ax.scatter(ok_lons,  ok_lats,  s=4,  c='lime',   alpha=0.7, label='Ocean seed')
    if bad_lons:
        ax.scatter(bad_lons, bad_lats, s=20, c='red', alpha=1.0, label='Land seed (mask=0)')

    # Draw Atlantic bbox
    rect = mpatches.Rectangle((LON_MIN, LAT_MIN),
                                LON_MAX - LON_MIN, LAT_MAX - LAT_MIN,
                                linewidth=1.5, edgecolor='yellow',
                                facecolor='none', label='Atlantic bbox')
    ax.add_patch(rect)

    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)
    ax.set_xlabel('Longitude (degrees)')
    ax.set_ylabel('Latitude (degrees)')
    ax.set_title(f'Atlantic seeds check  (n={len(seeds)}, on_land={len(on_land)})')
    ax.legend(loc='lower left', fontsize=8)
    ax.grid(True, linewidth=0.3, alpha=0.5)

    outpng = 'atlantic_seeds_check.png'
    plt.tight_layout()
    plt.savefig(outpng, dpi=150)
    print(f"\nSaved: {outpng}")

if __name__ == '__main__':
    main()
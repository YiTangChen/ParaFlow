#!/usr/bin/env python3
"""
Read MPAS-O trace binary files, convert XYZ (meters) to lat/lon,
and plot all traces on a world map using matplotlib.

Usage:
    python3 plot_traces_latlon.py [trace_dir] [output.png]

Defaults:
    trace_dir  : ./trace_outputs
    output.png : traces_latlon.png
"""

import sys
import os
import struct
import glob
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm


def read_trace_file(filepath):
    """Read a single binary trace file, return list of Nx3 numpy arrays."""
    traces = []
    total_size = os.path.getsize(filepath)
    with open(filepath, 'rb') as f:
        data = f.read()
    offset = 0
    while offset < total_size:
        if offset + 4 > total_size:
            break
        nPts = struct.unpack_from('<i', data, offset)[0]
        offset += 4
        if nPts <= 0 or nPts > 10_000_000:
            print(f"  Warning: suspicious nPts={nPts} at offset {offset-4}, stopping")
            break
        byte_count = nPts * 3 * 8
        if offset + byte_count > total_size:
            print(f"  Warning: not enough data for {nPts} points, stopping")
            break
        coords = np.frombuffer(data, dtype='<f8', count=nPts * 3, offset=offset)
        traces.append(coords.reshape(nPts, 3))
        offset += byte_count
    return traces


def read_all_traces(trace_dir):
    """Read all traces_*.bin files, return list of Nx3 arrays."""
    pattern = os.path.join(trace_dir, 'traces_*.bin')
    files = sorted(glob.glob(pattern))
    if not files:
        raise FileNotFoundError(f"No trace files found: {pattern}")
    all_traces = []
    for fp in files:
        print(f"  Reading {os.path.basename(fp)} ...", end=' ', flush=True)
        traces = read_trace_file(fp)
        print(f"{len(traces)} traces")
        all_traces.extend(traces)
    return all_traces


def xyz_to_latlon(coords):
    """
    Convert Cartesian (x, y, z) in meters to (lat, lon) in degrees.
    lat = arcsin(z / r),  lon = arctan2(y, x)
    """
    x, y, z = coords[:, 0], coords[:, 1], coords[:, 2]
    r = np.sqrt(x**2 + y**2 + z**2)
    lat = np.degrees(np.arcsin(np.clip(z / r, -1.0, 1.0)))
    lon = np.degrees(np.arctan2(y, x))
    return lat, lon


def plot_traces(all_traces, output_path):
    """Plot all traces on a lat/lon world map."""
    n_traces = len(all_traces)
    colors = cm.plasma(np.linspace(0, 1, n_traces))

    fig, ax = plt.subplots(figsize=(18, 9))

    # World map background lines
    ax.set_facecolor('#0a0a1a')
    fig.patch.set_facecolor('#0a0a1a')

    # Grid lines
    for lat_line in range(-90, 91, 30):
        ax.axhline(lat_line, color='#333355', linewidth=0.4, zorder=1)
    for lon_line in range(-180, 181, 30):
        ax.axvline(lon_line, color='#333355', linewidth=0.4, zorder=1)

    print(f"Plotting {n_traces} traces...")
    for i, trace in enumerate(all_traces):
        lat, lon = xyz_to_latlon(trace)

        # Detect longitude wrap-around and split into segments
        dlon = np.abs(np.diff(lon))
        breaks = np.where(dlon > 180)[0] + 1
        segments = np.split(np.arange(len(lon)), breaks)

        for seg in segments:
            if len(seg) < 2:
                continue
            ax.plot(lon[seg], lat[seg],
                    color=colors[i], linewidth=0.4, alpha=0.7, zorder=2)

    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)
    ax.set_xlabel('Longitude (°)', color='white', fontsize=12)
    ax.set_ylabel('Latitude (°)', color='white', fontsize=12)
    ax.set_title(f'MPAS-O Particle Traces ({n_traces} traces)', color='white', fontsize=14)
    ax.tick_params(colors='white')
    for spine in ax.spines.values():
        spine.set_edgecolor('#555577')

    # Colorbar to distinguish traces by index
    sm = plt.cm.ScalarMappable(cmap='plasma', norm=plt.Normalize(0, n_traces - 1))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, orientation='vertical', fraction=0.02, pad=0.01)
    cbar.set_label('Trace Index', color='white', fontsize=10)
    cbar.ax.yaxis.set_tick_params(color='white')
    plt.setp(cbar.ax.yaxis.get_ticklabels(), color='white')

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight',
                facecolor=fig.get_facecolor())
    print(f"Saved: {output_path}")
    plt.close()


def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    trace_dir = sys.argv[1] if len(sys.argv) > 1 else \
        os.path.join(script_dir, 'trace_outputs')
    output = sys.argv[2] if len(sys.argv) > 2 else \
        os.path.join(script_dir, 'traces_latlon.png')

    print(f"Reading traces from: {trace_dir}")
    all_traces = read_all_traces(trace_dir)
    total_pts = sum(len(t) for t in all_traces)
    print(f"Total: {len(all_traces)} traces, {total_pts:,} points")

    plot_traces(all_traces, output)
    print("Done.")


if __name__ == '__main__':
    main()

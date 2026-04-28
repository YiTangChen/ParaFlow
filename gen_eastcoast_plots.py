"""
Generate east coast seed distribution plots (surface and depth separately).
Output: png/presentation/eastcoast_surface_seeds.png
        png/presentation/eastcoast_depth_seeds.png
        png/presentation/eastcoast_compare.png
"""

import os
import struct
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import netCDF4 as nc
import cartopy.crs as ccrs
import cartopy.feature as cfeature

MESH_FILE       = "/global/cfs/cdirs/m4259/yqiu/MOPS_ab/climatology/ocean.EC30to60E2r2.210210.nc"
SURFACE_SEEDS   = "seeds/eastcoast_surface_seeds.bin"
DEPTH_SEEDS     = "seeds/eastcoast_depth_seeds.bin"
OUT_DIR         = "png/presentation"
os.makedirs(OUT_DIR, exist_ok=True)

EARTH_RADIUS = 6371220.0

# East Coast map extent [lon_min, lon_max, lat_min, lat_max]
EXTENT = [-85, -55, 22, 48]

# ── helpers ───────────────────────────────────────────────────────────────────
def load_mesh_latlon():
    ds  = nc.Dataset(MESH_FILE, 'r')
    lat = np.degrees(ds.variables['latCell'][:])
    lon = np.degrees(ds.variables['lonCell'][:])
    ds.close()
    lon = np.where(lon > 180, lon - 360, lon)
    # keep only cells inside the East Coast view (+ small margin)
    margin = 3
    mask = ((lat > EXTENT[2] - margin) & (lat < EXTENT[3] + margin) &
            (lon > EXTENT[0] - margin) & (lon < EXTENT[1] + margin))
    return lat[mask], lon[mask]

def load_seeds(path):
    with open(path, 'rb') as f:
        data = f.read()
    n  = len(data) // 24
    v  = struct.unpack(f'{3*n}d', data)
    xs, ys, zs = np.array(v[0::3]), np.array(v[1::3]), np.array(v[2::3])
    r     = np.sqrt(xs**2 + ys**2 + zs**2)
    lat   = np.degrees(np.arcsin(zs / r))
    lon   = np.degrees(np.arctan2(ys, xs))
    depth = EARTH_RADIUS - r          # meters below surface
    return lat, lon, depth

def make_base_map(ax, title):
    ax.set_extent(EXTENT, crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.LAND,       facecolor='#d4d4d4', zorder=1)
    ax.add_feature(cfeature.OCEAN,      facecolor='#cce5ff', zorder=0)
    ax.add_feature(cfeature.COASTLINE,  linewidth=0.8,        zorder=3)
    ax.add_feature(cfeature.STATES,     linewidth=0.4, edgecolor='#888888', zorder=2)
    ax.add_feature(cfeature.BORDERS,    linewidth=0.6, edgecolor='#666666', zorder=2)
    gl = ax.gridlines(draw_labels=True, linewidth=0.4, color='gray',
                      alpha=0.6, linestyle='--', zorder=2)
    gl.top_labels   = False
    gl.right_labels = False
    gl.xlabel_style = {'size': 9}
    gl.ylabel_style = {'size': 9}
    ax.set_title(title, fontsize=13, fontweight='bold', pad=10)

def save(fig, name):
    path = os.path.join(OUT_DIR, name)
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  saved: {path}")


# ── load data ─────────────────────────────────────────────────────────────────
print("Loading mesh (East Coast region)...")
mesh_lat, mesh_lon = load_mesh_latlon()
print(f"  {len(mesh_lat):,} cells in view")

print("Loading surface seeds...")
surf_lat, surf_lon, surf_depth = load_seeds(SURFACE_SEEDS)
print(f"  depth range: {surf_depth.min():.1f} – {surf_depth.max():.1f} m")

print("Loading depth seeds...")
deep_lat, deep_lon, deep_depth = load_seeds(DEPTH_SEEDS)
print(f"  depth range: {deep_depth.min():.1f} – {deep_depth.max():.1f} m")


# ── Plot 1: Surface seeds ─────────────────────────────────────────────────────
print("Plotting surface seeds...")
fig = plt.figure(figsize=(10, 8))
ax  = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
make_base_map(ax, "East Coast Surface Seeds\n(Depth: 0–10 m)")

ax.scatter(mesh_lon, mesh_lat, s=1.5, color='#888888', alpha=0.25,
           transform=ccrs.PlateCarree(), linewidths=0, rasterized=True,
           zorder=2, label='Ocean mesh cells')

ax.scatter(surf_lon, surf_lat, s=6, color='#e53935', alpha=0.55,
           transform=ccrs.PlateCarree(), linewidths=0, rasterized=True,
           zorder=4, label=f'Surface seeds ({len(surf_lat):,})')

ax.legend(loc='lower right', fontsize=10, framealpha=0.9, markerscale=3)
ax.text(0.02, 0.04,
        f"Lat: 25–45°N  |  Lon: 80–60°W\nDepth: {surf_depth.min():.0f}–{surf_depth.max():.0f} m",
        transform=ax.transAxes, fontsize=9, color='#333333',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.85))

save(fig, "eastcoast_surface_seeds.png")


# ── Plot 2: Depth seeds colored by depth value ────────────────────────────────
print("Plotting depth seeds...")
fig = plt.figure(figsize=(10, 8))
ax  = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
make_base_map(ax, "East Coast Deep Seeds\n(Depth: 200–2000 m)")

ax.scatter(mesh_lon, mesh_lat, s=1.5, color='#888888', alpha=0.25,
           transform=ccrs.PlateCarree(), linewidths=0, rasterized=True,
           zorder=2)

sc = ax.scatter(deep_lon, deep_lat, c=deep_depth, s=6, alpha=0.65,
                transform=ccrs.PlateCarree(), linewidths=0, rasterized=True,
                zorder=4, cmap='Blues_r', vmin=deep_depth.min(), vmax=deep_depth.max())

cbar = plt.colorbar(sc, ax=ax, orientation='vertical', pad=0.02, shrink=0.7)
cbar.set_label('Depth (m)', fontsize=10)
cbar.ax.invert_yaxis()   # shallow at top, deep at bottom

ax.text(0.02, 0.04,
        f"Lat: 25–45°N  |  Lon: 80–60°W\nDepth: {deep_depth.min():.0f}–{deep_depth.max():.0f} m\n{len(deep_lat):,} seeds",
        transform=ax.transAxes, fontsize=9, color='#333333',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.85))

save(fig, "eastcoast_depth_seeds.png")


# ── Plot 3: Side-by-side comparison ──────────────────────────────────────────
print("Plotting east coast comparison...")
fig, axes = plt.subplots(1, 2, figsize=(18, 8),
                         subplot_kw={'projection': ccrs.PlateCarree()})

# Surface
ax = axes[0]
make_base_map(ax, "Surface Seeds\n(0–10 m depth)")
ax.scatter(mesh_lon, mesh_lat, s=1.2, color='#888888', alpha=0.2,
           transform=ccrs.PlateCarree(), linewidths=0, rasterized=True, zorder=2)
ax.scatter(surf_lon, surf_lat, s=5, color='#e53935', alpha=0.55,
           transform=ccrs.PlateCarree(), linewidths=0, rasterized=True,
           zorder=4, label=f'{len(surf_lat):,} seeds')
ax.legend(loc='lower right', fontsize=10, framealpha=0.9, markerscale=3)
ax.text(0.02, 0.04, "Gulf Stream region\nSurface layer particles",
        transform=ax.transAxes, fontsize=9,
        bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.85))

# Depth
ax = axes[1]
make_base_map(ax, "Deep Seeds\n(200–2000 m depth)")
ax.scatter(mesh_lon, mesh_lat, s=1.2, color='#888888', alpha=0.2,
           transform=ccrs.PlateCarree(), linewidths=0, rasterized=True, zorder=2)
sc = ax.scatter(deep_lon, deep_lat, c=deep_depth, s=5, alpha=0.65,
                transform=ccrs.PlateCarree(), linewidths=0, rasterized=True,
                zorder=4, cmap='Blues_r', vmin=deep_depth.min(), vmax=deep_depth.max())
cbar = plt.colorbar(sc, ax=ax, orientation='vertical', pad=0.02, shrink=0.7)
cbar.set_label('Depth (m)', fontsize=10)
cbar.ax.invert_yaxis()
ax.text(0.02, 0.04, f"Gulf Stream region\nDeep ocean particles",
        transform=ax.transAxes, fontsize=9,
        bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.85))

fig.suptitle("East Coast Focused Seeding Strategy — US East Coast / Gulf Stream Region",
             fontsize=14, fontweight='bold', y=1.01)
save(fig, "eastcoast_compare.png")

print(f"\nAll east coast plots saved to {OUT_DIR}/")

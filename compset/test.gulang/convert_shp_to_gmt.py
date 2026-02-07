#!/usr/bin/env python3
"""
Convert fault trace shapefiles (lon/lat) to .gmt.txt files (km).

Uses a local Cartesian projection centered at a reference point.
At latitude phi:
  1 deg lon ≈ cos(phi) * 111.32 km
  1 deg lat ≈ 111.32 km

Usage:
  python3 convert_shp_to_gmt.py

Input:  shapefiles (1.shp, 2.shp, ...) in shp_dir
Output: ft1.gmt.txt, ft2.gmt.txt, ... in out_dir
"""

import numpy as np
import shapefile
import os

# --- User settings ---
shp_dir = 'central gulang-zhongwei_simplify'
out_dir = 'user_fault_geometry_input'
fault_ids = [1, 2, 3, 4, 5]
fault_names = ['ft1', 'ft2', 'ft3', 'ft4', 'ft5']

# Reference point (lon0, lat0) for local km projection.
# Set to None to auto-compute from centroid of all faults.
ref_lon = None
ref_lat = None

# --- End settings ---

# Collect all points to compute reference if needed
all_lons = []
all_lats = []
fault_data = {}

for fid, fname in zip(fault_ids, fault_names):
    sf = shapefile.Reader(os.path.join(shp_dir, str(fid)))
    pts = sf.shapes()[0].points
    lons = [p[0] for p in pts]
    lats = [p[1] for p in pts]
    fault_data[fname] = (lons, lats)
    all_lons.extend(lons)
    all_lats.extend(lats)

# Auto-compute reference point
if ref_lon is None:
    ref_lon = np.mean(all_lons)
if ref_lat is None:
    ref_lat = np.mean(all_lats)

print(f"Reference point: lon={ref_lon:.6f}, lat={ref_lat:.6f}")

# Conversion factors
km_per_deg_lat = 111.32
km_per_deg_lon = 111.32 * np.cos(np.radians(ref_lat))
print(f"Scale: 1 deg lon = {km_per_deg_lon:.4f} km, 1 deg lat = {km_per_deg_lat:.4f} km")

# Convert and write
os.makedirs(out_dir, exist_ok=True)

for fname in fault_names:
    lons, lats = fault_data[fname]
    x_km = [(lon - ref_lon) * km_per_deg_lon for lon in lons]
    y_km = [(lat - ref_lat) * km_per_deg_lat for lat in lats]

    outfile = os.path.join(out_dir, fname + '.gmt.txt')
    with open(outfile, 'w') as f:
        for x, y in zip(x_km, y_km):
            f.write(f"{x:.6f} {y:.6f}\n")

    print(f"{fname}: {len(lons)} pts, x=[{min(x_km):.1f}, {max(x_km):.1f}] km, y=[{min(y_km):.1f}, {max(y_km):.1f}] km -> {outfile}")

print("\nDone. Reference point saved for reproducing the projection.")
print(f"ref_lon = {ref_lon:.6f}")
print(f"ref_lat = {ref_lat:.6f}")

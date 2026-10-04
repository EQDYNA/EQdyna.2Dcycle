#!/usr/bin/env python3
"""Convert SAF Model A model-frame coordinates to CFM UTM and lon/lat.

Liu et al. (2022) built the model geometry from SCEC CFM 5.2 traces, which are
in UTM zone 11, with the published MATLAB chain
(archive/published/.../fault_geometry/Figure1_Prep_SAF_SJF_Trace.m):

    model_km = rotate( (utm - [x0, y0]) / 1e3 , theta )
    x0 = 3.7e5 m, y0 = 3.8e6 m, theta = 40 deg
    rotate:  b1 = a1*cos(t) - a2*sin(t);  b2 = a1*sin(t) + a2*cos(t)

This module inverts that exactly, so a model point maps back into the same UTM
coordinate system the CFM traces use -- the system a ground-motion modeller
needs to sample SCEC CVM and CFM.

Lon/lat is computed by inverting the published deg2utm.m (WGS84 ellipsoid,
ported line for line below) numerically, so model -> lon/lat -> model round
trips to the published forward chain. CFM's own datum is not stated in the
published scripts; if the CFM coordinates are NAD27, the lon/lat here can be
off by ~100-200 m. The UTM coordinates carry no such ambiguity: they are CFM's.

Usage:
    saf_model_to_geo.py <case_dir|stf.bin*> [-o fault_nodes_geo.csv]
    saf_model_to_geo.py --check        # validate against the published numbers
"""

from __future__ import annotations

import argparse
import csv
import os
import sys

import numpy as np

X0, Y0, THETA_DEG = 3.7e5, 3.8e6, 40.0
ZONE = 11


def model_to_utm(x_m, y_m):
    """Model frame (m) -> UTM zone 11 easting, northing (m)."""
    b1 = np.asarray(x_m, float) / 1e3
    b2 = np.asarray(y_m, float) / 1e3
    t = np.radians(THETA_DEG)
    a1 = b1 * np.cos(t) + b2 * np.sin(t)          # inverse of rotate.m
    a2 = -b1 * np.sin(t) + b2 * np.cos(t)
    return a1 * 1e3 + X0, a2 * 1e3 + Y0           # inverse of convert.m


def utm_to_model(e, n):
    """UTM zone 11 (m) -> model frame (m), the published forward chain."""
    a1 = (np.asarray(e, float) - X0) / 1e3
    a2 = (np.asarray(n, float) - Y0) / 1e3
    t = np.radians(THETA_DEG)
    return (a1 * np.cos(t) - a2 * np.sin(t)) * 1e3, (a1 * np.sin(t) + a2 * np.cos(t)) * 1e3


def deg2utm(lat, lon, zone=None):
    """Line-for-line port of the published deg2utm.m (Coticchia-Surace, WGS84).

    deg2utm.m picks the zone from the longitude, so points west of -120 deg
    drop into zone 10. CFM keeps the whole system in zone 11 (the NW end of the
    SAF sits at easting ~200 km), so the inverse passes zone=11 to hold the
    central meridian fixed. zone=None reproduces deg2utm.m exactly.
    """
    lat_d = np.asarray(lat, float)
    lon_d = np.asarray(lon, float)
    sa, sb = 6378137.0, 6356752.314245
    e2 = ((sa ** 2 - sb ** 2) ** 0.5) / sb
    e2c = e2 ** 2
    c = sa ** 2 / sb
    la = np.radians(lat_d)
    lo = np.radians(lon_d)
    huso = np.trunc(lon_d / 6 + 31) if zone is None else np.full_like(lon_d, float(zone))
    s = huso * 6 - 183
    ds = lo - np.radians(s)
    a = np.cos(la) * np.sin(ds)
    eps = 0.5 * np.log((1 + a) / (1 - a))
    nu = np.arctan(np.tan(la) / np.cos(ds)) - la
    v = (c / (1 + e2c * np.cos(la) ** 2) ** 0.5) * 0.9996
    ta = (e2c / 2) * eps ** 2 * np.cos(la) ** 2
    a1 = np.sin(2 * la)
    a2 = a1 * np.cos(la) ** 2
    j2 = la + a1 / 2
    j4 = (3 * j2 + a2) / 4
    j6 = (5 * j4 + a2 * np.cos(la) ** 2) / 3
    alfa = 0.75 * e2c
    beta = (5 / 3) * alfa ** 2
    gama = (35 / 27) * alfa ** 3
    bm = 0.9996 * c * (la - alfa * j2 + beta * j4 - gama * j6)
    x = eps * v * (1 + ta / 3) + 500000
    y = nu * v * (1 + ta) + bm
    return x, np.where(y < 0, 9999999 + y, y), huso


def utm_to_lonlat(e, n, iters=30):
    """Invert deg2utm numerically (Newton), zone 11. Returns lon, lat in degrees."""
    e = np.asarray(e, float)
    n = np.asarray(n, float)
    lat = 34.0 + (n - 3.76e6) / 111e3             # rough start, southern California
    lon = -117.0 + (e - 5e5) / (111e3 * np.cos(np.radians(lat)))
    h = 1e-7
    for _ in range(iters):
        x, y, _ = deg2utm(lat, lon, ZONE)
        fx, fy = x - e, y - n
        if max(np.abs(fx).max(), np.abs(fy).max()) < 1e-6:
            break
        xa, ya, _ = deg2utm(lat + h, lon, ZONE)
        xo, yo, _ = deg2utm(lat, lon + h, ZONE)
        j11, j12 = (xa - x) / h, (xo - x) / h     # d(x)/d(lat), d(x)/d(lon)
        j21, j22 = (ya - y) / h, (yo - y) / h
        det = j11 * j22 - j12 * j21
        lat = lat - (j22 * fx - j12 * fy) / det
        lon = lon - (-j21 * fx + j11 * fy) / det
    return lon, lat


def check() -> int:
    """Validate against numbers produced by the published MATLAB chain."""
    ok = True
    # observed_sliprates.m lon/lat -> model x; MATLAB R2024b values recorded in
    # saf_result_utils.OBSERVED_EQDYNA_X_KM
    sites = {"Carrizo": (-119.8145, 35.2591, -164.315662175388),
             "Mojave S": (-117.957, 34.487, 21.621033567916),
             "San Bernardino S": (-116.896, 34.025, 129.481611408865)}
    for name, (lon, lat, x_ref) in sites.items():
        e, n, zone = deg2utm(lat, lon)
        xm, ym = utm_to_model(e, n)
        d = abs(xm / 1e3 - x_ref)
        lon2, lat2 = utm_to_lonlat(*model_to_utm(xm, ym))
        back = max(abs(lon2 - lon), abs(lat2 - lat)) * 111e3
        print(f"  {name:17s} model x {xm/1e3:12.6f} km (MATLAB {x_ref:12.6f}, diff {d:.1e} km)"
              f" | round trip {back:.1e} m | zone {int(zone)}")
        ok &= d < 1e-6 and back < 1e-3 and int(zone) == ZONE
    # NW end of the SAF, west of -120 deg: must stay in zone 11 (CFM convention)
    e0, n0 = 201787.687, 3960338.182
    lon0, lat0 = utm_to_lonlat(e0, n0)
    e1, n1, _ = deg2utm(lat0, lon0, ZONE)
    err = max(abs(e1 - e0), abs(n1 - n0))
    print(f"  NW end (zone-11 easting 202 km): lon {float(lon0):.5f}, lat {float(lat0):.5f}, round trip {err:.1e} m")
    ok &= -121.0 < float(lon0) < -120.0 and err < 1e-3
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("path", nargs="?", help="case directory or stf.bin* file")
    ap.add_argument("-o", "--out", default="fault_nodes_geo.csv")
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    if a.check:
        sys.exit(check())
    if not a.path:
        ap.error("give a case directory or stf.bin file, or --check")
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from stf_read import STFCase
    nodes = STFCase(a.path).nodes
    e, n = model_to_utm(nodes["x"], nodes["y"])
    lon, lat = utm_to_lonlat(e, n)
    # strike from the model tangent, rotated back into UTM (degrees clockwise from north)
    te, tn = model_to_utm(nodes["x"] + nodes["tx"], nodes["y"] + nodes["ty"])
    strike = (np.degrees(np.arctan2(te - e, tn - n)) + 360.0) % 360.0
    with open(a.out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["node", "fault", "local_index", "x_model_m", "y_model_m",
                    "utm11_easting_m", "utm11_northing_m", "lon_deg", "lat_deg",
                    "tangent_azimuth_deg", "length_m"])
        for i in range(len(nodes)):
            w.writerow([i + 1, nodes["fault"][i], nodes["local"][i],
                        f"{nodes['x'][i]:.3f}", f"{nodes['y'][i]:.3f}",
                        f"{e[i]:.3f}", f"{n[i]:.3f}", f"{lon[i]:.7f}", f"{lat[i]:.7f}",
                        f"{strike[i]:.3f}", f"{nodes['length'][i]:.3f}"])
    print(f"wrote {a.out}: {len(nodes)} fault nodes")


if __name__ == "__main__":
    main()

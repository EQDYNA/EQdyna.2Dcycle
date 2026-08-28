#!/usr/bin/env python3
"""Export the Xianshuihe faults as meshgen input files.

Writes `user_fault_geometry_input/<name>.gmt.txt`, two columns (x y, km) in
the local rotated frame, one file per fault — the format `meshgen.py`
consumes via `loadFtLoc`.

The first mesh uses the THROUGH-GOING CHAIN only: the five faults that
succeed one another along strike. The two splay strands (ft3 = seg3 and
ft4 = seg5 in the 7-fault decomposition) are deferred — they run parallel
to ft2 rather than continuing it, so they need a hand-built multi-surface
topology like `subei.gmsh.lite`'s, not the chain topology `gulang.gmsh.lite`
uses. Mesh names are renumbered mf1..mf5 to keep them distinct from the
7-fault ids used in the README and figures.

    mesh fault   7-fault id   polylines      length
    mf1          ft1          seg1 + seg2     74 km
    mf2          ft2          seg0 + seg4    118 km
    mf3          ft5          seg6            54 km
    mf4          ft6          seg7            75 km
    mf5          ft7          seg8            63 km

Usage:
    python3 export_fault_geometry.py [--all]   # --all also writes the strands
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline

from plot_fault_trace import (assign_fault_ids, chain, polyline_length_km,
                              principal_frame, read_kml, rotate, to_local_km)

# 7-fault id -> mesh name, for the through-going chain only.
# All seven faults, exported under their own ids. ft3 and ft4 are the splay
# strands: ft3 runs above the chain, ft4 below it, so each divides one side
# into a lens plus the outer block (see meshgen.py).
CHAIN = {f"ft{i}": f"ft{i}" for i in range(1, 8)}
STRANDS: dict[str, str] = {}


def even_arclength(p: np.ndarray, n: int | None = None) -> np.ndarray:
    """Resample a polyline to EVEN arc-length spacing, same point count.

    The KML digitisation is wildly uneven -- ft1 has a 49 m segment next to a
    9.9 km one, a 200:1 ratio. gmsh subdivides each control-point segment
    independently and gives every segment at least one element, so a short
    segment produces a short element no matter what dx says. That is what puts
    a 49 m element beside 369 m neighbours.

    This evens the spacing WITHOUT densifying. Densifying is a different thing
    and is known to break the mesh: gmsh's tolerance merges or drops points
    that are too close together, leaving orphaned fault nodes with no mass
    (see the note in saf.gmsh.lite/meshgen.py). Point count is preserved here,
    so the density is unchanged and only the distribution improves.

    x stays strictly increasing, since even resampling of a monotone curve is
    monotone, which meshgen's CubicSpline requires.
    """
    seg = np.hypot(*np.diff(p, axis=0).T)
    s_cum = np.concatenate([[0.0], np.cumsum(seg)])
    n = len(p) if n is None else n
    s_new = np.linspace(0.0, s_cum[-1], n)
    return np.column_stack([CubicSpline(s_cum, p[:, 0], bc_type="natural")(s_new),
                            CubicSpline(s_cum, p[:, 1], bc_type="natural")(s_new)])


def main() -> None:
    here = Path(__file__).resolve().parent
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--all", action="store_true",
                    help="also export the two splay strands")
    args = ap.parse_args()

    segs = read_kml(here / "user_fault_geometry_input" / "xianshuihe_fault_trace.kml")
    lon0 = np.concatenate([s[:, 0] for s in segs]).mean()
    lat0 = np.concatenate([s[:, 1] for s in segs]).mean()
    local = [to_local_km(s, lon0, lat0) for s in segs]
    centre, theta = principal_frame(np.vstack(local))
    rotated = [rotate(p, centre, theta) for p in local]

    wanted = dict(CHAIN)
    if args.all:
        wanted.update(STRANDS)

    out_dir = here / "user_fault_geometry_input"
    print(f"origin lon0={lon0:.4f} lat0={lat0:.4f}  rotation {np.degrees(theta):.2f} deg")
    print(f"{'mesh':>12} {'fault':>6} {'npts':>5} {'len_km':>8} {'x_lo':>8} {'x_hi':>8}")
    for fid, group in assign_fault_ids(rotated):
        if fid not in wanted:
            continue
        pieces, _ = chain(rotated, group)
        p = np.vstack(pieces)
        seg0 = np.hypot(*np.diff(p, axis=0).T)
        q = even_arclength(p)
        if np.any(np.diff(q[:, 0]) <= 0):
            print(f"  {fid}: even-spacing resample broke x-monotonicity, keeping raw")
        else:
            seg1 = np.hypot(*np.diff(q, axis=0).T)
            print(f"  {fid}: spacing ratio {seg0.max()/seg0.min():7.1f} -> "
                  f"{seg1.max()/seg1.min():5.2f}   min {seg0.min():.3f} -> "
                  f"{seg1.min():.3f} km")
            p = q
        if np.any(np.diff(p[:, 0]) <= 0):
            raise SystemExit(f"error: {fid} is not strictly increasing in x after "
                             f"merge dedup; meshgen's CubicSpline would reject it")
        path = out_dir / f"{wanted[fid]}.gmt.txt"
        np.savetxt(path, p, fmt="%.6f")
        print(f"{wanted[fid]:>12} {fid:>6} {len(p):>5} {polyline_length_km(p):>8.1f} "
              f"{p[:, 0].min():>8.1f} {p[:, 0].max():>8.1f}")
    print(f"\nWrote to {out_dir}")


if __name__ == "__main__":
    main()

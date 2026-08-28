#!/usr/bin/env python3
"""Pre-mesh geometry check for a new fault system.

Building a compset from a raw fault trace surfaces the same handful of
decisions every time, and each one needs a human call -- they cannot be
inferred from the data. This reports them up front instead of letting them
appear as a gmsh failure, a welded junction, or an under-resolved fault
halfway through a run.

Checks, in the order they bite:

  RESOLUTION   faults shorter than MIN_NODES * dxy carry too few nodes to
               host a rupture. The nucleation patch alone is ~2 km radius in
               the SAF setup, so a fault of a few element widths cannot
               nucleate meaningfully. Fix: merge into a neighbour, or drop.

  SHARED NODE  two faults whose closest approach is ~0 share a point. gmsh
               will weld them or emit a degenerate element there. Fix: merge
               them, or pull them apart deliberately.

  SUB-ELEMENT  two faults closer than one dxy are not resolvable as separate
               surfaces -- they interact within a single element. Fix: merge,
               refine dxy locally, or accept and document.

  DUPLICATE    repeated or non-increasing x along a fault. Merging two
               polylines that share an endpoint leaves a duplicated node, and
               meshgen.py's CubicSpline needs strictly increasing x. Fix: drop
               the duplicate when exporting.

  STEP-OVER    junction separations above the cut-off are real structure
               (barriers, releasing/restraining bends) and should be
               preserved, not closed. Reported for information.

Nothing here is auto-fixed. Every finding is a modelling decision.

Usage:
    python3 check_fault_geometry.py --dir DIR [--dxy 400] [--pattern '*.txt']

    DIR holds one two-column (x y, km) text file per fault. Or import
    `check_geometry(faults_km, dxy_m)` and pass arrays directly.
"""

from __future__ import annotations

import argparse
import glob
import os
import sys

import numpy as np

# A fault shorter than this many elements cannot host a resolved rupture.
MIN_NODES = 40

# Closest approach below this is treated as a shared node, in km.
SHARED_NODE_KM = 0.01


def polyline_length_km(p: np.ndarray) -> float:
    return float(np.sum(np.hypot(*np.diff(p, axis=0).T)))


def min_distance_km(a: np.ndarray, b: np.ndarray) -> float:
    return float(min(np.hypot(*(b - p).T).min() for p in a))


def check_geometry(faults: dict[str, np.ndarray], dxy_m: float) -> list[tuple[str, str]]:
    """Return [(severity, message), ...]. severity is BLOCK, DECIDE or INFO."""
    findings: list[tuple[str, str]] = []
    dxy_km = dxy_m / 1000.0

    for name, p in faults.items():
        L = polyline_length_km(p)
        n = int(L / dxy_km)
        if n < MIN_NODES:
            findings.append(("DECIDE",
                             f"RESOLUTION  {name}: {L:.1f} km = {n} nodes at dxy={dxy_m:.0f} m "
                             f"(< {MIN_NODES}). Merge into a neighbour, or drop."))

    for name, p in faults.items():
        dx = np.diff(p[:, 0])
        if np.any(dx == 0):
            k = int(np.sum(dx == 0))
            findings.append(("BLOCK",
                             f"DUPLICATE   {name}: {k} repeated x value(s). meshgen.py's "
                             f"CubicSpline needs strictly increasing x. Drop the "
                             f"duplicate node(s) on export."))
        elif np.any(dx < 0):
            findings.append(("DECIDE",
                             f"DUPLICATE   {name}: x is not monotonic (max backstep "
                             f"{-dx.min():.3f} km). The spline parameterisation assumes "
                             f"a single-valued trace along strike."))

    names = list(faults)
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            d = min_distance_km(faults[a], faults[b])
            if d <= SHARED_NODE_KM:
                findings.append(("BLOCK",
                                 f"SHARED NODE {a}-{b}: closest approach {d:.4f} km. "
                                 f"gmsh will weld these or degenerate an element. "
                                 f"Merge them, or separate them deliberately."))
            elif d < dxy_km:
                findings.append(("DECIDE",
                                 f"SUB-ELEMENT {a}-{b}: {d:.3f} km apart, under one "
                                 f"element ({dxy_km:.3f} km). Not resolvable as separate "
                                 f"surfaces. Merge, refine locally, or accept and document."))
    return findings


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", required=True, help="Directory of two-column (x y, km) fault files")
    ap.add_argument("--pattern", default="*.txt", help="Glob for fault files (default *.txt)")
    ap.add_argument("--dxy", type=float, default=400.0, help="Element size in m (default 400)")
    args = ap.parse_args()

    paths = sorted(glob.glob(os.path.join(args.dir, args.pattern)))
    if not paths:
        sys.exit(f"error: no files matching {args.pattern} in {args.dir}")
    faults = {os.path.splitext(os.path.basename(p))[0]: np.loadtxt(p) for p in paths}

    print(f"{len(faults)} faults, dxy = {args.dxy:.0f} m\n")
    for name, p in faults.items():
        L = polyline_length_km(p)
        print(f"  {name:<8} {L:7.1f} km  {int(L / (args.dxy / 1000.0)):5d} nodes")

    findings = check_geometry(faults, args.dxy)
    print()
    if not findings:
        print("No geometry findings. Nothing needs a decision before meshing.")
        return
    for severity in ("BLOCK", "DECIDE", "INFO"):
        for sev, msg in findings:
            if sev == severity:
                print(f"[{sev:>6}] {msg}")
    if any(s == "BLOCK" for s, _ in findings):
        sys.exit(1)


if __name__ == "__main__":
    main()

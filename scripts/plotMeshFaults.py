#!/usr/bin/env python3
"""Mesh figures per fault segment and at every fault tip.

Two views that `plotMeshNearFault.py` does not give:

  --mode segments   one panel per fault, the whole fault end to end
  --mode tips       one panel per fault TIP, both ends of every fault

The tip view exists because `plotMeshNearFault.py` zooms on each fault's
MIDPOINT, so defects at the ends stay invisible -- and the ends are exactly
where an embedded fault's mesh misbehaves, since gmsh has to close the front
around a free tip. Every orphaned split node found on this project sat at a
tip (issue #1).

Cells failing the `checkMeshQuality.py` gates (interior angle < 20 or > 160
deg) are filled red in both views, so a bad cell is visible rather than
inferred from a table.

Usage:
    python3 plotMeshFaults.py <case_dir> [--mode segments|tips|both]
                              [--out PREFIX] [--zoom KM] [--dpi N]
"""

from __future__ import annotations

import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PolyCollection

# ---- frozen style ---------------------------------------------------------
# Keep these fixed so figures from different runs are directly comparable.
MM = 1.0 / 25.4
PAGE_W_MM = 190.0          # journal double-column width
STYLE = {
    "font.size": 11, "axes.titlesize": 11.5, "axes.labelsize": 12,
    "xtick.labelsize": 9, "ytick.labelsize": 9, "legend.fontsize": 10,
    "savefig.bbox": "tight", "pdf.fonttype": 42, "ps.fonttype": 42,
}
CELL_EDGE = "0.6"
CELL_LW = 0.5
BAD_FACE = "tab:red"
BAD_EDGE = "darkred"
FAULT_LW = 2.4
TIP_MS = 7
MIN_ANGLE_DEG = 20.0       # matches checkMeshQuality
MAX_ANGLE_DEG = 160.0
TIP_ZOOM_KM = 1.5


def _find(root, name):
    for d in (root, os.path.join(root, "fem_mesh_output")):
        p = os.path.join(d, name)
        if os.path.exists(p):
            return p
    raise SystemExit(f"error: {name} not found in {root} or its fem_mesh_output/")


def load(root):
    vert = np.loadtxt(_find(root, "vert.txt"))[:, :2]
    fac = np.loadtxt(_find(root, "fac.txt"), dtype=int)
    if fac.min() == 1:
        fac = fac - 1
    nsmp = np.loadtxt(_find(root, "nsmp.txt"), dtype=int)
    if nsmp.min() == 1:
        nsmp = nsmp - 1
    rows = [l.split() for l in open(_find(root, "meshGeneralInfo.txt")) if l.strip()]
    nfn = [int(x) for x in rows[1]]
    mx = max(nfn)
    faults = {i + 1: vert[nsmp[i * mx: i * mx + n, 0]] for i, n in enumerate(nfn)}
    return vert, fac, faults


def cell_angles(vert, fac):
    p = vert[fac]
    out = np.empty((len(fac), 4))
    for i in range(4):
        a = p[:, (i - 1) % 4] - p[:, i]
        b = p[:, (i + 1) % 4] - p[:, i]
        cos = (a * b).sum(1) / (np.linalg.norm(a, axis=1) * np.linalg.norm(b, axis=1) + 1e-30)
        out[:, i] = np.degrees(np.arccos(np.clip(cos, -1, 1)))
    return out


def draw(ax, vert, polys, bad, faults, xlim, ylim, colours):
    cen = polys.mean(axis=1)
    m = ((cen[:, 0] > xlim[0]) & (cen[:, 0] < xlim[1]) &
         (cen[:, 1] > ylim[0]) & (cen[:, 1] < ylim[1]))
    ax.add_collection(PolyCollection(polys[m], facecolors="none",
                                     edgecolors=CELL_EDGE, lw=CELL_LW))
    mb = m & bad
    if mb.any():
        ax.add_collection(PolyCollection(polys[mb], facecolors=BAD_FACE, alpha=0.5,
                                         edgecolors=BAD_EDGE, lw=1.2))
    for k, q in faults.items():
        sel = ((q[:, 0] > xlim[0] - 1) & (q[:, 0] < xlim[1] + 1) &
               (q[:, 1] > ylim[0] - 1) & (q[:, 1] < ylim[1] + 1))
        if sel.any():
            ax.plot(q[sel, 0], q[sel, 1], "-", lw=FAULT_LW, color=colours[(k - 1) % 10])
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_aspect("equal")
    ax.grid(alpha=0.2, lw=0.5)
    return int(mb.sum())


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("case_dir")
    ap.add_argument("--mode", choices=["segments", "tips", "both"], default="both")
    ap.add_argument("--out", default=None, help="Output prefix (default <case>/aPlots/mesh)")
    ap.add_argument("--zoom", type=float, default=TIP_ZOOM_KM, help="Tip half-window, km")
    ap.add_argument("--dpi", type=int, default=200)
    args = ap.parse_args()

    vert, fac, faults = load(args.case_dir)
    polys = vert[fac]
    A = cell_angles(vert, fac)
    bad = (A.min(axis=1) < MIN_ANGLE_DEG) | (A.max(axis=1) > MAX_ANGLE_DEG)
    print(f"{len(vert)} nodes, {len(fac)} cells, {len(faults)} faults, "
          f"{int(bad.sum())} cell(s) failing the angle gates")

    plt.rcParams.update(STYLE)
    colours = plt.cm.tab10(np.linspace(0, 1, 10))
    prefix = args.out or os.path.join(args.case_dir, "aPlots", "mesh")
    os.makedirs(os.path.dirname(os.path.abspath(prefix)), exist_ok=True)

    if args.mode in ("segments", "both"):
        n = len(faults)
        fig, axes = plt.subplots(n, 1, figsize=(PAGE_W_MM * MM, 42 * MM * n))
        axes = np.atleast_1d(axes)
        for i, (k, q) in enumerate(sorted(faults.items())):
            pad = 1.0
            nb = draw(axes[i], vert, polys, bad, faults,
                      (q[:, 0].min() - pad, q[:, 0].max() + pad),
                      (q[:, 1].min() - pad, q[:, 1].max() + pad), colours)
            axes[i].set_title(f"ft{k}: {len(q)} fault nodes, "
                              f"{np.hypot(*np.diff(q, axis=0).T).sum():.1f} km"
                              + (f", {nb} BAD CELL(S)" if nb else ""))
            axes[i].set_ylabel("y (km)")
        axes[-1].set_xlabel("x (km)")
        fig.suptitle("Mesh per fault segment (red = cell failing the angle gates)",
                     fontsize=13, y=0.999)
        out = f"{prefix}_segments.png"
        fig.savefig(out, dpi=args.dpi)
        print(f"Wrote {out}")
        plt.close(fig)

    if args.mode in ("tips", "both"):
        n = len(faults)
        fig, axes = plt.subplots(n, 2, figsize=(PAGE_W_MM * MM, 47 * MM * n))
        axes = np.atleast_2d(axes)
        for i, (k, q) in enumerate(sorted(faults.items())):
            for j, (tag, tip) in enumerate((("start", q[0]), ("end", q[-1]))):
                ax = axes[i, j]
                nb = draw(ax, vert, polys, bad, faults,
                          (tip[0] - args.zoom, tip[0] + args.zoom),
                          (tip[1] - args.zoom, tip[1] + args.zoom), colours)
                ax.plot(*tip, "o", ms=TIP_MS, mfc="none", mec="k", mew=1.8)
                ax.set_title(f"ft{k} {tag} ({tip[0]:.1f},{tip[1]:.1f})"
                             + (f"  {nb} BAD" if nb else ""))
            axes[i, 0].set_ylabel("y (km)")
        for ax in axes[-1]:
            ax.set_xlabel("x (km)")
        fig.suptitle("Mesh at every fault tip (circle = tip; red = cell failing the gates)",
                     fontsize=13, y=0.999)
        out = f"{prefix}_tips.png"
        fig.savefig(out, dpi=args.dpi)
        print(f"Wrote {out}")
        plt.close(fig)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Render the whole FE mesh with its faults: domain overview plus a zoom.

`plotMeshNearFault.py` zooms on each fault; this is its whole-domain
counterpart. Drawing 17k cell outlines across a 500+ km domain produces
solid grey, so the overview shades each cell by its element size
(sqrt of area) on a log scale instead. That shows what actually matters
at domain scale -- how the mesh grades from the fine near-fault resolution
out to the coarse boundary -- and stays readable at print size.

Outlines are drawn only in the zoom panel, where cells are large enough on
the page to be distinguishable.

Usage:
    python3 plotMeshDomain.py <case_dir> [--out PATH] [--zoom X0 X1 Y0 Y1]
                              [--edges] [--dpi 250]
"""

from __future__ import annotations

import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PolyCollection

MM = 1.0 / 25.4


def _find(root: str, name: str) -> str:
    for d in (root, os.path.join(root, "fem_mesh_output")):
        p = os.path.join(d, name)
        if os.path.exists(p):
            return p
    sys.exit(f"error: {name} not found in {root} or its fem_mesh_output/")


def load_mesh(root: str):
    vert = np.loadtxt(_find(root, "vert.txt"))
    # fac.txt and nsmp.txt are written 0-based by meshgen.py but 1-based by
    # some older cases, so shift only when the minimum is 1 -- the same test
    # plotMeshNearFault.py and checkMeshQuality.py use. Subtracting
    # unconditionally scrambles every cell (it inflated the total mesh area
    # 150x when this script first did it).
    fac = np.loadtxt(_find(root, "fac.txt"), dtype=int)
    if fac.min() == 1:
        fac = fac - 1
    nsmp = np.loadtxt(_find(root, "nsmp.txt"), dtype=int)
    if nsmp.min() == 1:
        nsmp = nsmp - 1
    with open(_find(root, "meshGeneralInfo.txt")) as f:
        rows = [l.split() for l in f if l.strip()]
    nfn = [int(x) for x in rows[1]]
    rng = [float(x) for x in rows[2]] if len(rows) > 2 else None
    return vert[:, :2], fac[:, :4], nsmp, nfn, rng


def fault_polylines(vert, nsmp, nfn):
    """One (n,2) array per fault, from the master-node column.

    nsmp is padded to maxftnode per fault, so take exactly the first nfn[i]
    rows of each block. Do NOT filter on id > 0: with 0-based indexing node 0
    is a real node.
    """
    maxftnode = max(nfn)
    return [vert[nsmp[i * maxftnode: i * maxftnode + n, 0]]
            for i, n in enumerate(nfn)]


def quad_size_km(polys: np.ndarray) -> np.ndarray:
    """Representative element size: sqrt of the shoelace area, per cell."""
    x, y = polys[:, :, 0], polys[:, :, 1]
    area = 0.5 * np.abs(np.sum(x * np.roll(y, -1, axis=1) - np.roll(x, -1, axis=1) * y, axis=1))
    return np.sqrt(area)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("case_dir")
    ap.add_argument("--out", default=None, help="Output PNG (default <case>/aPlots/meshDomain.png)")
    ap.add_argument("--zoom", nargs=4, type=float, metavar=("X0", "X1", "Y0", "Y1"),
                    default=None, help="Zoom window in km (default: 40 km box on fault 1)")
    ap.add_argument("--edges", action="store_true", help="Also outline cells in the overview")
    ap.add_argument("--dpi", type=int, default=250)
    args = ap.parse_args()

    vert, fac, nsmp, nfn, rng = load_mesh(args.case_dir)
    polys = vert[fac]
    size = quad_size_km(polys)
    faults = fault_polylines(vert, nsmp, nfn)

    print(f"{len(vert)} nodes, {len(fac)} cells, {len(nfn)} faults, nfnodes={nfn}")
    print(f"element size km: min={size.min():.3f} med={np.median(size):.3f} max={size.max():.3f}")

    if args.zoom is None:
        c = faults[0].mean(axis=0)
        args.zoom = [c[0] - 20, c[0] + 20, c[1] - 20, c[1] + 20]
    x0, x1, y0, y1 = args.zoom

    plt.rcParams.update({"font.size": 13, "axes.titlesize": 15, "axes.labelsize": 14,
                         "xtick.labelsize": 12, "ytick.labelsize": 12,
                         "legend.fontsize": 12, "savefig.bbox": "tight"})
    fig = plt.figure(figsize=(190 * MM, 205 * MM))
    gs = fig.add_gridspec(2, 1, height_ratios=[1.0, 1.0], hspace=0.30)
    ax_a, ax_b = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])

    pc = PolyCollection(polys, array=size, cmap="viridis",
                        norm=matplotlib.colors.LogNorm(vmax=size.max()),
                        edgecolors=("0.5" if args.edges else "none"),
                        linewidths=0.08 if args.edges else 0.0)
    ax_a.add_collection(pc)
    cb = fig.colorbar(pc, ax=ax_a, pad=0.02, fraction=0.035)
    cb.set_label("element size (km)")
    cols = plt.cm.tab10(np.linspace(0, 1, 10))
    for i, f in enumerate(faults):
        ax_a.plot(f[:, 0], f[:, 1], "-", lw=2.0, color="w", zorder=4)
        ax_a.plot(f[:, 0], f[:, 1], "-", lw=1.2, color=cols[i % 10], zorder=5,
                  label=f"ft{i + 1} ({nfn[i]})")
    if rng:
        ax_a.set_xlim(rng[0], rng[1]); ax_a.set_ylim(rng[2], rng[3])
    else:
        ax_a.autoscale_view()
    ax_a.set_aspect("equal"); ax_a.set_xlabel("x (km)"); ax_a.set_ylabel("y (km)")
    ax_a.set_title("Whole domain, shaded by element size")
    ax_a.legend(loc="lower left", frameon=False, ncol=3, fontsize=10)
    ax_a.text(0.012, 0.97, "(a)", transform=ax_a.transAxes, fontsize=16,
              fontweight="bold", va="top")

    cen = polys.mean(axis=1)
    m = (cen[:, 0] > x0) & (cen[:, 0] < x1) & (cen[:, 1] > y0) & (cen[:, 1] < y1)
    ax_b.add_collection(PolyCollection(polys[m], facecolors="none",
                                       edgecolors="0.45", linewidths=0.4))
    for i, f in enumerate(faults):
        ax_b.plot(f[:, 0], f[:, 1], "-", lw=2.4, color=cols[i % 10], zorder=5)
    ax_b.set_xlim(x0, x1); ax_b.set_ylim(y0, y1)
    ax_b.set_aspect("equal"); ax_b.set_xlabel("x (km)"); ax_b.set_ylabel("y (km)")
    ax_b.set_title(f"Zoom: {x1 - x0:.0f} x {y1 - y0:.0f} km, {int(m.sum())} cells")
    ax_b.text(0.012, 0.97, "(b)", transform=ax_b.transAxes, fontsize=16,
              fontweight="bold", va="top")

    fig.suptitle(f"{os.path.basename(os.path.abspath(args.case_dir))}: "
                 f"{len(vert)} nodes, {len(fac)} quads", fontsize=15, y=0.985)
    out = args.out or os.path.join(args.case_dir, "aPlots", "meshWithFaults.png")
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    fig.savefig(out, dpi=args.dpi)
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Parse and plot the Xianshuihe fault trace from the 1:1M active-fault KML.

Reads `user_fault_geometry_input/xianshuihe_fault_trace.kml` (Chinese
1:1,000,000 active fault database, 鲜水河断裂), projects it into a local
Cartesian frame, and plots both the raw digitisation and the proposed
EQdyna fault decomposition.

The KML holds 9 separate polylines across 3 Placemarks. They are not one
continuous line: the trace splits into sub-parallel strands between about
30.0 and 30.5 N (the Yalahe / Selaha / Zheduotang strands near Kangding).
Endpoint connectivity, in km:

    seg1 -> seg2  0.26     seg2 -> seg0  1.19     seg0 -> seg4  0.00
    seg5 -> seg6  1.62     seg7 -> seg8  4.41

seg4 and seg5 both leave the trace near (101.95-102.0 E, 29.9-30.07 N) and
rejoin near (101.52-101.56 E, 30.43-30.54 N). Neither gives a gap-free
through-going path: seg4 joins seg0 exactly (0.00 km) but stops 12.6 km short
of seg6, while seg5 reaches seg6 (1.6 km) but starts 19.8 km off seg0. ft1 is
therefore routed via seg4 (smaller maximum gap) and seg5 is carried as the
parallel strand. The remaining gaps are reported, and drawn dotted with an
x, rather than silently bridged -- they are a real decision the mesh build
has to make, not a plotting artifact.

Usage:
    python3 plot_fault_trace.py [--kml PATH] [--output PATH]
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

R_EARTH_KM = 6371.0

# Proposed EQdyna fault decomposition: name -> (KML polyline indices, colour).
# ft1 is routed via seg4 rather than seg5: seg4 joins seg0 at 0.00 km and
# leaves a 12.6 km gap into seg6, against 19.8 km if seg5 is used instead.
# seg5 is then the parallel strand -- note it is the strand that actually
# reaches seg6 (1.6 km), so the two together bracket the Kangding splay zone
# and the gap has to be closed when the mesh is built.
FAULT_DECOMPOSITION = {
    "ft1 through-going (segs 1,2,0,4,6,7,8)": ([1, 2, 0, 4, 6, 7, 8], "k"),
    "ft2 parallel strand (seg 5)":            ([5], "tab:red"),
    "ft3 short splay (seg 3)":                ([3], "tab:blue"),
}


def read_kml(path: Path) -> list[np.ndarray]:
    """Return one (n, 2) lon/lat array per <coordinates> block, in file order."""
    txt = path.read_text(encoding="utf-8")
    segs = []
    for placemark in re.findall(r"<Placemark>(.*?)</Placemark>", txt, re.S):
        for block in re.findall(r"<coordinates>(.*?)</coordinates>", placemark, re.S):
            pts = [tuple(map(float, tok.split(",")[:2])) for tok in block.split()]
            segs.append(np.array(pts))
    if not segs:
        raise SystemExit(f"error: no <coordinates> found in {path}")
    return segs


def to_local_km(lonlat: np.ndarray, lon0: float, lat0: float) -> np.ndarray:
    """Equirectangular projection about (lon0, lat0), in km."""
    x = np.radians(lonlat[:, 0] - lon0) * R_EARTH_KM * np.cos(np.radians(lat0))
    y = np.radians(lonlat[:, 1] - lat0) * R_EARTH_KM
    return np.column_stack([x, y])


def principal_frame(points: np.ndarray):
    """Rotation that puts the trace's principal axis along +x. Returns (centre, theta)."""
    centre = points.mean(axis=0)
    axis = np.linalg.svd(points - centre, full_matrices=False)[2][0]
    return centre, float(np.arctan2(axis[1], axis[0]))


def rotate(points: np.ndarray, centre: np.ndarray, theta: float) -> np.ndarray:
    d = points - centre
    return np.column_stack([d[:, 0] * np.cos(-theta) - d[:, 1] * np.sin(-theta),
                            d[:, 0] * np.sin(-theta) + d[:, 1] * np.cos(-theta)])


GAP_TOLERANCE_KM = 2.0


def chain(rotated: list[np.ndarray], indices: list[int]):
    """Order polylines along strike and orient each along +x.

    Returns (pieces, gaps): `pieces` is the list of oriented polylines in
    along-strike order, `gaps` is [(end_of_piece, start_of_next, distance_km)]
    for every junction wider than GAP_TOLERANCE_KM.

    The pieces are deliberately NOT concatenated into one array. Joining them
    would draw a straight line across every gap and make a discontinuous
    routing look like a continuous fault. Do not sort the points by x either --
    the trace is not monotonic in x, and sorting scrambles the geometry.
    """
    ordered = sorted(indices, key=lambda i: rotated[i][:, 0].mean())
    pieces = []
    for i in ordered:
        seg = rotated[i]
        pieces.append(seg[::-1] if seg[0, 0] > seg[-1, 0] else seg)
    gaps = []
    for a, b in zip(pieces, pieces[1:]):
        d = float(np.hypot(*(b[0] - a[-1])))
        if d > GAP_TOLERANCE_KM:
            gaps.append((a[-1], b[0], d))
    return pieces, gaps


def polyline_length_km(points: np.ndarray) -> float:
    return float(np.sum(np.hypot(*np.diff(points, axis=0).T)))


def main() -> None:
    here = Path(__file__).resolve().parent
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--kml", type=Path,
                    default=here / "user_fault_geometry_input" / "xianshuihe_fault_trace.kml")
    ap.add_argument("--output", type=Path, default=here / "xianshuihe_fault_trace.png")
    args = ap.parse_args()

    segs = read_kml(args.kml)
    lon0 = np.concatenate([s[:, 0] for s in segs]).mean()
    lat0 = np.concatenate([s[:, 1] for s in segs]).mean()
    local = [to_local_km(s, lon0, lat0) for s in segs]
    centre, theta = principal_frame(np.vstack(local))
    rotated = [rotate(p, centre, theta) for p in local]

    print(f"kml      : {args.kml}")
    print(f"origin   : lon0={lon0:.4f} lat0={lat0:.4f}")
    print(f"rotation : {np.degrees(theta):.2f} deg (principal axis -> +x)")
    print(f"{'seg':>4} {'npts':>5} {'len_km':>8}   start(lon,lat) -> end(lon,lat)")
    for i, (s, p) in enumerate(zip(segs, local)):
        print(f"{i:>4} {len(s):>5} {polyline_length_km(p):>8.1f}   "
              f"({s[0,0]:.3f},{s[0,1]:.3f}) -> ({s[-1,0]:.3f},{s[-1,1]:.3f})")

    # Publication settings: figure sized to a 190 mm double-column text width,
    # with point sizes chosen to stay legible at that printed size.
    MM = 1.0 / 25.4
    plt.rcParams.update({
        "font.size": 13, "axes.titlesize": 15, "axes.labelsize": 14,
        "xtick.labelsize": 12.5, "ytick.labelsize": 12.5, "legend.fontsize": 12.5,
        "axes.linewidth": 1.0, "xtick.major.width": 1.0, "ytick.major.width": 1.0,
        "xtick.direction": "out", "ytick.direction": "out",
        "savefig.bbox": "tight", "pdf.fonttype": 42, "ps.fonttype": 42,
    })
    fig = plt.figure(figsize=(190 * MM, 245 * MM))
    gs = fig.add_gridspec(2, 1, height_ratios=[2.35, 1.0], hspace=0.30)
    ax_map, ax_loc = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])

    colours = plt.cm.tab10(np.linspace(0, 1, 10))
    for i, s in enumerate(segs):
        ax_map.plot(s[:, 0], s[:, 1], "-o", ms=4.5, lw=2.2, color=colours[i % 10])
        ax_map.annotate(f"seg{i}", (s[len(s) // 2, 0], s[len(s) // 2, 1]),
                        fontsize=13, fontweight="bold", color=colours[i % 10],
                        xytext=(7, 7), textcoords="offset points")
    ax_map.set_xlabel("Longitude ($^\\circ$E)")
    ax_map.set_ylabel("Latitude ($^\\circ$N)")
    ax_map.set_aspect(1.0 / np.cos(np.radians(lat0)))
    ax_map.grid(alpha=0.25, lw=0.7)
    ax_map.annotate("double strand\n(Yalahe / Selaha / Zheduotang)",
                    (101.35, 30.10), xytext=(100.72, 29.72), fontsize=12.5, ha="center",
                    bbox=dict(fc="w", ec="0.45", lw=1.0, alpha=0.95, pad=0.45),
                    arrowprops=dict(arrowstyle="->", color="0.25", lw=1.3))
    ax_map.text(0.015, 0.985, "(a)", transform=ax_map.transAxes, fontsize=16,
                fontweight="bold", va="top", ha="left")
    ax_map.set_title("Xianshuihe fault, 1:1M active-fault database", pad=10)

    print(f"\nfault decomposition (gaps wider than {GAP_TOLERANCE_KM} km are flagged):")
    for name, (indices, colour) in FAULT_DECOMPOSITION.items():
        pieces, gaps = chain(rotated, indices)
        total = sum(polyline_length_km(p) for p in pieces)
        for k, piece in enumerate(pieces):
            ax_loc.plot(piece[:, 0], piece[:, 1], "-o", ms=4.5, lw=2.2, color=colour,
                        label=f"{name}   {total:.0f} km" if k == 0 else None)
        # Stagger the labels: consecutive gaps sit close enough along strike
        # that fixed-offset text overlaps at publication font sizes.
        for k, (a, b, d) in enumerate(gaps):
            ax_loc.plot([a[0], b[0]], [a[1], b[1]], ":", lw=2.2, color=colour)
            ax_loc.plot([a[0], b[0]], [a[1], b[1]], "x", ms=10, mew=2.4, color=colour)
            dy = 16 if k % 2 == 0 else -26
            ax_loc.annotate(f"{d:.0f} km gap", ((a[0] + b[0]) / 2, (a[1] + b[1]) / 2),
                            fontsize=12.5, color=colour, fontweight="bold",
                            xytext=(0, dy), textcoords="offset points", ha="center")
        print(f"  {name:<42} {total:6.1f} km, {len(pieces)} piece(s), {len(gaps)} gap(s)"
              + ("".join(f"\n      gap {d:.1f} km at x={a[0]:.0f}" for a, b, d in gaps)))

    ax_loc.set_xlabel(f"Along-strike distance (km), frame rotated {np.degrees(theta):.1f}$^\\circ$")
    ax_loc.set_ylabel("Fault-normal (km)")
    ax_loc.grid(alpha=0.25, lw=0.7)
    ax_loc.legend(loc="upper center", bbox_to_anchor=(0.5, -0.30), frameon=False,
                  handlelength=2.4, borderaxespad=0.0)
    ax_loc.text(0.015, 0.97, "(b)", transform=ax_loc.transAxes, fontsize=16,
                fontweight="bold", va="top", ha="left")
    ax_loc.set_title("Proposed EQdyna fault decomposition", pad=10)

    # Fill panel (a) to the same width as panel (b) without distorting the map:
    # keep the true geographic aspect and widen the longitude range to match the
    # axes box, rather than letting the tall lat range squeeze the panel narrow.
    # Use the slot the gridspec ALLOTTED the axes (original=True), not the box
    # the equal-aspect constraint shrank it to -- the latter just reproduces the
    # current narrow shape and the panel never widens.
    fig.canvas.draw()
    slot = ax_map.get_position(original=True)
    fig_w, fig_h = fig.get_size_inches()
    box_ratio = (slot.width * fig_w) / (slot.height * fig_h)
    aspect = 1.0 / np.cos(np.radians(lat0))
    y0, y1 = ax_map.get_ylim()
    dx_needed = box_ratio * (y1 - y0) * aspect
    lon_c = 0.5 * sum(ax_map.get_xlim())
    ax_map.set_xlim(lon_c - dx_needed / 2, lon_c + dx_needed / 2)

    fig.savefig(args.output, dpi=300)
    fig.savefig(args.output.with_suffix(".pdf"))
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()

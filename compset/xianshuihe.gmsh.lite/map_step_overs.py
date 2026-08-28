#!/usr/bin/env python3
"""Map and classify step-overs along the Xianshuihe fault trace.

Qiao et al. (2022, EPSL 596, 117799) use step-over density as a proxy for
fault structural maturity, counting step-overs wider than 1% of the fault
length. They report 4/350 = 0.01 per km for the Xianshuihe against
1/550 = 0.002 for the Yushu-Ganzi fault, and argue the Xianshuihe is the
less mature of the two. Step-overs are therefore a feature to preserve in
the model, not a digitisation gap to close.

For each junction along the through-going trace this reports:

  gap     along-strike separation (negative where the segments overlap)
  width   fault-normal separation, measured against the LOCAL strike at the
          end of the incoming segment rather than the global rotation
  sense   left- or right-stepping, taken looking along the direction of travel
  kind    releasing or restraining. The Xianshuihe is left-lateral
          (sinistral), so a LEFT step opens (releasing, pull-apart) and a
          RIGHT step closes (restraining). This is the opposite of the more
          familiar dextral convention.

Usage:
    python3 map_step_overs.py [--kml PATH] [--output PATH]
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from plot_fault_trace import (FAULT_DECOMPOSITION, chain, polyline_length_km,
                              principal_frame, read_kml, rotate, to_local_km)

# Qiao et al. count step-overs wider than this fraction of the fault length.
STEP_WIDTH_FRACTION = 0.01

# Sinistral fault: left step opens, right step closes.
SENSE_KIND = {"left": "releasing", "right": "restraining"}


def local_unit(seg: np.ndarray, at: str, n: int = 4) -> np.ndarray:
    """Unit vector along the local trend at one end of a polyline."""
    pts = seg[-n:] if at == "end" else seg[:n]
    d = pts[-1] - pts[0]
    return d / np.linalg.norm(d)


def step_overs(pieces: list[np.ndarray]) -> list[dict]:
    """Geometry of each junction between along-strike-ordered pieces."""
    out = []
    for i, (a, b) in enumerate(zip(pieces, pieces[1:])):
        u = local_unit(a, "end")
        perp = np.array([-u[1], u[0]])          # +perp is left of travel
        v = b[0] - a[-1]
        across = float(v @ perp)
        sense = "left" if across > 0 else "right"
        out.append({
            "index": i + 1,
            "a_end": a[-1], "b_start": b[0],
            "gap_km": float(v @ u),
            "width_km": abs(across),
            "separation_km": float(np.hypot(*v)),
            "sense": sense,
            "kind": SENSE_KIND[sense],
        })
    return out


def strand_overlaps(rotated: list[np.ndarray], pairs: list[tuple[int, int]]) -> list[dict]:
    """Fault-normal separation where two strands run side by side.

    The junctions between consecutive segments capture only the along-strike
    breaks. The structurally important step-overs on this fault are the
    parallel strands of the Kangding splay, which overlap along strike for
    tens of km with a separation of several km -- geometry a junction-only
    analysis misses entirely.
    """
    out = []
    for i, j in pairs:
        a, b = rotated[i], rotated[j]
        lo = max(a[:, 0].min(), b[:, 0].min())
        hi = min(a[:, 0].max(), b[:, 0].max())
        if hi <= lo:
            continue
        xs = np.linspace(lo, hi, 200)
        sa, sb = np.argsort(a[:, 0]), np.argsort(b[:, 0])
        sep = np.interp(xs, b[sb, 0], b[sb, 1]) - np.interp(xs, a[sa, 0], a[sa, 1])
        sense = "left" if sep.mean() > 0 else "right"
        out.append({
            "pair": (i, j), "x_lo": float(lo), "x_hi": float(hi),
            "overlap_km": float(hi - lo),
            "min_sep_km": float(np.abs(sep).min()),
            "max_sep_km": float(np.abs(sep).max()),
            "sense": sense, "kind": SENSE_KIND[sense],
            "xs": xs, "sep": sep,
        })
    return out


def main() -> None:
    here = Path(__file__).resolve().parent
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--kml", type=Path,
                    default=here / "user_fault_geometry_input" / "xianshuihe_fault_trace.kml")
    ap.add_argument("--output", type=Path, default=here / "xianshuihe_step_overs.png")
    args = ap.parse_args()

    segs = read_kml(args.kml)
    lon0 = np.concatenate([s[:, 0] for s in segs]).mean()
    lat0 = np.concatenate([s[:, 1] for s in segs]).mean()
    local = [to_local_km(s, lon0, lat0) for s in segs]
    centre, theta = principal_frame(np.vstack(local))
    rotated = [rotate(p, centre, theta) for p in local]

    main_idx = FAULT_DECOMPOSITION["ft1 through-going (segs 1,2,0,4,6,7,8)"][0]
    pieces, _ = chain(rotated, main_idx)
    total_km = sum(polyline_length_km(p) for p in pieces)
    threshold = STEP_WIDTH_FRACTION * total_km

    steps = step_overs(pieces)
    counted = [s for s in steps if s["width_km"] > threshold]

    print(f"through-going trace: {total_km:.1f} km over {len(pieces)} segments")
    print(f"Qiao et al. criterion: step-over wider than {STEP_WIDTH_FRACTION:.0%} "
          f"of length = {threshold:.1f} km\n")
    print(f"{'#':>2} {'x (km)':>8} {'gap':>7} {'width':>7} {'sep':>7}  {'sense':<6} "
          f"{'kind':<12} counted")
    for s in steps:
        print(f"{s['index']:>2} {s['a_end'][0]:>8.1f} {s['gap_km']:>7.1f} "
              f"{s['width_km']:>7.1f} {s['separation_km']:>7.1f}  {s['sense']:<6} "
              f"{s['kind']:<12} {'yes' if s in counted else 'no'}")
    # Parallel strands of the Kangding splay. seg3/seg4/seg5 run side by side,
    # so the independent steps are 3->4 and 4->5; 3->5 is the sum of the two
    # and is reported but not counted, to avoid double counting one structure.
    overlaps = strand_overlaps(rotated, [(3, 4), (4, 5), (3, 5)])
    print(f"\nparallel strand overlaps (Kangding splay):")
    print(f"{'pair':>7} {'x_lo':>8} {'x_hi':>8} {'overlap':>8} {'min_sep':>8} "
          f"{'max_sep':>8}  {'kind':<12} counted")
    counted_overlaps = []
    for o in overlaps:
        indep = o["pair"] != (3, 5)
        is_counted = indep and o["max_sep_km"] > threshold
        if is_counted:
            counted_overlaps.append(o)
        note = "yes" if is_counted else ("no (3->4 + 4->5)" if not indep else "no")
        print(f"{o['pair'][0]}-{o['pair'][1]:<5} {o['x_lo']:>8.1f} {o['x_hi']:>8.1f} "
              f"{o['overlap_km']:>8.1f} {o['min_sep_km']:>8.1f} {o['max_sep_km']:>8.1f}  "
              f"{o['kind']:<12} {note}")

    n_total = len(counted) + len(counted_overlaps)
    density = n_total / total_km
    print(f"\ncounted step-overs: {len(counted)} junction + "
          f"{len(counted_overlaps)} strand = {n_total}")
    print(f"step-over density : {n_total}/{total_km:.0f} = {density:.4f} per km")
    print(f"Qiao et al. report: 4/350 = 0.0114 per km (Xianshuihe), "
          f"1/550 = 0.0018 per km (Yushu-Ganzi)")
    if n_total < 4:
        print("\nNOTE: this 1:1M regional compilation resolves fewer step-overs than\n"
              "Qiao et al. count from detailed mapping. Junction offsets here are\n"
              f"0.0-3.4 km, all below the {threshold:.1f} km cut-off; only the splay\n"
              "strands exceed it. Reproducing their maturity metric needs a finer\n"
              "source (Chevalier et al. 2018 / Xu et al. 2016, or 1:250k mapping).")

    # ---- figure -----------------------------------------------------------
    MM = 1.0 / 25.4
    plt.rcParams.update({
        "font.size": 13, "axes.titlesize": 14, "axes.labelsize": 14,
        "xtick.labelsize": 12, "ytick.labelsize": 12, "legend.fontsize": 11.5,
        "savefig.bbox": "tight", "pdf.fonttype": 42, "ps.fonttype": 42,
    })
    fig = plt.figure(figsize=(190 * MM, 230 * MM))
    gs = fig.add_gridspec(3, 1, height_ratios=[1.0, 0.85, 1.15], hspace=0.42)
    ax_a, ax_b, ax_c = (fig.add_subplot(gs[0]), fig.add_subplot(gs[1]),
                        fig.add_subplot(gs[2]))

    # (a) along-strike overview
    for k, p in enumerate(pieces):
        ax_a.plot(p[:, 0], p[:, 1], "-", lw=2.6, color="0.15",
                  label="through-going trace" if k == 0 else None)
    ax_a.plot(rotated[4][:, 0], rotated[4][:, 1], "-", lw=2.4, color="tab:orange",
              label="seg4 strand")
    ax_a.plot(rotated[5][:, 0], rotated[5][:, 1], "-", lw=2.4, color="tab:red",
              label="seg5 strand")
    ax_a.plot(rotated[3][:, 0], rotated[3][:, 1], "-", lw=2.4, color="tab:blue",
              label="seg3 strand")
    for o in counted_overlaps:
        ax_a.axvspan(o["x_lo"], o["x_hi"], color="tab:red", alpha=0.10, lw=0)
    for s_ in steps:
        ax_a.plot(s_["a_end"][0], s_["a_end"][1], "v", ms=9, color="0.35")
        ax_a.annotate(f"{s_['index']}", (s_["a_end"][0], s_["a_end"][1]), fontsize=11,
                      color="0.25", xytext=(0, 11), textcoords="offset points",
                      ha="center")
    ax_a.set_xlabel("Along-strike distance (km)")
    ax_a.set_ylabel("Fault-normal (km)")
    ax_a.grid(alpha=0.25, lw=0.7)
    ax_a.legend(loc="lower right", frameon=False, ncol=2)
    ax_a.text(0.012, 0.96, "(a)", transform=ax_a.transAxes, fontsize=16,
              fontweight="bold", va="top")
    ax_a.set_title("Trace, strands, and junction step-overs (triangles)")

    # (b) separation profiles -- this is the quantitative step-over map
    for o in overlaps:
        indep = o["pair"] != (3, 5)
        ax_b.plot(o["xs"], np.abs(o["sep"]), "-" if indep else "--", lw=2.4,
                  label=f"seg{o['pair'][0]}-seg{o['pair'][1]}"
                        f"{'' if indep else ' (sum of the other two)'}")
    ax_b.axhline(threshold, color="0.35", ls=":", lw=2.0)
    ax_b.annotate(f"{STEP_WIDTH_FRACTION:.0%} of fault length = {threshold:.1f} km",
                  (ax_b.get_xlim()[0], threshold), fontsize=11.5, color="0.3",
                  xytext=(6, 6), textcoords="offset points")
    for s_ in steps:
        ax_b.plot(s_["a_end"][0], s_["width_km"], "v", ms=9, color="0.35")
    ax_b.set_xlabel("Along-strike distance (km)")
    ax_b.set_ylabel("Step width (km)")
    ax_b.grid(alpha=0.25, lw=0.7)
    ax_b.legend(loc="upper right", frameon=False)
    ax_b.text(0.012, 0.96, "(b)", transform=ax_b.transAxes, fontsize=16,
              fontweight="bold", va="top")
    ax_b.set_title("Strand separation vs along-strike distance")

    # (c) map view of the splay zone
    for i_, seg in enumerate(segs):
        ax_c.plot(seg[:, 0], seg[:, 1], "-o", ms=3.5, lw=2.2,
                  color={3: "tab:blue", 4: "tab:orange", 5: "tab:red"}.get(i_, "0.15"))
        if i_ in (3, 4, 5):
            ax_c.annotate(f"seg{i_}", (seg[len(seg) // 2, 0], seg[len(seg) // 2, 1]),
                          fontsize=13, fontweight="bold",
                          color={3: "tab:blue", 4: "tab:orange", 5: "tab:red"}[i_],
                          xytext=(8, 0), textcoords="offset points")
    splay = np.vstack([segs[3], segs[4], segs[5]])
    padx, pady = 0.18, 0.12
    ax_c.set_xlim(splay[:, 0].min() - padx, splay[:, 0].max() + padx)
    ax_c.set_ylim(splay[:, 1].min() - pady, splay[:, 1].max() + pady)
    ax_c.set_aspect(1.0 / np.cos(np.radians(lat0)))
    ax_c.set_xlabel("Longitude ($^\\circ$E)")
    ax_c.set_ylabel("Latitude ($^\\circ$N)")
    ax_c.grid(alpha=0.25, lw=0.7)
    ax_c.text(0.012, 0.97, "(c)", transform=ax_c.transAxes, fontsize=16,
              fontweight="bold", va="top")
    ax_c.set_title("Kangding splay zone, map view")

    fig.suptitle(f"Xianshuihe step-overs: {len(counted)} junction + "
                 f"{len(counted_overlaps)} strand wider than {threshold:.1f} km "
                 f"(density {density:.4f} / km)", fontsize=15, y=0.995)
    fig.savefig(args.output, dpi=300)
    fig.savefig(args.output.with_suffix(".pdf"))
    print(f"\nWrote {args.output}")


if __name__ == "__main__":
    main()

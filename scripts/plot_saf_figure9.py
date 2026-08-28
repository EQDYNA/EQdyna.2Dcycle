#!/usr/bin/env python3
"""Special-event slip/rupture-time panels for an EQdyna.2Dcycle SAF/NSJF run.

Port of the Zenodo MATLAB script
`archive/published/zenodo_software/dunyuliu/dunyuliu-Multicycle_dynamic_SSAF_NSJF-9e645cf/
result/Figure9_Special_Events.m` (do not modify that file; it is the parity
oracle). Reproduces Figure 9 of Liu et al. (2022, JGR-SE, 10.1029/2021JB023420):
a 10-row stack of (1) slip + rupture-time profiles for 8 named characteristic
events (1857 Fort Tejon-like, 1812 Wrightwood-like, Cajon Pass gate-crossing,
single-strand events), (2) a fault-strength (normal-stress-derived) profile
for one of them, and (3) the fault trace map.

The 8 event cycle IDs (absolute, 1-based, matching `ic(1)==1` in the .m) and
their subplot ROW assignment are hardcoded in the .m via if/elseif on the
cycle index -- NOT sorted by id or time. Ported verbatim as EVENT_ROWS below;
row 9 (fault strength) and row 10 (fault geometry) are fixed extras.

Magnitude convention: lockdepth=22e3 m, rig=3500*3500*3000=3.675e10 Pa --
identical to Figure6_Plot_Magnitude_Frequency.m, see plot_saf_figure6.py's
docstring for the rig/lockdepth discrepancy with Figure5 / paleo_site_stats.py
and scripts/plotRuptureDynamics.

Item ported verbatim, not "fixed": rupture times > 500 s are set to NaN
before plotting only (a display clamp for un-ruptured/never-triggered nodes,
not a physical bound); does not affect slip or magnitude.

Usage:
    python3 plot_saf_figure9.py [case_dir] [--output-dir DIR]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from saf_result_utils import (
    fault_node_weights_km,
    load_saf_full_case,
    moment_and_magnitude,
)

LOCKDEPTH_M = 22.0e3
RIG_PA = 3500.0 * 3500.0 * 3000.0
ALPHA = 1.0

MIUS = 0.5  # friction coefficient for Model A, used in the strength panel

RUPTURE_TIME_NAN_CUTOFF_S = 500.0

# (subplot row 1-based, cycle id 1-based) verbatim from the if/elseif chain
# in Figure9_Special_Events.m lines 112-141.
EVENT_ROWS = [
    (1, 856),
    (2, 1041),
    (3, 2047),
    (4, 355),
    (5, 2704),
    (6, 3677),
    (7, 3904),
    (8, 3909),
]
STRENGTH_PANEL_ROW = 9
STRENGTH_PANEL_CYCLE = 2704
GEOMETRY_PANEL_ROW = 10

EVENT_LABELS = {4: "1857-like", 5: "1857-like", 6: "1812-like", 8: "1918-like"}


def compute(case_dir: Path):
    full = load_saf_full_case(case_dir)
    if full.cycle_start != 1:
        sys.exit(f"error: Figure9 porting assumes cycle_start==1 (ic(1)==1 in "
                 f"the .m); got cycle_start={full.cycle_start}.")

    weights = fault_node_weights_km(full.fault_blocks)

    events = {}
    for row, cyc in EVENT_ROWS:
        idx = cyc - full.cycle_start
        if idx < 0 or idx >= full.total_events:
            sys.exit(f"error: cycle {cyc} (row {row}) out of range "
                     f"[1, {full.cycle_end}] for this case")
        slip_row = full.slip_m[idx]
        rupt_row = full.rupture_time_s[idx].copy()
        rupt_row[rupt_row > RUPTURE_TIME_NAN_CUTOFF_S] = np.nan
        mom, mag = moment_and_magnitude(
            slip_row[None, :], weights, RIG_PA, LOCKDEPTH_M, ALPHA
        )
        events[row] = {
            "cycle": cyc, "slip": slip_row, "rupt": rupt_row,
            "mom": float(mom[0]), "mag": float(mag[0]),
        }

    strength_idx = STRENGTH_PANEL_CYCLE - full.cycle_start
    ns_row = full.normal_stress[strength_idx]
    strength = -ns_row * MIUS / 1.0e6  # MPa

    return full, events, strength


def plot(full, events: dict, strength: np.ndarray, output_path: Path) -> None:
    blocks = full.fault_blocks
    fig, axes = plt.subplots(10, 1, figsize=(8, 12), sharex=True)
    colors = ["tab:blue", "tab:red", "k"]
    panel_letters = "abcdefghij"

    for row, cyc in EVENT_ROWS:
        ax = axes[row - 1]
        ev = events[row]
        for b, c in zip(blocks, colors):
            seg_slip = ev["slip"][b.index_start:b.index_stop]
            ax.plot(b.x_km, seg_slip, color=c, linewidth=1)
            ax.fill_between(b.x_km, seg_slip, color=c, alpha=0.3)
        ax.set_xlim(-300, 220)
        ax2 = ax.twinx()
        for b, c in zip(blocks, colors):
            seg_rupt = ev["rupt"][b.index_start:b.index_stop]
            ax2.plot(b.x_km, seg_rupt, color=c, linewidth=2)
        ax2.set_ylim(0, 100)
        if cyc == 2704:
            ax.set_ylabel("Slip (m)")
            ax2.set_ylabel("Time (s)")
        ax.text(0.025, 0.85, f"({panel_letters[row - 1]})", transform=ax.transAxes, fontsize=12)
        ax.text(0.025, 0.55, f"{ev['mag']:.2f}", transform=ax.transAxes, fontsize=12)
        if row in EVENT_LABELS:
            ax.text(0.85, 0.85, EVENT_LABELS[row], transform=ax.transAxes, fontsize=12)

    ax = axes[STRENGTH_PANEL_ROW - 1]
    for b, c in zip(blocks, colors):
        seg_strength = strength[b.index_start:b.index_stop]
        ax.plot(b.x_km, seg_strength, color=c, linewidth=2)
    ax.set_ylabel("Strength (MPa)")
    ax.text(0.025, 0.85, f"({panel_letters[STRENGTH_PANEL_ROW - 1]})",
            transform=ax.transAxes, fontsize=12)

    ax = axes[GEOMETRY_PANEL_ROW - 1]
    for b, c in zip(blocks, colors):
        ax.plot(b.x_km, b.y_km, color=c, linewidth=2)
    ax.set_ylim(-20, 70)
    ax.set_xlim(-300, 220)
    ax.set_xlabel("NW-SE (km)")
    ax.set_ylabel("SW-NE (km)")
    ax.text(0.025, 0.85, f"({panel_letters[GEOMETRY_PANEL_ROW - 1]})",
            transform=ax.transAxes, fontsize=12)

    fig.tight_layout()
    fig.savefig(output_path, dpi=150, facecolor="white")
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("case_dir", nargs="?", default=".", help="Case directory (default: cwd)")
    ap.add_argument("--output-dir", default=None, help="Output directory (default: <case_dir>/aPlots)")
    args = ap.parse_args()

    case_dir = Path(args.case_dir).resolve()
    output_dir = Path(args.output_dir) if args.output_dir else case_dir / "aPlots"
    output_dir.mkdir(parents=True, exist_ok=True)

    full, events, strength = compute(case_dir)
    print(f"case   : {case_dir}")
    print(f"cycles : {full.total_events}  (cycle_start={full.cycle_start}, cycle_end={full.cycle_end})")
    print(f"{'row':>3} {'cycle':>6} {'mag':>6} {'mom (N*m)':>14}")
    for row, cyc in EVENT_ROWS:
        ev = events[row]
        print(f"{row:>3} {cyc:>6} {ev['mag']:>6.2f} {ev['mom']:>14.6e}")

    out = output_dir / "Figure9_SpecialEvents.png"
    plot(full, events, strength, out)
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()

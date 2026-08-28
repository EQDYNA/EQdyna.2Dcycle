#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from saf_result_utils import geologic_rate_profile, load_saf_case, modeled_slip_rate_mm_per_yr

# NOTE (parity audit, see report): geologic_rate_profile()'s segment-boundary constant
# OBSERVED_EQDYNA_X_KM in saf_result_utils.py does not match the reference MATLAB
# (Figure3_Sliprates_UCERF3.m + observed_sliprates.m, which compute these x-positions via
# deg2utm/convert/rotate of the UCERF3 site lon/lat). Running the MATLAB oracle gives
# loc_eqdyna(1:3,1) = [-164.315662175388, 21.621033567916, 129.481611408865] km, vs.
# [-158.576402497732, -69.7459994869569, -17.4212282973009] currently hardcoded --
# shifting the Carrizo/Mojave-S/San-Bernardino-S segment boundaries along the SSAF trace
# by 10s of km and misclassifying ~35% of fault-3 nodes into the wrong observed-rate
# bucket. saf_result_utils.py is owned by another agent; this file cannot fix it, so
# flagging it here where the wrong values are consumed. Do not treat this plot's
# "observed" band as parity-verified until that constant is corrected upstream.


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot a Figure 3 style SAF slip-rate comparison directly from an EQdyna case."
    )
    parser.add_argument("case_dir", nargs="?", default="work/safpub", help="Case directory to analyze")
    parser.add_argument(
        "--output",
        default=None,
        help="Output image path. Defaults to <case_dir>/aPlots/Figure3_sliprates.png",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    case_dir = Path(args.case_dir)
    output = Path(args.output) if args.output else case_dir / "aPlots" / "Figure3_sliprates.png"
    output.parent.mkdir(parents=True, exist_ok=True)

    case_data = load_saf_case(case_dir)
    modeled_rate = modeled_slip_rate_mm_per_yr(case_data)

    fig = plt.figure(figsize=(7, 5), constrained_layout=False)
    ax_rate = fig.add_axes([0.1, 0.45, 0.8, 0.45])
    ax_geom = fig.add_axes([0.1, 0.2, 0.8, 0.2], sharex=ax_rate)

    first_model = True
    fault_colors = {"sjfn": "#1f77b4", "sjfs": "#d62728", "ssaf": "#222222"}
    fault_labels = {"sjfn": "NSJF north", "sjfs": "NSJF south", "ssaf": "SSAF"}

    for block in case_data.fault_blocks:
        rates = geologic_rate_profile(block)
        block_modeled = modeled_rate[block.index_start:block.index_stop]
        color = fault_colors[block.name]

        ax_rate.fill_between(
            block.x_km,
            rates[:, 0],
            rates[:, 1],
            color=color,
            alpha=0.12,
            linewidth=0.0,
        )
        ax_rate.plot(
            block.x_km,
            rates[:, 2],
            color=color,
            linestyle="--",
            linewidth=1.4,
            label=f"{fault_labels[block.name]} observed" if first_model else None,
        )
        ax_rate.plot(
            block.x_km,
            block_modeled,
            color=color,
            linewidth=2.0,
            label=f"{fault_labels[block.name]} modeled",
        )

        ax_geom.plot(block.x_km, block.y_km, color=color, linewidth=2.0)
        first_model = False

    ax_rate.set_ylabel("Slip rate (mm/yr)")
    ax_rate.set_xlim(-300.0, 220.0)
    ax_rate.set_ylim(bottom=0.0)
    ax_rate.tick_params(axis="x", labelbottom=False)
    ax_rate.set_title(f"SAF long-term slip rates: {case_dir.name}")
    ax_rate.legend(loc="upper right", fontsize=8, ncol=2)

    ax_geom.set_xlabel("NW-SE (km)")
    ax_geom.set_ylabel("SW-NE (km)")
    ax_geom.set_xlim(-300.0, 220.0)
    ax_geom.set_ylim(-20.0, 70.0)
    ax_geom.set_aspect("equal", adjustable="box")

    fig.patch.set_facecolor("white")
    fig.savefig(output, dpi=200)
    plt.close(fig)

    print(f"Wrote {output}")


if __name__ == "__main__":
    main()

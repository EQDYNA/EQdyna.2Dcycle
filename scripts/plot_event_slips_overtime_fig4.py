#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from saf_result_utils import load_saf_case, select_time_window


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot Figure 4/7/8 style SAF slip distributions directly from an EQdyna case."
    )
    parser.add_argument("case_dir", nargs="?", default=".", help="Case directory to analyze (default: cwd)")
    parser.add_argument(
        "--tstart",
        nargs="+",
        type=float,
        default=None,
        help="Window start times in kyr. If omitted, auto-populate 0, 3, 6, ... up to total simulated kyr.",
    )
    parser.add_argument("--duration", type=float, default=3.0, help="Window duration in kyr (default 3)")
    parser.add_argument("--threshold", type=float, default=1.0, help="Minimum event max-slip in meters")
    parser.add_argument("--scale", type=float, default=30.0, help="Slip scaling factor for filled polygons")
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Output directory. Defaults to <case_dir>/aPlots",
    )
    return parser.parse_args()


def _auto_tstart_kyr(case_dir: Path, duration_kyr: float) -> list[float]:
    """Auto-populate tstart=[0, duration, 2·duration, ...] covering all simulated years.

    Reads interval.txt<icstart> (sum of interseismic durations in years) from the case
    dir. One window if <= duration kyr of data. Ceil division for longer runs.
    """
    import math
    total_yr = 0.0
    # take the first interval.txt<icstart> we can find
    candidates = sorted(case_dir.glob("interval.txt*")) + sorted((case_dir / "aRawSimuData").glob("interval.txt*")) if (case_dir / "aRawSimuData").exists() else sorted(case_dir.glob("interval.txt*"))
    for f in candidates:
        try:
            vals = np.loadtxt(f)
            total_yr += float(np.atleast_1d(vals).sum())
        except Exception:
            continue
    total_kyr = total_yr / 1000.0
    nwin = max(1, math.ceil(total_kyr / duration_kyr))
    return [i * duration_kyr for i in range(nwin)]


def plot_window(
    case_data,
    case_dir: Path,
    output_dir: Path,
    tstart_kyr: float,
    duration_kyr: float,
    threshold_m: float,
    scale: float,
) -> Path:
    tend_kyr = tstart_kyr + duration_kyr
    event_indices, cycle_ids = select_time_window(case_data, tstart_kyr, tend_kyr)

    fig, ax = plt.subplots(figsize=(7, 7))
    saf_colors = {"sjfn": "#1f77b4", "sjfs": "#d62728", "ssaf": "#222222"}
    cmap = plt.get_cmap("tab10")
    colors = {
        block.name: saf_colors.get(block.name, cmap(i % 10))
        for i, block in enumerate(case_data.fault_blocks)
    }
    x_label_anchor = max(block.x_km.max() for block in case_data.fault_blocks) - 35.0
    shown = 0

    for index, cycle_id in zip(event_indices, cycle_ids):
        slip = case_data.slips_m[index]
        max_slip = float(np.max(slip))
        if max_slip <= threshold_m:
            continue
        shown += 1
        time_kyr = case_data.event_times_kyr[index]

        for block in case_data.fault_blocks:
            block_slip = slip[block.index_start:block.index_stop]
            xx = np.concatenate([block.x_km, block.x_km[::-1]])
            yy = np.concatenate(
                [
                    block_slip / scale + time_kyr,
                    np.full(block.count, time_kyr)[::-1],
                ]
            )
            ax.fill(xx, yy, color=colors[block.name], alpha=0.28, linewidth=0.0)

        ax.text(
            x_label_anchor + 10.0 * (shown % 3),
            time_kyr,
            str(int(cycle_id)),
            fontsize=6,
            fontweight="bold",
        )

    baseline_y = tstart_kyr - 0.4
    for block in case_data.fault_blocks:
        ax.plot(block.x_km, block.y_km / 150.0 + baseline_y, color=colors[block.name], linewidth=2.0)

    ymin = tstart_kyr - 0.5
    ymax = tend_kyr + 0.2
    x_min = min(block.x_km.min() for block in case_data.fault_blocks) - 20.0
    x_max = max(block.x_km.max() for block in case_data.fault_blocks) + 20.0
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(ymin, ymax)
    ax.set_xlabel("NW-SE (km)")
    ax.set_ylabel("Time (kyrs)")
    ax.set_title(f"SAF slip distributions: {case_dir.name}, {tstart_kyr:g}-{tend_kyr:g} kyr")

    scale_x = x_max - 40.0
    ax.plot([scale_x, scale_x], [ymin + 0.05, ymin + 0.05 + 5.0 / scale], color="k", linewidth=3.0)
    ax.text(scale_x + 5.0, ymin + 0.05 + 2.5 / scale, "slip = 5 m", va="center")

    fig.patch.set_facecolor("white")
    output = output_dir / f"Figure4_7_8_slipdist_{tstart_kyr:g}_{tend_kyr:g}kyr.png"
    fig.savefig(output, dpi=200)
    plt.close(fig)
    return output


def main() -> None:
    args = parse_args()
    case_dir = Path(args.case_dir)
    output_dir = Path(args.output_dir) if args.output_dir else case_dir / "aPlots"
    output_dir.mkdir(parents=True, exist_ok=True)
    case_data = load_saf_case(case_dir)

    tstart_list = args.tstart if args.tstart is not None else _auto_tstart_kyr(case_dir, args.duration)
    print(f"tstart windows (kyr): {tstart_list}  duration={args.duration} kyr")

    outputs = []
    for tstart_kyr in tstart_list:
        output = plot_window(
            case_data=case_data,
            case_dir=case_dir,
            output_dir=output_dir,
            tstart_kyr=tstart_kyr,
            duration_kyr=args.duration,
            threshold_m=args.threshold,
            scale=args.scale,
        )
        outputs.append(output)

    for output in outputs:
        print(f"Wrote {output}")


if __name__ == "__main__":
    main()

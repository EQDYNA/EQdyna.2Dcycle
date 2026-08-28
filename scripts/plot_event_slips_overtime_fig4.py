#!/usr/bin/env python3

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from saf_result_utils import load_saf_case, select_time_window


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot Figure 4/7/8 style slip distributions directly from an EQdyna case."
    )
    parser.add_argument("case_dir", nargs="?", default=".", help="Case directory to analyze (default: cwd)")
    parser.add_argument(
        "--tstart",
        nargs="+",
        type=float,
        default=None,
        help=(
            "Window start times in kyr. If omitted, auto-populate duration, 2*duration, ... "
            "up to total simulated kyr, skipping the initial [0, duration) window (see "
            "_auto_tstart_kyr docstring: this mirrors the paper's own Fig 4 panel choice, "
            "tstart=[3,6,9,12] kyr for Model A / duration=3kyr, which always skips the first "
            "window)."
        ),
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


def _auto_tstart_kyr(case_data, duration_kyr: float) -> list[float]:
    """Auto-populate tstart windows covering the simulated run, skipping the first window.

    The reference MATLAB (Figure4_7_8_Plot_Slip_distributions4sequence.m) never plots
    from t=0: its header comment hardcodes the paper's own panel choices --
    tstart=[3,6,9,12] kyr for Model A at duration=3kyr (covering 3-15 kyr of the 15.167
    kyr run), tstart=[4,6] for Model B at duration=2kyr, tstart=[6,9] for Model C at
    duration=3kyr. In each case the *first* duration-length window is skipped (spin-up /
    less representative early cycles) and the ones after are shown back-to-back.

    This reproduces that convention as a formula instead of hardcoding one model's
    numbers: with n_full = floor(total_kyr / duration_kyr) complete windows available,
    return duration_kyr*[1, 2, ..., n_full-1] (i.e. skip window 0, keep the rest). For
    n_full <= 1 (not even one full window's worth of extra data after the burn-in) fall
    back to a single window starting at 0, since there is nothing to skip past on a short
    run. total_kyr is taken from case_data.event_times_kyr, already loaded by
    load_saf_case with correct icstart-tag ordering, rather than re-globbing interval.txt*
    (which double-counted files present in both case_dir and case_dir/aRawSimuData).
    """
    import math
    total_kyr = float(case_data.event_times_kyr[-1] - case_data.event_times_kyr[0])
    n_full = int(math.floor(total_kyr / duration_kyr))
    if n_full <= 1:
        return [0.0]
    return [i * duration_kyr for i in range(1, n_full)]


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
    # Label anti-collision: the reference MATLAB staggers labels across only 3 x-slots
    # (text(152 + mod(ntag,3)*15, ...)), which visibly overlaps whenever events cluster
    # in time (seen on a real run: several events <0.01 kyr apart in the same slot land
    # on top of each other). This is cosmetic-only -- it changes label placement, never
    # which events are shown, their order, or the slip polygons -- so it does not touch
    # parity. Use more slots (6) and nudge close-in-time labels within the same slot
    # apart vertically so they no longer overprint.
    n_slots = 8
    slot_spacing = 12.0
    last_label_y = [None] * n_slots

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

        slot = shown % n_slots
        min_gap_kyr = 0.025 * duration_kyr
        label_y = time_kyr
        if last_label_y[slot] is not None and label_y - last_label_y[slot] < min_gap_kyr:
            label_y = last_label_y[slot] + min_gap_kyr
        last_label_y[slot] = label_y
        ax.text(
            x_label_anchor + slot_spacing * slot,
            label_y,
            str(int(cycle_id)),
            fontsize=5,
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
    # case_dir.name is '' for "." or a trailing-slash path (Path(".").name == ""),
    # which rendered as a stray ", " in the title. Resolve to an absolute path first
    # so cwd-relative invocations (the documented default) still get a real name.
    case_label = case_dir.resolve().name
    # The script is not SAF-specific -- it reads whatever fault blocks the case
    # has, and is used for Xianshuihe and others.
    title = f"Slip distributions: {tstart_kyr:g}-{tend_kyr:g} kyr"
    if case_label:
        title = f"Slip distributions: {case_label}, {tstart_kyr:g}-{tend_kyr:g} kyr"
    ax.set_title(title)

    scale_x = x_max - 40.0
    ax.plot([scale_x, scale_x], [ymin + 0.05, ymin + 0.05 + 5.0 / scale], color="k", linewidth=3.0)
    ax.text(scale_x + 5.0, ymin + 0.05 + 2.5 / scale, "slip = 5 m", va="center")

    fig.patch.set_facecolor("white")
    output = output_dir / f"Figure4_7_8_slipdist_{tstart_kyr:g}_{tend_kyr:g}kyr.png"
    fig.savefig(output, dpi=200)
    plt.close(fig)
    return output, shown


def main() -> None:
    args = parse_args()
    case_dir = Path(args.case_dir)
    output_dir = Path(args.output_dir) if args.output_dir else case_dir / "aPlots"
    output_dir.mkdir(parents=True, exist_ok=True)
    case_data = load_saf_case(case_dir)

    tstart_list = args.tstart if args.tstart is not None else _auto_tstart_kyr(case_data, args.duration)
    print(f"tstart windows (kyr): {tstart_list}  duration={args.duration} kyr")

    # Event times start the clock at the FIRST event, so the last event sits at
    # sum(intervals) - intervals[0]. On a run whose first interseismic period is
    # long, that is well short of the total simulated time, and a window chosen
    # from the total lands past the end of the data.
    t_max = float(case_data.event_times_kyr[-1]) if len(case_data.event_times_kyr) else 0.0
    print(f"events span 0 to {t_max:.2f} kyr ({len(case_data.event_times_kyr)} events)")

    outputs, empty = [], []
    for tstart_kyr in tstart_list:
        output, shown = plot_window(
            case_data=case_data,
            case_dir=case_dir,
            output_dir=output_dir,
            tstart_kyr=tstart_kyr,
            duration_kyr=args.duration,
            threshold_m=args.threshold,
            scale=args.scale,
        )
        outputs.append((output, shown))
        if shown == 0:
            empty.append((tstart_kyr, output))

    for output, shown in outputs:
        print(f"Wrote {output}  ({shown} events drawn)")

    # R29: never render an empty result without saying so. This figure used to
    # save axes, fault traces and a scale bar with no events on them and exit 0,
    # so a window past the end of the data -- or a threshold above every event's
    # slip -- looked like a finished plot.
    if empty:
        print()
        for tstart_kyr, output in empty:
            beyond = tstart_kyr > t_max
            why = (f"window starts at {tstart_kyr:g} kyr but the events end at {t_max:.2f} kyr"
                   if beyond else
                   f"no event in the window exceeds --threshold {args.threshold} m")
            print(f"EMPTY: {output.name} has no events -- {why}")
        print("Nothing was plotted. Pick a window inside the event range, or lower "
              "--threshold.")
        sys.exit(1)


if __name__ == "__main__":
    main()

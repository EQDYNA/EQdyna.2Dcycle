#!/usr/bin/env python3
"""Magnitude-frequency statistics for an EQdyna.2Dcycle SAF/NSJF run.

Port of the Zenodo MATLAB script
`archive/published/zenodo_software/dunyuliu/dunyuliu-Multicycle_dynamic_SSAF_NSJF-9e645cf/
result/Figure6_Plot_Magnitude_Frequency.m` (do not modify that file; it is the
parity oracle). Reproduces Figure 6 of Liu et al. (2022, JGR-SE,
10.1029/2021JB023420): (a) cumulative moment release, (b) magnitude histogram
+ CDF of all events, (c) magnitudes > 6.6 vs prehistoric/historical Scharer &
Yule (2020) reconstructions.

Magnitude convention (verbatim from the MATLAB, NOT the same as
Figure5_Plot_Recurrene_Stats.m / paleo_site_stats.py):
    lockdepth = 22e3 m,  rig (mu) = 3500*3500*3000 = 3.675e10 Pa
Figure5 instead uses lockdepth=14e3, rig=3464*3464*2670=3.2039e10 Pa, but
Figure5's script never actually emits a magnitude (Table 2 has none), so
those constants are dead weight there. Figure6 AND Figure9 (which computes
the paper's named-event magnitudes, e.g. Mw 7.9 for the 1857 Fort Tejon
analog) both use the 22 km / 3.675e10 Pa pair -- that is the paper's actual
published magnitude convention. `scripts/plotRuptureDynamics` currently uses
22 km (matches) but rig=2670*3464**2=3.2039e10 Pa (matches Figure5, NOT
Figure6/9), which understates every magnitude by
2/3*log10(3.675e10/3.2039e10) = 0.040 Mw relative to the paper's own Fig 6/9
convention. Not changed here (out of scope / owned by another agent); see
the port report.

Moment/magnitude kernel: see saf_result_utils.fault_node_weights_km and
moment_and_magnitude, shared verbatim with plot_saf_figure9.py (both port
the same mom/len accumulation loop from the .m files).

Two items ported verbatim, not "fixed":
  * Panel (c) filters `magrec(i) > 6.6`, not 6.5 as the paper text/figure
    caption implies ("magnitudes >6.5"). This is what the .m literally does;
    ported as-is.
  * `slipthreshold = 0.5` is assigned in the MATLAB (line 6) but never used
    anywhere in the rest of the script -- dead code, not ported.

Scharer & Yule (2020) Table S2 reconstructions (magrec_obs1, magrec_obs2,
count_obs1, count_obs2) are hardcoded directly in the .m (comment: "Recontructed
following Kate's 10/09/2021 email figure", 10/13/2021) -- there is no
underlying data file in the archive. Carried across verbatim below with the
same citation, since there is nothing else to port from.

Usage:
    python3 plot_saf_figure6.py [case_dir] [--output-dir DIR]
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

# --- constants, verbatim from Figure6_Plot_Magnitude_Frequency.m -----------
LOCKDEPTH_M = 22.0e3
RIG_PA = 3500.0 * 3500.0 * 3000.0
ALPHA = 1.0

MAG_BIN_EDGES_ALL = np.linspace(5.5, 8.5, 31)          # magrange
MAG_BIN_CENTERS_ALL = np.linspace(5.55, 8.45, 30)      # magrangeplot

MAG_THRESHOLD_C = 6.6  # magrec(i) > 6.6, verbatim (not 6.5)
MAG_BIN_EDGES_C = np.linspace(6.5, 8.5, 21)             # magrange1
MAG_BIN_CENTERS_C = np.linspace(6.55, 8.45, 20) - 0.03  # magrangeplot1

# Scharer & Yule (2020), Table S2 reconstruction as hardcoded in
# Figure6_Plot_Magnitude_Frequency.m lines 251-268 ("Recontructed following
# Kate's 10/09/2021 email figure", comment dated 10/13/2021). No source data
# file exists in the archive; this is the entirety of what was published.
MAG_OBS_BIN_LEFT = np.linspace(6.52, 8.02, 16)   # mag_obs_plot1
MAG_OBS_BIN_RIGHT = MAG_OBS_BIN_LEFT + 0.01      # mag_obs_plot2
COUNT_OBS_HISTORICAL = np.array([1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0])  # count_obs1
COUNT_OBS_PALEO = np.array([0, 0, 0, 0, 0, 2, 7, 3, 4, 7, 4, 2, 4, 0, 1, 0])       # count_obs2
MAGREC_OBS_HISTORICAL = np.array([6.4, 6.5, 6.9, 7.2, 7.3, 7.5, 7.8])  # magrec_obs1
MAGREC_OBS_PALEO = np.array(  # magrec_obs2
    [7.0, 7.0,
     7.1, 7.1, 7.1, 7.1, 7.1, 7.1, 7.1,
     7.2, 7.2, 7.2,
     7.3, 7.3, 7.3, 7.3,
     7.4, 7.4, 7.4, 7.4, 7.4, 7.4, 7.4,
     7.5, 7.5, 7.5, 7.5,
     7.6, 7.6,
     7.7, 7.7, 7.7, 7.7,
     7.9]
)


def histogram_matlab(values: np.ndarray, edges: np.ndarray) -> np.ndarray:
    """Reproduce the .m's `for k=2:nm: if v>=edges(k-1) && v<edges(k)` binning
    -- right-open intervals [edges[k-1], edges[k]), values outside
    [edges[0], edges[-1]) are simply not counted (no clamping to end bins)."""
    idx = np.digitize(values, edges, right=False)  # edges[i-1] <= v < edges[i] -> idx=i
    counts = np.zeros(edges.size - 1, dtype=int)
    valid = (idx >= 1) & (idx <= edges.size - 1)
    np.add.at(counts, idx[valid] - 1, 1)
    return counts


def ecdf(values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Empirical CDF matching MATLAB cdfplot's stairstep curve."""
    x = np.sort(values)
    y = np.arange(1, x.size + 1) / x.size
    return x, y


def compute(case_dir: Path):
    full = load_saf_full_case(case_dir)
    if full.cycle_start != 1:
        sys.exit(f"error: Figure6 porting assumes cycle_start==1 (ic(1)==1 in "
                 f"the .m); got cycle_start={full.cycle_start}. Indexing "
                 f"would need re-deriving.")

    weights = fault_node_weights_km(full.fault_blocks)
    mom, mag = moment_and_magnitude(full.slip_m, weights, RIG_PA, LOCKDEPTH_M, ALPHA)
    moaccum = np.cumsum(mom)
    timestamp = full.event_times_kyr
    maxsliprec = full.slip_m.max(axis=1)

    recurrence = np.diff(timestamp)  # kyr, all ic(2)-1 system-wide intervals
    meanrecurr = float(recurrence.mean())
    stand = float(recurrence.std(ddof=1))
    cov = stand / meanrecurr
    meanslip = float(maxsliprec.mean())
    stdslip = float(maxsliprec.std(ddof=1))

    count_all = histogram_matlab(mag, MAG_BIN_EDGES_ALL)
    n_total = full.cycle_end  # ic(2), matches the .m's denominator exactly

    magrec1 = mag[mag > MAG_THRESHOLD_C]
    count_c = histogram_matlab(magrec1, MAG_BIN_EDGES_C)

    return {
        "full": full, "mom": mom, "mag": mag, "moaccum": moaccum,
        "timestamp": timestamp, "recurrence": recurrence,
        "meanrecurr": meanrecurr, "stand": stand, "cov": cov,
        "meanslip": meanslip, "stdslip": stdslip,
        "count_all": count_all, "n_total": n_total,
        "magrec1": magrec1, "count_c": count_c,
    }


def plot(results: dict, output_path: Path) -> None:
    mag = results["mag"]
    magrec1 = results["magrec1"]

    fig, axes = plt.subplots(1, 3, figsize=(10, 4))

    ax = axes[0]
    ax.plot(results["timestamp"], results["moaccum"])
    ax.set_xlabel("Thousand Years")
    ax.set_ylabel("Cumulative Moment Release (N*m)")
    ax.text(0.025, 0.9, "(a)", transform=ax.transAxes, fontsize=12)

    ax = axes[1]
    ax.bar(MAG_BIN_CENTERS_ALL, results["count_all"], width=0.1, color="tab:blue")
    ax.set_ylabel("Number of Events")
    ax2 = ax.twinx()
    x, y = ecdf(mag)
    ax2.step(x, y, where="post", color="tab:blue")
    ax2.set_xlim(5, 8.5)
    ax.set_xlabel("Magnitude")
    ax2.set_ylabel("Frequency")
    ax.text(0.025, 0.9, "(b)", transform=ax.transAxes, fontsize=12)

    ax = axes[2]
    ax.bar(MAG_BIN_CENTERS_C, results["count_c"], width=0.1, color="tab:blue", alpha=0.3, label="model")
    ax.bar(MAG_OBS_BIN_LEFT, COUNT_OBS_PALEO, width=0.1, color="tab:red", alpha=0.3, label="paleo obs.")
    ax.bar(MAG_OBS_BIN_RIGHT, COUNT_OBS_HISTORICAL, width=0.1, color="k", alpha=0.3, label="historical obs.")
    ax.set_ylabel("Number of Events")
    ax2 = ax.twinx()
    x, y = ecdf(magrec1)
    ax2.step(x, y, where="post", color="tab:blue", linewidth=2)
    x, y = ecdf(MAGREC_OBS_PALEO)
    ax2.step(x, y, where="post", color="tab:red", linewidth=2)
    x, y = ecdf(MAGREC_OBS_HISTORICAL)
    ax2.step(x, y, where="post", color="k", linewidth=2)
    ax2.set_xlim(6.0, 8.5)
    ax.set_xlabel("Magnitude")
    ax2.set_ylabel("Frequency")
    ax.legend(fontsize=7, loc="lower right")
    ax.text(0.025, 0.9, "(c)", transform=ax.transAxes, fontsize=12)

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

    r = compute(case_dir)
    full = r["full"]
    print(f"case      : {case_dir}")
    print(f"cycles    : {full.total_events}  (cycle_start={full.cycle_start}, "
          f"cycle_end={full.cycle_end})")
    print(f"meanrecurr = {r['meanrecurr']:.6f} kyr")
    print(f"stand      = {r['stand']:.6f} kyr")
    print(f"COV        = {r['cov']:.6f}")
    print(f"meanslip   = {r['meanslip']:.6f} m")
    print(f"stdslip    = {r['stdslip']:.6f} m")
    print(f"n_total (ic(2)) = {r['n_total']}")
    print(f"count (all, sum)  = {int(r['count_all'].sum())}")
    print(f"count (mag>6.6, sum) = {int(r['count_c'].sum())}")

    out = output_dir / "Figure6_MagFreq.png"
    plot(r, out)
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()

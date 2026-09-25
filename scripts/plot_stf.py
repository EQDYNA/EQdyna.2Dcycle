#!/usr/bin/env python3
"""Plot the source time functions of one event from stf.bin output.

Four panels:
  (a) moment-rate function of the whole event
  (b) record section: slip-rate time functions at ~40 nodes along strike,
      rupture front and nucleation overlaid
  (c) slip-rate time functions at the nucleation node and the peak-slip nodes
  (d) final slip along strike, per fault

Usage:
    plot_stf.py <case_dir|stf.bin*> <eqid> [<eqid> ...] [--out DIR]
"""

from __future__ import annotations

import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from stf_read import STFCase  # noqa: E402

MM = 1.0 / 25.4
STYLE = {"font.size": 10, "axes.titlesize": 10.5, "axes.labelsize": 10,
         "xtick.labelsize": 9, "ytick.labelsize": 9, "legend.fontsize": 8.5,
         "savefig.bbox": "tight"}
FAULT_COLOURS = plt.cm.tab10(np.arange(10))      # fixed order, as the other figures
N_TRACES = 40                                     # traces in the record section
TRACE_GAIN = 2.0                                  # peak trace width, in trace spacings


def plot_event(case: STFCase, eqid: int, out_dir: str) -> str:
    ev = case.event(eqid)
    t, mdot = case.moment_rate(eqid)
    x_km = ev["x"] / 1e3
    trap = getattr(np, "trapezoid", None) or np.trapz
    m0 = trap(mdot, t) if len(t) > 1 else 0.0
    mw = 2.0 / 3.0 * np.log10(m0 * 1e7) - 10.7 if m0 > 0 else float("nan")

    plt.rcParams.update(STYLE)
    fig, ax = plt.subplots(4, 1, figsize=(170 * MM, 270 * MM), layout="constrained",
                           gridspec_kw={"height_ratios": [1, 2.6, 1.2, 1]})

    # (a) moment rate
    ax[0].plot(t, mdot, color="0.15", lw=1.6)
    ax[0].set_ylabel(r"$\dot M_0$ (N m/s)")
    ax[0].set_title(f"(a) Event {eqid}: Mw {mw:.2f}, t = {ev['time_yr']:.1f} yr in the "
                    f"sequence, {ev['nout']} nodes, {ev['t_end']:.1f} s", loc="left")
    ax[0].set_xlim(0, t[-1])
    ax[0].set_xlabel("time (s)")

    faults = np.unique(ev["fault"])
    rt = ev["rupture_time"]
    ok = rt < 999
    xlim = (float(x_km.min()), float(x_km.max()))

    # (b) record section: each trace is one node's slip-rate function, drawn at
    # its along-strike position with amplitude to the right, time upward
    rs = ax[1]
    order_all = np.argsort(x_km)
    step = max(1, len(order_all) // N_TRACES)
    sel = order_all[::step]
    spacing = (x_km.max() - x_km.min()) / max(len(sel), 1)
    gain = TRACE_GAIN * spacing / max(float(ev["slip_rate"].max()), 1e-9)
    for i in sel:
        tr = x_km[i] + ev["slip_rate"][i] * gain
        rs.fill_betweenx(t, x_km[i], tr, color="0.25", lw=0)
        rs.plot(tr, t, color="0.05", lw=0.4)
    rs.plot(x_km[ok], rt[ok], ".", ms=1.2, color="C3", label="rupture front")
    rs.plot(ev["nuc_x"] / 1e3, 0, marker="*", ms=11, color="C3", ls="none",
            label="nucleation")
    rs.plot([], [], color="0.05", lw=1,
            label=f"1 m/s = {gain:.2f} km")
    rs.set_ylim(0, t[-1])
    rs.set_ylabel("time (s)")
    rs.set_xlabel("along-strike x (km)")
    rs.set_title(f"(b) Record section: slip-rate functions at every {step}th node", loc="left")
    rs.legend(loc="upper right", frameon=True)
    rs.set_xlim(xlim)

    # (c) slip-rate functions at a few informative nodes
    final = ev["slip"][:, -1]
    pick = [int(np.argmin(np.abs(ev["node"] - (ev["nuc_node"] - 1))))]
    pick += [int(i) for i in np.argsort(final)[::-1][:3] if int(i) not in pick]
    greys = ["0.05", "0.35", "0.55", "0.72"]
    for c, i in zip(greys, pick):
        lab = "nucleation node" if i == pick[0] else f"x = {x_km[i]:.1f} km"
        ax[2].plot(t, ev["slip_rate"][i], color=c, lw=1.4, label=lab)
    ax[2].set_ylabel("slip rate (m/s)")
    ax[2].set_xlim(0, t[-1])
    ax[2].set_xlabel("time (s)")
    ax[2].set_title("(c) Slip-rate time functions", loc="left")
    ax[2].legend(loc="upper right", ncol=2, frameon=False)

    # (d) final slip along strike
    for f in faults:
        m = ev["fault"] == f
        order = np.argsort(x_km[m])
        ax[3].plot(x_km[m][order], final[m][order], lw=1.8,
                   color=FAULT_COLOURS[(int(f) - 1) % 10], label=f"ft{int(f)}")
    ax[3].set_ylabel("slip (m)")
    ax[3].set_xlabel("along-strike x (km)")
    ax[3].set_title("(d) Final slip", loc="left")
    if len(faults) > 1:
        ax[3].legend(loc="upper right", ncol=min(len(faults), 7), frameon=False)
    ax[3].set_xlim(xlim)
    for a in ax:
        a.grid(alpha=0.25, lw=0.5)

    os.makedirs(out_dir, exist_ok=True)
    out = os.path.join(out_dir, f"stf_event{eqid}.png")
    fig.savefig(out, dpi=200)
    plt.close(fig)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("path")
    ap.add_argument("eqid", type=int, nargs="+")
    ap.add_argument("--out", default=None, help="output directory (default <case>/aPlots)")
    a = ap.parse_args()
    case = STFCase(a.path)
    base = a.path if os.path.isdir(a.path) else os.path.dirname(os.path.abspath(a.path))
    out_dir = a.out or os.path.join(base, "aPlots")
    missing = [e for e in a.eqid if e not in case.owner]
    if missing:
        sys.exit(f"event(s) not in the STF file: {missing}; available {case.eqids[0]}..{case.eqids[-1]}")
    for e in a.eqid:
        if case.meta(e)["nout"] == 0:
            print(f"event {e}: no node exceeded the slip-rate threshold; nothing to plot")
            continue
        print(f"wrote {plot_event(case, e, out_dir)}")


if __name__ == "__main__":
    main()

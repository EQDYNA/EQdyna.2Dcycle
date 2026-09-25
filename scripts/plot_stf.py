#!/usr/bin/env python3
"""Plot the source time functions of one event from stf.bin output.

Five panels:
  (a) moment-rate function of the whole event
  (b) slip rate in space and time (along-strike x vs time, log colour),
      rupture front overlaid
  (c) record section: slip-rate time functions at ~40 nodes along strike
  (d) slip-rate time functions at the nucleation node and the peak-slip nodes
  (e) final slip along strike, per fault

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
from matplotlib.colors import LogNorm

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from stf_read import STFCase  # noqa: E402

MM = 1.0 / 25.4
STYLE = {"font.size": 10, "axes.titlesize": 10.5, "axes.labelsize": 10,
         "xtick.labelsize": 9, "ytick.labelsize": 9, "legend.fontsize": 8.5,
         "savefig.bbox": "tight"}
FAULT_COLOURS = plt.cm.tab10(np.arange(10))      # fixed order, as the other figures
VMIN_PLOT = 1.0e-3                                # m/s, bottom of the colour scale
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
    fig, ax = plt.subplots(5, 1, figsize=(170 * MM, 320 * MM), layout="constrained",
                           gridspec_kw={"height_ratios": [1, 2.2, 2.2, 1.2, 1]})

    # (a) moment rate
    ax[0].plot(t, mdot, color="0.15", lw=1.6)
    ax[0].set_ylabel(r"$\dot M_0$ (N m/s)")
    ax[0].set_title(f"(a) Event {eqid}: Mw {mw:.2f}, t = {ev['time_yr']:.1f} yr in the "
                    f"sequence, {ev['nout']} nodes, {ev['t_end']:.1f} s", loc="left")
    ax[0].set_xlim(0, t[-1])
    ax[0].set_xlabel("time (s)")

    # (b) space-time slip rate, one strip per fault so overlapping strands stay apart
    faults = np.unique(ev["fault"])
    v = np.maximum(ev["slip_rate"], VMIN_PLOT * 0.5)
    vmax = max(float(v.max()), VMIN_PLOT * 10)
    for f in faults:
        m = ev["fault"] == f
        order = np.argsort(x_km[m])
        xs, vs = x_km[m][order], v[m][order]
        if xs.size < 2:
            continue
        pc = ax[1].pcolormesh(xs, t, vs.T, shading="nearest", cmap="Oranges",
                              norm=LogNorm(VMIN_PLOT, vmax), rasterized=True)
    rt = ev["rupture_time"]
    ok = rt < 999
    ax[1].plot(x_km[ok], rt[ok], ".", ms=1.5, color="0.1", label="rupture front")
    ax[1].plot(ev["nuc_x"] / 1e3, 0, marker="*", ms=11, color="0.1", ls="none",
               label="nucleation")
    ax[1].set_ylim(0, t[-1])
    ax[1].set_ylabel("time (s)")
    ax[1].set_xlabel("along-strike x (km)")
    ax[1].set_title("(b) Slip rate", loc="left")
    ax[1].legend(loc="upper right", frameon=True)
    cb = fig.colorbar(pc, ax=ax[1], pad=0.01, fraction=0.03)
    cb.set_label("slip rate (m/s)")

    # (c) record section: each trace is one node's slip-rate function, drawn at
    # its along-strike position with amplitude to the right, time upward
    rs = ax[2]
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
    rs.set_title(f"(c) Record section: slip-rate functions at every {step}th node", loc="left")
    rs.legend(loc="upper right", frameon=True)
    rs.set_xlim(ax[1].get_xlim())

    # (d) slip-rate functions at a few informative nodes
    final = ev["slip"][:, -1]
    pick = [int(np.argmin(np.abs(ev["node"] - (ev["nuc_node"] - 1))))]
    pick += [int(i) for i in np.argsort(final)[::-1][:3] if int(i) not in pick]
    greys = ["0.05", "0.35", "0.55", "0.72"]
    for c, i in zip(greys, pick):
        lab = "nucleation node" if i == pick[0] else f"x = {x_km[i]:.1f} km"
        ax[3].plot(t, ev["slip_rate"][i], color=c, lw=1.4, label=lab)
    ax[3].set_ylabel("slip rate (m/s)")
    ax[3].set_xlim(0, t[-1])
    ax[3].set_xlabel("time (s)")
    ax[3].set_title("(d) Slip-rate time functions", loc="left")
    ax[3].legend(loc="upper right", ncol=2, frameon=False)

    # (e) final slip along strike
    for f in faults:
        m = ev["fault"] == f
        order = np.argsort(x_km[m])
        ax[4].plot(x_km[m][order], final[m][order], lw=1.8,
                   color=FAULT_COLOURS[(int(f) - 1) % 10], label=f"ft{int(f)}")
    ax[4].set_ylabel("slip (m)")
    ax[4].set_xlabel("along-strike x (km)")
    ax[4].set_title("(e) Final slip", loc="left")
    if len(faults) > 1:
        ax[4].legend(loc="upper right", ncol=min(len(faults), 7), frameon=False)
    ax[4].set_xlim(ax[1].get_xlim())
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

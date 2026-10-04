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
N_TRACES = 120                                    # traces across the whole fault system
TRACE_GAIN = 2.0                                  # width of REF_SLIP_RATE, in trace spacings
REF_SLIP_RATE = 2.0                               # m/s; fixed across events for comparability


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
    ax[0].set_title(f"(a) {case.name}, event {eqid}: Mw {mw:.2f}, t = {ev['time_yr']:.1f} yr "
                    f"in the sequence, {ev['nout']} nodes, {ev['t_end']:.1f} s", loc="left")
    ax[0].set_xlim(0, t[-1])
    ax[0].set_xlabel("time (s)")

    faults = np.unique(ev["fault"])
    rt = ev["rupture_time"]
    ok = rt < 999
    # Whole fault system, not just the nodes that slipped: every node of the
    # table, with zero slip rate where the event did not reach.
    allx = case.nodes["x"] / 1e3
    xlim = (float(allx.min()), float(allx.max()))
    nall = len(allx)
    full_v = np.zeros((nall, len(t)), dtype=float)
    full_v[ev["node"]] = ev["slip_rate"]
    full_slip = np.zeros(nall)
    full_slip[ev["node"]] = ev["slip"][:, -1]
    all_fault = case.nodes["fault"]

    # (b) record section: each trace is one node's slip-rate function, drawn at
    # its along-strike position with amplitude to the right, time upward
    rs = ax[1]
    # Traces at a fixed spacing along each fault (same nodes for every event),
    # coloured by fault so overlapping strands stay distinguishable, plus the
    # nucleation node so a small rupture is never missed between traces.
    step = max(1, nall // N_TRACES)
    sel = []
    for f in np.unique(all_fault):
        idx = np.where(all_fault == f)[0]
        sel += list(idx[np.argsort(allx[idx])][::step])
    nuc = int(ev["nuc_node"]) - 1
    if 0 <= nuc < nall and nuc not in sel:
        sel.append(nuc)
    spacing = (xlim[1] - xlim[0]) / max(nall // step, 1)
    # Fixed amplitude scale, identical for every event of a case, so record
    # sections are comparable: REF_SLIP_RATE spans TRACE_GAIN trace spacings.
    gain = TRACE_GAIN * spacing / REF_SLIP_RATE
    for i in sel:
        col = FAULT_COLOURS[(int(all_fault[i]) - 1) % 10]
        if full_v[i].max() > 0:
            tr = allx[i] + full_v[i] * gain
            rs.fill_betweenx(t, allx[i], tr, color=col, alpha=0.55, lw=0)
            rs.plot(tr, t, color=col, lw=0.5)
        else:
            rs.plot([allx[i], allx[i]], [t[0], t[-1]], color="0.88", lw=0.3, zorder=0)
    rs.plot(x_km[ok], rt[ok], ".", ms=1.2, color="C3", label="rupture front")
    rs.plot(ev["nuc_x"] / 1e3, 0, marker="*", ms=11, color="C3", ls="none",
            label="nucleation")
    for f in np.unique(all_fault):
        rs.plot([], [], color=FAULT_COLOURS[(int(f) - 1) % 10], lw=3, label=f"ft{int(f)}")
    rs.plot([], [], color="0.05", lw=1, label=f"1 m/s = {gain:.2f} km")
    rs.set_ylim(0, t[-1])
    rs.set_ylabel("time (s)")
    rs.set_xlabel("along-strike x (km)")
    rs.set_title(f"(b) Record section: every {step}th node per fault, ~{spacing:.1f} km apart (grey = no slip)", loc="left")
    rs.legend(loc="upper right", frameon=True, ncol=2, fontsize=7.5)
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
    all_faults = np.unique(all_fault)
    for f in all_faults:
        m = all_fault == f
        order = np.argsort(allx[m])
        xs, ss = allx[m][order], full_slip[m][order]
        ax[3].plot(xs, np.zeros_like(xs), color="0.8", lw=0.8, zorder=0)
        ax[3].plot(xs, np.where(ss > 0, ss, np.nan), lw=1.8,
                   color=FAULT_COLOURS[(int(f) - 1) % 10], label=f"ft{int(f)}")
    ax[3].set_ylabel("slip (m)")
    ax[3].set_xlabel("along-strike x (km)")
    ax[3].set_title("(d) Final slip", loc="left")
    if len(all_faults) > 1:
        ax[3].legend(loc="upper right", ncol=min(len(all_faults), 7), frameon=False)
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

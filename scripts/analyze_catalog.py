#!/usr/bin/env python3
"""Analyze an aPlots/catalog.csv produced by `CATALOG=1 plotRuptureDynamics`.

Computes:
- Gutenberg-Richter b-value (Aki 1965 max-likelihood) above a chosen Mc.
- Magnitude-frequency distribution (cumulative + linear fit).
- Cycle-id vs magnitude scatter (warm-up / drift visible).
- Nucleation along strike vs cycle (color = magnitude).

Saves a 3-panel figure to <case>/aPlots/catalog_analysis.png.

Usage:
    python3 scripts/analyze_catalog.py [case_dir] [--mc 6.0]
"""
from __future__ import annotations

import argparse
import math
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def load_catalog(path: str) -> np.ndarray:
    return np.genfromtxt(path, delimiter=",", names=True, comments="#",
                         dtype=None, encoding=None)


def fit_b_lsq(bin_centers: np.ndarray, log10_cum: np.ndarray,
              mc: float, m_max: float) -> tuple[float, float, int]:
    """Least-squares slope of log10(N≥M) vs M in the [mc, m_max] window.
    Robust to characteristic-fault tails — matches what people visually fit."""
    sel = (bin_centers >= mc) & (bin_centers <= m_max) & np.isfinite(log10_cum)
    if sel.sum() < 3:
        return float("nan"), float("nan"), int(sel.sum())
    x, y = bin_centers[sel], log10_cum[sel]
    slope, intercept = np.polyfit(x, y, 1)
    # residual std → slope std
    resid = y - (slope * x + intercept)
    sx2   = np.sum((x - x.mean()) ** 2)
    s_b   = float(np.sqrt(np.sum(resid ** 2) / max(1, len(x) - 2) / sx2))
    return -float(slope), s_b, int(sel.sum())


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("case_dir", nargs="?", default=".")
    ap.add_argument("--mc", type=float, default=None,
                    help="Magnitude of completeness (default: auto = max bin in MFD)")
    ap.add_argument("--mmax", type=float, default=None,
                    help="Upper magnitude for LSQ fit window (default: max M with N>=5). "
                         "For characteristic-fault MFDs, set below the characteristic bump "
                         "(e.g. 7.0 for SAF) so the small-event b-value isn't blended in.")
    args = ap.parse_args()

    case_dir = os.path.abspath(args.case_dir)
    cat_path = os.path.join(case_dir, "aPlots", "catalog.csv")
    if not os.path.exists(cat_path):
        sys.exit(f"ERROR: {cat_path} not found. Run CATALOG=1 plotRuptureDynamics first.")

    cat = load_catalog(cat_path)
    eq      = cat["eqId"].astype(int)
    mag     = cat["magnitude"]
    nuc_x   = cat["nuc_x_km"]
    rup_dur = cat["rup_dur_s"]
    valid   = mag > 0                 # zero-magnitude = no rupture
    eq, mag, nuc_x, rup_dur = eq[valid], mag[valid], nuc_x[valid], rup_dur[valid]

    # MFD: cumulative count log10(N >= M) at uniform M bins
    bin_edges = np.arange(np.floor(mag.min() * 10) / 10.0,
                          np.ceil(mag.max() * 10) / 10.0 + 0.01, 0.1)
    cum = np.array([np.sum(mag >= m) for m in bin_edges], dtype=float)

    # Auto Mc = lowest non-empty bin; user can override
    mc = float(args.mc) if args.mc is not None else float(bin_edges[np.argmax(cum)])
    # Auto Mmax = highest bin where cum >= 5 events (statistical floor)
    if args.mmax is not None:
        mmax = float(args.mmax)
    else:
        idx = np.where(cum >= 5)[0]
        mmax = float(bin_edges[idx[-1]]) if idx.size else float(mag.max())

    log10_cum = np.where(cum > 0, np.log10(cum), np.nan)
    b_lsq, sb_lsq, n_lsq = fit_b_lsq(bin_edges, log10_cum, mc, mmax)

    print(f"events       : {len(mag)}")
    print(f"M range      : {mag.min():.2f} .. {mag.max():.2f}")
    print(f"Mc           : {mc:.2f}     Mmax (LSQ window) : {mmax:.2f}")
    print(f"b-value      : {b_lsq:.3f} ± {sb_lsq:.3f}   (LSQ on M=[{mc:.1f},{mmax:.1f}], {n_lsq} bins)")
    if len(mag) < 1000:
        print(f"NOTE         : only {len(mag)} events — b-value still noisy. Need ~1000+ for stable fit.")

    fig, ax = plt.subplots(1, 3, figsize=(18, 5))

    # Panel A: MFD with both fits
    ax[0].semilogy(bin_edges, cum, "o-", color="#1f77b4", label="cumulative N(M≥M)")
    if not math.isnan(b_lsq):
        # LSQ: line in [mc, mmax] window
        sel = (bin_edges >= mc) & (bin_edges <= mmax)
        x_l = bin_edges[sel]
        y_l = log10_cum[sel]
        intercept = np.nanmean(y_l + b_lsq * x_l)
        m_line = np.linspace(mc, mmax + 0.5, 50)
        ax[0].semilogy(m_line, 10 ** (intercept - b_lsq * m_line),
                       "--", color="red",
                       label=f"LSQ fit  b={b_lsq:.2f}±{sb_lsq:.2f}  (M∈[{mc:.1f},{mmax:.1f}])")
    ax[0].axvline(mc, color="grey", linestyle=":", label=f"Mc={mc:.2f}")
    ax[0].axvline(mmax, color="grey", linestyle="--", alpha=0.5, label=f"Mmax={mmax:.2f}")
    ax[0].set_xlabel("Magnitude")
    ax[0].set_ylabel("Cumulative count")
    ax[0].set_title("Magnitude-frequency distribution")
    ax[0].grid(True, which="both", alpha=0.3)
    ax[0].legend()

    # Panel B: cycle id vs magnitude
    ax[1].scatter(eq, mag, s=12, c=mag, cmap="viridis", alpha=0.8)
    ax[1].set_xlabel("event id (cycle)")
    ax[1].set_ylabel("magnitude")
    ax[1].set_title("Magnitude over cycles (warm-up / drift)")
    ax[1].grid(True, alpha=0.3)

    # Panel C: nucleation along strike vs cycle, color = magnitude
    sc = ax[2].scatter(eq, nuc_x, s=12 + 6 * (mag - mag.min()),
                       c=mag, cmap="viridis", alpha=0.8)
    fig.colorbar(sc, ax=ax[2], label="M")
    ax[2].set_xlabel("event id (cycle)")
    ax[2].set_ylabel("nucleation x (km)")
    ax[2].set_title("Nucleation along strike (size = M)")
    ax[2].grid(True, alpha=0.3)

    out = os.path.join(case_dir, "aPlots", "catalog_analysis.png")
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()

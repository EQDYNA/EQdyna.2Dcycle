#!/usr/bin/env python3
"""Resample gulang fault .gmt control points to >= 3 km spacing.

Reads ft{i}.gmt.txt.orig (creating it from current .gmt.txt on first run),
fits a per-fault smoothing spline x(s), y(s) parameterised by arc length,
and writes back ft{i}.gmt.txt with control points at uniform arc-length
step >= 3 km. Endpoints are preserved exactly.

Run this once before meshgen.py to replace the dense (~30 m median)
original digitisation with a SAF-like coarse control set.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
from scipy.interpolate import UnivariateSpline

HERE = Path(__file__).parent
INPUT_DIR = HERE / "user_fault_geometry_input"

MIN_SPACING_KM = 3.0
SMOOTH_SIGMA_KM = 0.010  # ~10 m noise floor


def preprocess(name: str) -> tuple[int, int, float]:
    raw = INPUT_DIR / f"{name}.gmt.txt"
    backup = INPUT_DIR / f"{name}.gmt.txt.orig"

    if not backup.exists():
        backup.write_bytes(raw.read_bytes())

    data = np.loadtxt(backup)
    x, y = data[:, 0], data[:, 1]
    n_orig = len(x)

    s = np.concatenate([[0.0], np.cumsum(np.hypot(np.diff(x), np.diff(y)))])
    L = s[-1]

    smooth = n_orig * SMOOTH_SIGMA_KM ** 2
    spl_x = UnivariateSpline(s, x, k=3, s=smooth)
    spl_y = UnivariateSpline(s, y, k=3, s=smooth)

    n_seg = max(1, int(np.floor(L / MIN_SPACING_KM)))
    s_new = np.linspace(0.0, L, n_seg + 1)
    x_new = spl_x(s_new)
    y_new = spl_y(s_new)
    x_new[0], y_new[0] = x[0], y[0]
    x_new[-1], y_new[-1] = x[-1], y[-1]

    with open(raw, "w") as f:
        for xi, yi in zip(x_new, y_new):
            f.write(f"{xi:.6f} {yi:.6f}\n")

    step = L / n_seg
    return n_orig, len(x_new), step


def main() -> None:
    print(f"Preprocessing gulang control points to >= {MIN_SPACING_KM} km spacing")
    print(f"Source : {INPUT_DIR}/ft{{i}}.gmt.txt.orig (auto-backed up on first run)")
    print(f"Output : {INPUT_DIR}/ft{{i}}.gmt.txt (overwritten)")
    print()
    print(f"{'fault':6s}  {'n_orig':>7s}  {'n_new':>5s}  {'step_km':>7s}")
    for i in range(1, 6):
        n_orig, n_new, step = preprocess(f"ft{i}")
        print(f"ft{i:<5d}  {n_orig:7d}  {n_new:5d}  {step:7.3f}")


if __name__ == "__main__":
    main()

#! /usr/bin/env python3
"""Compare the 5 totalop variables over strike for a chosen cycle id,
across an arbitrary list of case dirs plus the Pangaea reference.

Usage:
  compare_cycle_over_strike.py --cycle 1 --model A \
      --case work/paper.saf.A.omp1 --case work/paper.saf.A.omp2 --case work/paper.saf.A.omp4

  # skip the Pangaea overlay
  compare_cycle_over_strike.py --cycle 1 --case dir1 --case dir2 --no-reference

Output: one PNG with 5 rows (shear, normal, slip, slip_rate, rupture_time) and
every case plotted as a separate line along the fault strike.
"""

from __future__ import annotations

import argparse
import os
import re
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


_FIX_EXP_RE = re.compile(r'(\d)([+-]\d{3,})')  # insert 'E' when Fortran dropped it at 3-digit exp


def _robust_loadtxt(path: str) -> np.ndarray:
    """np.loadtxt but tolerates Fortran scientific notation missing 'E'
    (e.g. '0.8384675-101' → '0.8384675E-101')."""
    rows = []
    with open(path) as f:
        for line in f:
            s = _FIX_EXP_RE.sub(r'\1E\2', line)
            parts = s.split()
            if not parts:
                continue
            rows.append([float(x) for x in parts])
    return np.asarray(rows, dtype=np.float64)

ROOT = os.environ.get('EQDYNA2DCYCLEROOT') or os.path.abspath(
    os.path.join(os.path.dirname(__file__), os.pardir))

PANGAEA = {
    'A': os.path.join(ROOT, 'archive/published/pangaea_results/work_vis7_fs0.5'),
    'B': os.path.join(ROOT, 'archive/published/pangaea_results/work_vis4.2_fs0.3'),
    'C': os.path.join(ROOT, 'archive/published/pangaea_results/work_vis12_fs0.7'),
}

COL_NAMES = ['shear (MPa)', 'normal (MPa)', 'slip (m)', 'slip_rate (m/s)', 'rupture_time (s)']
COL_SCALE = [1e-6, 1e-6, 1.0, 1.0, 1.0]


def read_segmented_totalop(case_dir: str) -> dict[int, np.ndarray]:
    """Return {icstart: raw (nrows, 5)} from totalop.txt<ic> in case_dir or its aRawSimuData/."""
    search_dirs = [case_dir, os.path.join(case_dir, 'aRawSimuData')]
    segments = {}
    for d in search_dirs:
        if not os.path.isdir(d):
            continue
        for fname in os.listdir(d):
            if not fname.startswith('totalop.txt'):
                continue
            try:
                icstart = int(fname[len('totalop.txt'):])
            except ValueError:
                continue
            segments[icstart] = os.path.join(d, fname)
    if not segments:
        raise FileNotFoundError(f'no totalop.txt<ic> files under {case_dir}')
    return {ic: _robust_loadtxt(p) for ic, p in sorted(segments.items())}


def extract_cycle(segments: dict[int, np.ndarray], cycle: int, totftnode: int) -> np.ndarray | None:
    starts = sorted(segments.keys())
    for i, ic in enumerate(starts):
        arr = segments[ic]
        n = arr.shape[0] // totftnode
        nxt = starts[i + 1] if i + 1 < len(starts) else ic + n
        if ic <= cycle < nxt and cycle < ic + n:
            k = cycle - ic
            return arr[k * totftnode:(k + 1) * totftnode, :]
    return None


def infer_totftnode(ref_segs: dict | None, fallback: int = 2242) -> int:
    if ref_segs and len(ref_segs) >= 2:
        starts = sorted(ref_segs.keys())
        first_len = starts[1] - starts[0]
        return ref_segs[starts[0]].shape[0] // first_len
    return fallback


def load_strike(case_dir: str, totftnode: int) -> np.ndarray:
    for candidate in ('nsmpTanLen.txt', 'nsmpnv.txt'):
        p = os.path.join(case_dir, candidate)
        if os.path.isfile(p):
            raw = np.loadtxt(p)
            if raw.shape[0] >= totftnode:
                # col3 = segment length in m → cumulative km
                return np.cumsum(raw[:totftnode, 2]) / 1e3
    return np.arange(totftnode, dtype=float)


def main(argv):
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument('--case', action='append', required=True,
                   help='case dir; repeat to overlay multiple')
    p.add_argument('--cycle', type=int, required=True, help='absolute cycle id')
    p.add_argument('--model', choices=['A', 'B', 'C'], default='A',
                   help='Pangaea reference model (default A)')
    p.add_argument('--no-reference', action='store_true', help='skip Pangaea overlay')
    p.add_argument('--out', default=None, help='output png path')
    p.add_argument('--label', action='append', default=None,
                   help='custom label per case (repeat); defaults to basename')
    args = p.parse_args(argv)

    # resolve cases
    case_dirs = [os.path.abspath(c) for c in args.case]
    labels = args.label if args.label else [os.path.basename(c) for c in case_dirs]
    if len(labels) != len(case_dirs):
        print('ERROR: number of --label entries must match --case entries', file=sys.stderr)
        return 2

    # reference
    ref_segs = None
    if not args.no_reference:
        try:
            ref_segs = read_segmented_totalop(PANGAEA[args.model])
        except FileNotFoundError as e:
            print(f'WARN: {e}; plotting without reference', file=sys.stderr)

    nft = infer_totftnode(ref_segs)
    print(f'totftnode = {nft}')

    # pull cycle from each case
    case_arrs = []
    for cdir, lab in zip(case_dirs, labels):
        try:
            segs = read_segmented_totalop(cdir)
            arr = extract_cycle(segs, args.cycle, nft)
            if arr is None:
                print(f'  {lab}: cycle {args.cycle} NOT FOUND (skipping)')
                continue
            case_arrs.append((lab, arr))
            print(f'  {lab}: cycle {args.cycle} OK ({arr.shape[0]} rows)')
        except FileNotFoundError as e:
            print(f'  {lab}: {e} (skipping)')

    ref_arr = None if ref_segs is None else extract_cycle(ref_segs, args.cycle, nft)
    if ref_arr is None and ref_segs is not None:
        print(f'WARN: Pangaea cycle {args.cycle} not found; no reference overlay')

    if not case_arrs and ref_arr is None:
        print('ERROR: nothing to plot', file=sys.stderr)
        return 1

    # pick a strike axis from the first existing case (all should share mesh)
    strike_src = case_dirs[0] if case_dirs else PANGAEA[args.model]
    strike = load_strike(strike_src, nft)

    fig, axes = plt.subplots(5, 1, figsize=(12, 14), sharex=True)
    colors = plt.get_cmap('tab10').colors
    for col, (ax, name, scale) in enumerate(zip(axes, COL_NAMES, COL_SCALE)):
        for i, (lab, arr) in enumerate(case_arrs):
            ax.plot(strike, arr[:, col] * scale, lw=1.2,
                    color=colors[i % len(colors)], label=lab)
        if ref_arr is not None:
            ax.plot(strike, ref_arr[:, col] * scale, lw=1.4,
                    color='k', ls='--', alpha=0.8,
                    label=f'Pangaea Model {args.model}')
        ax.set_ylabel(name)
        ax.grid(alpha=0.3)
        if col == 0:
            ax.legend(loc='best', fontsize=8, ncol=2)
    axes[-1].set_xlabel('along-strike cumulative distance (km)' if strike[-1] > strike[0] + 1
                        else 'fault node index')
    fig.suptitle(f'cycle {args.cycle} — model {args.model} comparison across cases')
    fig.tight_layout()

    out = args.out or f'compare_cycle_{args.cycle}_model{args.model}.png'
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f'wrote {out}')
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))

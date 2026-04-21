#!/usr/bin/env python3
"""Plot near-fault mesh for one or more cases for visual comparison.

Usage: plotMeshNearFault.py <case_dir> [<case_dir> ...] [--out OUT] [--zoom R_KM]

Produces one PNG per case showing a zoom window around each fault's
center with all quad cells drawn, fault master nodes highlighted.
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection


def read_Cmesh(root):
    with open(os.path.join(root, 'FE_Global.txt')) as f:
        return int(f.readline().split('#')[0].strip())


def read_nfn(root):
    with open(os.path.join(root, 'meshGeneralInfo.txt')) as f:
        next(f)
        return list(map(int, f.readline().split()))


def load_case(root):
    cmesh = read_Cmesh(root)
    to_km = 1.0e-3 if cmesh == 2 else 1.0   # coords to km
    v = np.loadtxt(os.path.join(root, 'vert.txt')) * to_km
    f = np.loadtxt(os.path.join(root, 'fac.txt'), dtype=int)
    if f.min() == 1:
        f = f - 1
    nsmp = np.loadtxt(os.path.join(root, 'nsmp.txt'), dtype=int)
    if nsmp.min() == 1:
        nsmp = nsmp - 1
    nfn = read_nfn(root)
    return v, f, nsmp, nfn, cmesh


def plot_one(root, zoom_km, out_path):
    v, f, nsmp, nfn, cmesh = load_case(root)
    maxftnode = max(nfn)
    ncells = f.shape[0]
    # Build vertex coords for each cell (ncells, 4, 2)
    quads = v[f]

    nFt = len(nfn)
    fig, axes = plt.subplots(1, nFt, figsize=(6 * nFt, 6), squeeze=False)
    for iFt, n in enumerate(nfn):
        master = nsmp[iFt * maxftnode: iFt * maxftnode + n, 0]
        slave = nsmp[iFt * maxftnode: iFt * maxftnode + n, 1]
        mc = v[master]
        sc = v[slave]
        # zoom center: fault midpoint
        mid = mc[n // 2]
        xlim = (mid[0] - zoom_km, mid[0] + zoom_km)
        ylim = (mid[1] - zoom_km, mid[1] + zoom_km)

        # Only cells with all 4 verts inside the zoom box (for speed)
        in_box = (
            (quads[:, :, 0] >= xlim[0]).all(axis=1) &
            (quads[:, :, 0] <= xlim[1]).all(axis=1) &
            (quads[:, :, 1] >= ylim[0]).all(axis=1) &
            (quads[:, :, 1] <= ylim[1]).all(axis=1)
        )
        sel = quads[in_box]

        ax = axes[0, iFt]
        pc = PolyCollection(sel, facecolors='none', edgecolors='0.4',
                            linewidths=0.4)
        ax.add_collection(pc)
        # fault trace
        ax.plot(mc[:, 0], mc[:, 1], 'r-', linewidth=1.5, label='fault master')
        ax.plot(mc[:, 0], mc[:, 1], 'r.', markersize=2.5)
        # slave (usually coincident, offset for visual?)
        ax.plot(sc[:, 0], sc[:, 1], 'b.', markersize=1.5, alpha=0.5,
                label='fault slave')

        ax.set_xlim(xlim); ax.set_ylim(ylim)
        ax.set_aspect('equal')
        ax.set_xlabel('x (km)'); ax.set_ylabel('y (km)')
        ax.set_title(f'Fault {iFt + 1} ({n} nodes); {sel.shape[0]} cells shown')
        ax.grid(True, alpha=0.3)
        if iFt == 0:
            ax.legend(fontsize=8, loc='best')

    label = os.path.basename(root.rstrip('/'))
    fig.suptitle(f'{label}  (C_mesh={cmesh})  zoom ±{zoom_km} km',
                 fontsize=13)
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    plt.close()
    print(f'Wrote {out_path}')


def main():
    zoom = 10.0
    out_dir = '.'
    args = []
    it = iter(sys.argv[1:])
    for a in it:
        if a == '--zoom':
            zoom = float(next(it))
        elif a == '--out':
            out_dir = next(it)
        else:
            args.append(a)
    if not args:
        print(__doc__); sys.exit(1)
    os.makedirs(out_dir, exist_ok=True)
    for case in args:
        name = os.path.basename(case.rstrip('/'))
        out = os.path.join(out_dir, f'meshNearFault_{name}.png')
        plot_one(case, zoom, out)


if __name__ == '__main__':
    main()

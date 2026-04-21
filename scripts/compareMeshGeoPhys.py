#!/usr/bin/env python3
"""Compare fault geometry, tangent/normal vectors, and loading between two cases.

Usage: compareMeshGeoPhys.py <case_A_dir> <case_B_dir> [out_dir]

Each case dir must contain:
  - vert.txt, nsmp.txt, meshGeneralInfo.txt
  - Either Rate_direction.txt (C_mesh=2) or nsmpGeoPhys.txt (C_mesh=3)
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def read_FE_Global_Cmesh(root):
    with open(os.path.join(root, 'FE_Global.txt')) as f:
        return int(f.readline().split('#')[0].strip())


def read_nfnodes(root):
    with open(os.path.join(root, 'meshGeneralInfo.txt')) as f:
        next(f)
        return list(map(int, f.readline().split()))


def load_case(root):
    cmesh = read_FE_Global_Cmesh(root)
    to_m = 1.0 if cmesh == 2 else 1.0e3
    v = np.loadtxt(os.path.join(root, 'vert.txt')) * to_m
    nsmp = np.loadtxt(os.path.join(root, 'nsmp.txt'), dtype=int)
    if nsmp.min() == 1:
        nsmp = nsmp - 1
    nfn = read_nfnodes(root)
    maxftnode = max(nfn)

    faults = []
    for i, n in enumerate(nfn):
        master = nsmp[i * maxftnode: i * maxftnode + n, 0]
        faults.append({'coords': v[master], 'n': n})

    # Loading rate + angle
    if cmesh == 2:
        rd = np.loadtxt(os.path.join(root, 'Rate_direction.txt'))
        for i, n in enumerate(nfn):
            blk = rd[i * maxftnode: i * maxftnode + n]
            faults[i]['rate'] = blk[:, 0]
            faults[i]['angle'] = blk[:, 1]
            # tangent from chord of coords (for plotting vectors)
            faults[i]['tx'], faults[i]['ty'] = _chord_tangent(faults[i]['coords'])
    else:
        gp = np.loadtxt(os.path.join(root, 'nsmpGeoPhys.txt'))
        for i, n in enumerate(nfn):
            blk = gp[i * maxftnode: i * maxftnode + n]
            faults[i]['tx'] = blk[:, 0]
            faults[i]['ty'] = blk[:, 1]
            faults[i]['rate'] = blk[:, 5]
            faults[i]['angle'] = blk[:, 6]
    return {'cmesh': cmesh, 'faults': faults, 'nfn': nfn}


def _chord_tangent(coords):
    tx = np.zeros(len(coords))
    ty = np.zeros(len(coords))
    for i in range(1, len(coords) - 1):
        d = coords[i + 1] - coords[i - 1]
        L = np.linalg.norm(d)
        tx[i], ty[i] = d[0] / L, d[1] / L
    for i in (0, -1):
        j = 1 if i == 0 else -2
        d = coords[i] - coords[j] if i == -1 else coords[j] - coords[i]
        L = np.linalg.norm(d)
        tx[i], ty[i] = d[0] / L, d[1] / L
    return tx, ty


def plot_fault_traces(caseA, labelA, caseB, labelB, out_path, n_vec=20):
    fig, ax = plt.subplots(1, 1, figsize=(12, 10))
    colors_A = ['C0', 'C1', 'C2', 'C3']
    colors_B = ['C4', 'C5', 'C6', 'C7']
    for i, ft in enumerate(caseA['faults']):
        c = ft['coords'] / 1e3  # m -> km
        ax.plot(c[:, 0], c[:, 1], '-', color=colors_A[i % 4],
                label=f'{labelA} Ft{i + 1}', linewidth=1.2)
        # unit tangent vectors (downsampled)
        step = max(1, len(c) // n_vec)
        idx = np.arange(0, len(c), step)
        ax.quiver(c[idx, 0], c[idx, 1], ft['tx'][idx], ft['ty'][idx],
                  color=colors_A[i % 4], scale=30, width=0.003, alpha=0.7)
    for i, ft in enumerate(caseB['faults']):
        c = ft['coords'] / 1e3
        ax.plot(c[:, 0], c[:, 1], '--', color=colors_B[i % 4],
                label=f'{labelB} Ft{i + 1}', linewidth=1.2)
        step = max(1, len(c) // n_vec)
        idx = np.arange(0, len(c), step)
        ax.quiver(c[idx, 0], c[idx, 1], ft['tx'][idx], ft['ty'][idx],
                  color=colors_B[i % 4], scale=30, width=0.003, alpha=0.7)
    ax.set_xlabel('x (km)')
    ax.set_ylabel('y (km)')
    ax.set_title('Fault traces and unit tangent vectors')
    ax.set_aspect('equal')
    ax.legend(loc='best', fontsize=10)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    plt.close()
    print(f'Wrote {out_path}')


def plot_rate_angle_vs_node(caseA, labelA, caseB, labelB, out_path):
    nft = max(len(caseA['faults']), len(caseB['faults']))
    fig, axes = plt.subplots(nft, 2, figsize=(14, 4 * nft), squeeze=False)
    for i in range(nft):
        ax_r = axes[i, 0]
        ax_a = axes[i, 1]
        if i < len(caseA['faults']):
            ft = caseA['faults'][i]
            ax_r.plot(ft['rate'], 'C0-', label=labelA, linewidth=1.5)
            ax_a.plot(ft['angle'], 'C0-', label=labelA, linewidth=1.5)
        if i < len(caseB['faults']):
            ft = caseB['faults'][i]
            ax_r.plot(ft['rate'], 'C1--', label=labelB, linewidth=1.5)
            ax_a.plot(ft['angle'], 'C1--', label=labelB, linewidth=1.5)
        ax_r.set_title(f'Fault {i + 1}: loading rate')
        ax_r.set_xlabel('along-strike node ID')
        ax_r.set_ylabel('rate (s$^{-1}$)')
        ax_r.legend(fontsize=9)
        ax_r.grid(True, alpha=0.3)
        ax_a.set_title(f'Fault {i + 1}: angle (loading vs fault tangent)')
        ax_a.set_xlabel('along-strike node ID')
        ax_a.set_ylabel('angle (deg)')
        ax_a.legend(fontsize=9)
        ax_a.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    plt.close()
    print(f'Wrote {out_path}')


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(1)
    rootA = sys.argv[1]
    rootB = sys.argv[2]
    out_dir = sys.argv[3] if len(sys.argv) > 3 else '.'
    os.makedirs(out_dir, exist_ok=True)

    labelA = os.path.basename(rootA.rstrip('/'))
    labelB = os.path.basename(rootB.rstrip('/'))
    print(f'Loading {labelA} ({rootA}) ...')
    A = load_case(rootA)
    print(f'  C_mesh={A["cmesh"]}, nfnodes={A["nfn"]}')
    print(f'Loading {labelB} ({rootB}) ...')
    B = load_case(rootB)
    print(f'  C_mesh={B["cmesh"]}, nfnodes={B["nfn"]}')

    plot_fault_traces(A, labelA, B, labelB,
                      os.path.join(out_dir, 'faultTraces_tangents.png'))
    plot_rate_angle_vs_node(A, labelA, B, labelB,
                            os.path.join(out_dir, 'loadingRate_angle_vs_node.png'))


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""Mesh quality checks for EQdyna.2Dcycle cases.

Usage: checkMeshQuality.py <case_dir> [<case_dir> ...]

Runs a suite of checks on each case's mesh files:
  - Split-node corner cells (master-only / slave-only / mixed)
  - Mixed master+slave cells (topology bug)
  - Fault-node coverage (each fault node appears in >=1 cell)
  - Quad cell geometry (area, aspect ratio, min interior angle)
  - Orphan / degenerate cells
"""
import os
import sys
import numpy as np


def _load_case(root):
    v = np.loadtxt(os.path.join(root, 'vert.txt'))
    f = np.loadtxt(os.path.join(root, 'fac.txt'), dtype=int)
    if f.min() == 1:
        f = f - 1
    nsmp = np.loadtxt(os.path.join(root, 'nsmp.txt'), dtype=int)
    if nsmp.min() == 1:
        nsmp = nsmp - 1
    with open(os.path.join(root, 'meshGeneralInfo.txt')) as fh:
        next(fh)
        nfn = list(map(int, fh.readline().split()))
    return v, f, nsmp, nfn


def _masters_slaves(nsmp, nfn):
    maxftnode = max(nfn)
    out = []
    for i, n in enumerate(nfn):
        blk = nsmp[i * maxftnode: i * maxftnode + n]
        out.append((set(blk[:, 0].tolist()), set(blk[:, 1].tolist())))
    return out


def checkSplitNodes(v, f, nsmp, nfn):
    """Corner-cell split-node classification per fault.

    Returns list of dicts per fault with counts of master-only, slave-only,
    mixed, and cells touching fault but with no ft node.
    """
    ms = _masters_slaves(nsmp, nfn)
    results = []
    for iFt, (master, slave) in enumerate(ms):
        m_edge = s_edge = 0
        m_corner = s_corner = 0
        mixed = 0
        for cell in f:
            nm = sum(1 for x in cell if x in master)
            ns = sum(1 for x in cell if x in slave)
            if nm > 0 and ns > 0:
                mixed += 1
            elif nm == 1 and ns == 0:
                m_corner += 1
            elif ns == 1 and nm == 0:
                s_corner += 1
            elif nm >= 2 and ns == 0:
                m_edge += 1
            elif ns >= 2 and nm == 0:
                s_edge += 1
        results.append({
            'fault': iFt + 1, 'nfnodes': nfn[iFt],
            'm_edge': m_edge, 's_edge': s_edge,
            'm_corner': m_corner, 's_corner': s_corner,
            'mixed': mixed,
        })
    return results


def checkFaultNodeCoverage(f, nsmp, nfn):
    """Each fault node (master and slave) should appear in >=1 cell."""
    all_cell_nodes = set(f.flatten().tolist())
    ms = _masters_slaves(nsmp, nfn)
    results = []
    for iFt, (master, slave) in enumerate(ms):
        m_orphan = master - all_cell_nodes
        s_orphan = slave - all_cell_nodes
        results.append({
            'fault': iFt + 1,
            'master_orphans': len(m_orphan), 'slave_orphans': len(s_orphan),
        })
    return results


def _quad_geom(pts):
    """Return (area, aspect_ratio, min_angle_deg) for one quad (4x2)."""
    p = pts
    # Shoelace area for quad
    x, y = p[:, 0], p[:, 1]
    area = 0.5 * abs(x[0] * (y[1] - y[3]) + x[1] * (y[2] - y[0]) +
                     x[2] * (y[3] - y[1]) + x[3] * (y[0] - y[2]))
    # Edge lengths
    e = np.linalg.norm(np.diff(np.vstack([p, p[:1]]), axis=0), axis=1)
    aspect = e.max() / max(e.min(), 1e-30)
    # Interior angles
    angs = []
    for i in range(4):
        a = p[(i - 1) % 4] - p[i]
        b = p[(i + 1) % 4] - p[i]
        cos = np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b) + 1e-30)
        angs.append(np.degrees(np.arccos(np.clip(cos, -1.0, 1.0))))
    return area, aspect, min(angs)


def checkCellGeometry(v, f):
    """Return summary stats + counts of bad cells."""
    areas = np.zeros(len(f))
    aspects = np.zeros(len(f))
    min_angs = np.zeros(len(f))
    degenerate = 0
    for i, cell in enumerate(f):
        pts = v[cell]
        a, r, mn = _quad_geom(pts)
        areas[i], aspects[i], min_angs[i] = a, r, mn
        if a <= 0.0:
            degenerate += 1
    return {
        'ncells': len(f),
        'area_min': areas.min(), 'area_median': np.median(areas),
        'area_max': areas.max(),
        'aspect_median': np.median(aspects), 'aspect_p99': np.percentile(aspects, 99),
        'aspect_max': aspects.max(),
        'min_angle_min': min_angs.min(), 'min_angle_median': np.median(min_angs),
        'bad_angle_lt20': int((min_angs < 20).sum()),
        'bad_aspect_gt10': int((aspects > 10).sum()),
        'degenerate': degenerate,
    }


def report(root):
    print(f'\n=== {root} ===')
    v, f, nsmp, nfn = _load_case(root)
    print(f'  {len(v)} nodes, {len(f)} cells, {len(nfn)} fault(s), nfnodes={nfn}')

    print('\n  [split-node classification per fault]')
    for r in checkSplitNodes(v, f, nsmp, nfn):
        flag = '' if r['mixed'] == 0 else '  **MIXED CELLS — BUG**'
        print(f"    Ft{r['fault']} (n={r['nfnodes']}): "
              f"edge m/s={r['m_edge']}/{r['s_edge']}  "
              f"corner m/s={r['m_corner']}/{r['s_corner']}  "
              f"mixed={r['mixed']}{flag}")

    print('\n  [fault-node coverage]')
    for r in checkFaultNodeCoverage(f, nsmp, nfn):
        flag = '' if (r['master_orphans'] + r['slave_orphans']) == 0 else '  **ORPHAN**'
        print(f"    Ft{r['fault']}: master_orphans={r['master_orphans']}  "
              f"slave_orphans={r['slave_orphans']}{flag}")

    print('\n  [cell geometry]')
    g = checkCellGeometry(v, f)
    print(f"    area:   min={g['area_min']:.3e}  med={g['area_median']:.3e}  "
          f"max={g['area_max']:.3e}")
    print(f"    aspect: med={g['aspect_median']:.2f}  p99={g['aspect_p99']:.2f}  "
          f"max={g['aspect_max']:.2f}")
    print(f"    min interior angle: min={g['min_angle_min']:.1f}  "
          f"med={g['min_angle_median']:.1f}")
    print(f"    bad cells: angle<20°={g['bad_angle_lt20']}  "
          f"aspect>10={g['bad_aspect_gt10']}  degenerate={g['degenerate']}")


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    for case in sys.argv[1:]:
        report(case)


if __name__ == '__main__':
    main()

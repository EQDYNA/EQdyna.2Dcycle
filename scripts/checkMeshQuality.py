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

# Hard-failure thresholds. A mesh violating any of these must not drive a run.
MIN_ANGLE_DEG = 20.0     # sliver corner
MAX_ANGLE_DEG = 160.0    # near-flat corner; the old report never looked at this
MAX_ASPECT = 10.0
# Two faults closer than this many elements: gmsh cannot recombine between
# them, leaves triangles, and fac.txt drops them (issue #1).
MIN_FAULT_GAP_ELEMENTS = 3.0


def _find(root, name):
    """Locate a mesh file in the case dir or in fem_mesh_output/.

    meshgen.py writes to fem_mesh_output/ and only run.sh copies the files up
    to the case directory, so a freshly meshed case has them in exactly one
    place -- and it is not the case root. Reading the root unconditionally made
    this script fail with a bare FileNotFoundError on any case that had been
    meshed but not yet run, which is precisely when a mesh check is wanted.
    The .msh and fac.txt lookups further down already searched both.
    """
    for d in (root, os.path.join(root, 'fem_mesh_output')):
        cand = os.path.join(d, name)
        if os.path.exists(cand):
            return cand
    raise SystemExit(f"error: {name} not found in {root} or {root}/fem_mesh_output")


def _load_case(root):
    v = np.loadtxt(_find(root, 'vert.txt'))
    f = np.loadtxt(_find(root, 'fac.txt'), dtype=int)
    if f.min() == 1:
        f = f - 1
    nsmp = np.loadtxt(_find(root, 'nsmp.txt'), dtype=int)
    if nsmp.min() == 1:
        nsmp = nsmp - 1
    with open(_find(root, 'meshGeneralInfo.txt')) as fh:
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
    return area, aspect, min(angs), max(angs)


def checkCellGeometry(v, f):
    """Return summary stats + counts of bad cells."""
    areas = np.zeros(len(f))
    aspects = np.zeros(len(f))
    min_angs = np.zeros(len(f))
    max_angs = np.zeros(len(f))
    degenerate = 0
    for i, cell in enumerate(f):
        pts = v[cell]
        a, r, mn, mx = _quad_geom(pts)
        areas[i], aspects[i], min_angs[i], max_angs[i] = a, r, mn, mx
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
        # The MAXIMUM interior angle matters as much as the minimum: a cell
        # with a near-180 deg corner is a sliver even when its smallest angle
        # looks fine, and the old report never looked at it.
        'max_angle_max': float(max_angs.max()),
        'bad_angle_gt160': int((max_angs > MAX_ANGLE_DEG).sum()),
    }


def checkElementCount(root):
    """fac.txt must hold every 2D element the .msh does.

    meshGenLib writes only quads. Where gmsh cannot recombine a region it
    leaves triangles, and those are dropped without a word, leaving holes in
    the FE mesh -- which is how split nodes end up orphaned (issue #1). This
    is invisible to every other check here, because they all read fac.txt.
    """
    msh = None
    for d in (root, os.path.join(root, 'fem_mesh_output')):
        cand = os.path.join(d, 'eqdynaMesh.msh')
        if os.path.exists(cand):
            msh = cand
            break
    if msh is None:
        return None
    try:
        import meshio
    except ImportError:
        return None
    m = meshio.read(msh)
    counts = {}
    for cb in m.cells:
        if cb.type in ('quad', 'triangle'):
            counts[cb.type] = counts.get(cb.type, 0) + len(cb.data)
    fac = None
    for d in (root, os.path.join(root, 'fem_mesh_output')):
        cand = os.path.join(d, 'fac.txt')
        if os.path.exists(cand):
            fac = cand
            break
    n_fac = sum(1 for _ in open(fac)) if fac else 0
    return {'quad': counts.get('quad', 0), 'triangle': counts.get('triangle', 0),
            'fac_rows': n_fac,
            'dropped': counts.get('quad', 0) + counts.get('triangle', 0) - n_fac}


def report(root):
    """Print the report and return a list of hard failures."""
    failures = []
    print(f'\n=== {root} ===')
    v, f, nsmp, nfn = _load_case(root)
    print(f'  {len(v)} nodes, {len(f)} cells, {len(nfn)} fault(s), nfnodes={nfn}')

    print('\n  [split-node classification per fault]')
    for r in checkSplitNodes(v, f, nsmp, nfn):
        flag = '' if r['mixed'] == 0 else '  **MIXED CELLS — BUG**'
        if r['mixed']:
            failures.append(f"Ft{r['fault']}: {r['mixed']} mixed cells")
        print(f"    Ft{r['fault']} (n={r['nfnodes']}): "
              f"edge m/s={r['m_edge']}/{r['s_edge']}  "
              f"corner m/s={r['m_corner']}/{r['s_corner']}  "
              f"mixed={r['mixed']}{flag}")

    print('\n  [fault-node coverage]')
    for r in checkFaultNodeCoverage(f, nsmp, nfn):
        flag = '' if (r['master_orphans'] + r['slave_orphans']) == 0 else '  **ORPHAN**'
        if r['master_orphans'] or r['slave_orphans']:
            failures.append(
                f"Ft{r['fault']}: {r['master_orphans']} master / "
                f"{r['slave_orphans']} slave orphan(s) - the split node has no "
                f"element to transmit traction")
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
    print(f"    max interior angle: {g['max_angle_max']:.1f}")
    print(f"    bad cells: angle<20°={g['bad_angle_lt20']}  angle>160°={g['bad_angle_gt160']}  "
          f"aspect>10={g['bad_aspect_gt10']}  degenerate={g['degenerate']}")
    for key, label in (('degenerate', 'degenerate cells'),
                       ('bad_angle_lt20', 'cells with an interior angle < 20 deg'),
                       ('bad_angle_gt160', f'cells with an interior angle > {MAX_ANGLE_DEG:.0f} deg'),
                       ('bad_aspect_gt10', 'cells with aspect ratio > 10')):
        if g[key]:
            failures.append(f"{g[key]} {label}")

    ec = checkElementCount(root)
    print('\n  [element export]')
    if ec is None:
        print('    no .msh or meshio unavailable, cannot verify export')
    else:
        print(f"    .msh: {ec['quad']} quads + {ec['triangle']} triangles;  "
              f"fac.txt: {ec['fac_rows']} rows")
        if ec['dropped']:
            print(f"    **{ec['dropped']} ELEMENT(S) DROPPED FROM fac.txt**")
            failures.append(
                f"{ec['dropped']} element(s) present in the .msh but missing from "
                f"fac.txt - the FE mesh has holes")
        if ec['triangle']:
            failures.append(
                f"{ec['triangle']} triangle(s) in the mesh; gmsh could not "
                f"recombine there, usually because two faults pass closer than "
                f"one element - refine dxy in that region")

    return failures


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    allFailures = []
    for case in sys.argv[1:]:
        for f in report(case):
            allFailures.append(f'{case}: {f}')
    print()
    if allFailures:
        print(f'FAILED: {len(allFailures)} mesh problem(s)')
        for f in allFailures:
            print(f'  - {f}')
        print('\nThese are hard failures. A mesh with orphaned split nodes, '
              'dropped\nelements or degenerate cells must not be used for a run.')
        sys.exit(1)
    print('PASS: no mesh problems found')


if __name__ == '__main__':
    main()

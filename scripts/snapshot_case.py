#!/usr/bin/env python3
"""Snapshot a running case into a self-contained directory safe to post-process.

Two problems this solves, both of which silently corrupt plots otherwise:

1. **Partial cycles.** The binary appends to totalop.txt/interval.txt while it
   runs, so a naive read catches a half-written block and the loaders either
   raise "does not contain a complete event block" or, worse, plot garbage.
   Every segment is truncated here to a whole number of cycles, using
   min(totalop rows // totftnode, interval lines).

2. **Restart segments.** Restarting from `binaryop` writes to
   totalop.txt{icstart} with the NEW icstart, so a run restarted at cycle 73
   leaves totalop.txt1 and totalop.txt73 side by side. plotRuptureDynamics
   reads one icstart only and would silently report just that segment. The
   segments are concatenated here, in icstart order, into a single
   totalop.txt1 / interval.txt1 covering the whole sequence, and FE_Global +
   user_defined_params are rewritten to match.

Usage:
    python3 snapshot_case.py <case_dir> [<dest_dir>]     # default <case_dir>_snap
"""

from __future__ import annotations

import argparse
import glob
import os
import re
import shutil

MESH_AND_INPUT = ['FE_Global.txt', 'FE_Model_Geometry.txt', 'FE_Fault_Geometry.txt',
                  'meshGeneralInfo.txt', 'nsmp.txt', 'vert.txt', 'fac.txt',
                  'nsmpTanLen.txt', 'nsmpGeoPhys.txt', 'user_defined_params.py']


def totftnode(case: str) -> int:
    for d in (case, os.path.join(case, 'fem_mesh_output')):
        p = os.path.join(d, 'meshGeneralInfo.txt')
        if os.path.exists(p):
            rows = [l.split() for l in open(p) if l.strip()]
            return sum(int(x) for x in rows[1])
    raise SystemExit(f'error: meshGeneralInfo.txt not found under {case}')


def find(case: str, name: str) -> str | None:
    for d in (case, os.path.join(case, 'fem_mesh_output'),
              os.path.join(case, 'aRawSimuData')):
        p = os.path.join(d, name)
        if os.path.exists(p):
            return p
    return None


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('case_dir')
    ap.add_argument('dest_dir', nargs='?', default=None)
    ap.add_argument('--quiet', action='store_true')
    args = ap.parse_args()

    case = os.path.abspath(args.case_dir.rstrip('/'))
    dest = os.path.abspath(args.dest_dir) if args.dest_dir else case + '_snap'
    n = totftnode(case)
    os.makedirs(dest, exist_ok=True)

    for f in MESH_AND_INPUT:
        p = find(case, f)
        if p:
            shutil.copy2(p, os.path.join(dest, f))

    # Collect the segments in icstart order, dropping incomplete trailing cycles.
    segs = []
    for root in (case, os.path.join(case, 'aRawSimuData')):
        for tp in glob.glob(os.path.join(root, 'totalop.txt*')):
            tag = re.search(r'totalop\.txt(\d+)$', tp)
            ip = os.path.join(root, f'interval.txt{tag.group(1)}') if tag else None
            if tag and ip and os.path.exists(ip):
                segs.append((int(tag.group(1)), tp, ip))
    segs.sort()
    if not segs:
        raise SystemExit(f'error: no totalop/interval pairs under {case}')

    out_t = os.path.join(dest, 'totalop.txt1')
    out_i = os.path.join(dest, 'interval.txt1')
    total = 0
    with open(out_t, 'w') as ft, open(out_i, 'w') as fi:
        for tag, tp, ip in segs:
            rows = sum(1 for _ in open(tp))
            ivl = [l for l in open(ip) if l.strip()]
            nb = min(rows // n, len(ivl))
            if nb == 0:
                continue
            with open(tp) as fh:
                for i, line in enumerate(fh):
                    if i >= nb * n:
                        break
                    ft.write(line)
            fi.writelines(ivl[:nb])
            total += nb
            if not args.quiet:
                print(f'  segment {tag}: {nb} complete cycles')

    # cyclelog carries (first_cycle, last_cycle), and load_saf_case() derives
    # expected_events = cyclelog[1] - cyclelog[0] + 1 from it, then truncates the
    # slip array to that many events. Writing "total total" here therefore made
    # every loader-based plot (figure3, figure4) see exactly ONE event while the
    # catalogue-based ones saw all of them -- an empty Figure 4 with no error.
    with open(os.path.join(dest, 'cyclelog.txt1'), 'w') as f:
        f.write(f' 1 {total}\n')

    # Point the metadata at the merged, single-segment sequence.
    fg = os.path.join(dest, 'FE_Global.txt')
    if os.path.exists(fg):
        lines = open(fg).read().splitlines()
        if len(lines) >= 7:
            lines[6] = f'1 {total}'
            open(fg, 'w').write('\n'.join(lines) + '\n')
    up = os.path.join(dest, 'user_defined_params.py')
    if os.path.exists(up):
        s = re.sub(r'^par\.icstart, par\.icend = .*$',
                   f'par.icstart, par.icend = 1, {total}',
                   open(up).read(), flags=re.M)
        open(up, 'w').write(s)

    if not args.quiet:
        print(f'  {total} cycles -> {dest}')
    else:
        print(f'{dest}: {total} cycles')


if __name__ == '__main__':
    main()

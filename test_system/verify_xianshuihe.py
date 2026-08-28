#!/usr/bin/env python3
"""Solver regression test against the frozen 5-cycle xianshuihe reference.

Unpacks the frozen inputs (mesh + per-node loading + FE files), runs 5 cycles
single-threaded, and compares against the stored output.

This pins the solver, not the mesher: the mesh and loading are frozen inputs,
so a gmsh version change or a GSRM re-download cannot move the answer.

It must run with OMP_NUM_THREADS=1 (R30). The OpenMP reductions in qdct3 and
hrglss accumulate under $OMP ATOMIC, so the summation order follows thread
scheduling and the result differs in the last bit between thread counts; on
this very case, 3 threads diverges from 1 thread by cycle 3 (interval 35 yr
against 38). The test sets the variable itself rather than trusting the
caller's environment.

Usage:
    python3 test_system/verify_xianshuihe.py [--keep] [--tol 0.0]
"""

from __future__ import annotations

import argparse
import gzip
import os
import shutil
import subprocess
import sys
import tarfile
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
REF = os.path.join(HERE, "reference.results", "xianshuihe.gmsh.lite")
NCYCLE = 5


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--keep", action="store_true", help="keep the scratch case dir")
    ap.add_argument("--tol", type=float, default=0.0,
                    help="absolute tolerance (default 0.0: bit-identical)")
    args = ap.parse_args()

    if not os.path.exists(os.path.join(REF, "inputs.tar.gz")):
        print(f"SKIP: no reference under {REF}")
        return 0

    version = open(os.path.join(ROOT, "VERSION")).read().strip()
    exe = os.path.join(ROOT, "bin", f"run_eqdyna2d_{version}")
    if not os.path.exists(exe):
        print(f"SKIP: {exe} not built (run 'cd src && make')")
        return 0

    case = tempfile.mkdtemp(prefix="xsh_verify_")
    try:
        with tarfile.open(os.path.join(REF, "inputs.tar.gz")) as t:
            t.extractall(case)
        shutil.copy2(exe, case)

        env = dict(os.environ, OMP_NUM_THREADS="1", GFORTRAN_UNBUFFERED_ALL="1")
        print(f"running {NCYCLE} cycles, 1 thread, {os.path.basename(exe)} ...")
        r = subprocess.run([f"./run_eqdyna2d_{version}"], cwd=case, env=env,
                           capture_output=True, text=True, timeout=3600)
        if r.returncode != 0:
            print(f"FAIL: binary exited {r.returncode}\n{r.stdout[-2000:]}")
            return 1

        got_p = os.path.join(case, "totalop.txt1")
        if not os.path.exists(got_p):
            print("FAIL: no totalop.txt1 produced")
            return 1
        got = np.loadtxt(got_p)
        with gzip.open(os.path.join(REF, "totalop.txt1.gz"), "rt") as f:
            want = np.loadtxt(f)

        if got.shape != want.shape:
            print(f"FAIL: shape {got.shape} against reference {want.shape}")
            return 1
        d = np.abs(got - want).max()
        iv_got = np.atleast_1d(np.loadtxt(os.path.join(case, "interval.txt1")))
        iv_want = np.atleast_1d(np.loadtxt(os.path.join(REF, "interval.txt1")))
        iv_ok = iv_got.shape == iv_want.shape and np.array_equal(iv_got, iv_want)

        print(f"  intervals: {' '.join(f'{v:.0f}' for v in iv_got)}")
        print(f"  reference: {' '.join(f'{v:.0f}' for v in iv_want)}")
        print(f"  totalop max abs diff: {d:.3e}")

        if d <= args.tol and iv_ok:
            print("PASS")
            return 0
        if not iv_ok:
            print("FAIL: recurrence intervals differ from the reference")
        else:
            print(f"FAIL: totalop differs by {d:.3e} (tolerance {args.tol})")
        print("If this change was meant to alter results, regenerate the "
              "reference and say why in the commit message.")
        return 1
    finally:
        if args.keep:
            print(f"case kept at {case}")
        else:
            shutil.rmtree(case, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())

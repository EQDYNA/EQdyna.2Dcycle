#!/usr/bin/env python3
"""Regression test for the source-time-function (STF) output.

Runs the first 3 cycles of the frozen xianshuihe reference with outputSTF on,
single-threaded, and checks:

  1. read-only: totalop.txt1 and interval.txt1 equal the frozen reference
     bit for bit, so turning the output on did not touch the solution;
  2. consistency: each event's final slip, shear stress and rupture time in
     stf.bin equal totalop.txt to float32 precision;
  3. kinematics: the time integral of the signed slip rate equals the signed
     slip (relative error < 1e-5);
  4. completeness: every node that ruptured (totalop rupture time < 999) is
     present in the STF record;
  5. sequence time: each event's time_yr equals the cumulative sum of
     interval.txt;
  6. moment: M0 from the STF equals the catalogue moment computed the
     plotRuptureDynamics way (relative error < 1e-4).

Usage:
    python3 test_system/verify_stf.py [--keep]
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
NCYCLE = 3
sys.path.insert(0, os.path.join(ROOT, "scripts"))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--keep", action="store_true")
    a = ap.parse_args()
    from stf_read import STFCase, MU_DEFAULT, WIDTH_DEFAULT

    version = open(os.path.join(ROOT, "VERSION")).read().strip()
    exe = os.path.join(ROOT, "bin", f"run_eqdyna2d_{version}")
    if not os.path.exists(exe):
        print(f"SKIP: {exe} not built")
        return 0

    case = tempfile.mkdtemp(prefix="stf_verify_")
    fails = []
    try:
        with tarfile.open(os.path.join(REF, "inputs.tar.gz")) as t:
            t.extractall(case)
        shutil.copy2(exe, case)
        fg = os.path.join(case, "FE_Global.txt")
        lines = open(fg).read().splitlines()
        lines[6] = f"1 {NCYCLE}"
        open(fg, "w").write("\n".join(lines) + "\n1 1 0.001\n")

        env = dict(os.environ, OMP_NUM_THREADS="1", GFORTRAN_UNBUFFERED_ALL="1")
        print(f"running {NCYCLE} cycles with STF on, 1 thread ...")
        r = subprocess.run([f"./run_eqdyna2d_{version}"], cwd=case, env=env,
                           capture_output=True, text=True, timeout=3600)
        if r.returncode != 0:
            print(f"FAIL: binary exited {r.returncode}\n{r.stdout[-1500:]}")
            return 1

        tot = np.loadtxt(os.path.join(case, "totalop.txt1"))
        iv = np.atleast_1d(np.loadtxt(os.path.join(case, "interval.txt1")))
        with gzip.open(os.path.join(REF, "totalop.txt1.gz"), "rt") as f:
            ref = np.loadtxt(f)
        ref_iv = np.atleast_1d(np.loadtxt(os.path.join(REF, "interval.txt1")))
        n = tot.shape[0]

        # 1 read-only
        if not (np.array_equal(tot, ref[:n]) and np.array_equal(iv, ref_iv[:len(iv)])):
            fails.append("totalop/interval differ from the frozen reference with STF on")
        print(f"  1 read-only: totalop max diff vs reference "
              f"{np.abs(tot - ref[:n]).max():.1e}")

        c = STFCase(case)
        N = c.files[0].nnode
        if c.eqids != list(range(1, NCYCLE + 1)):
            fails.append(f"event ids {c.eqids}, expected 1..{NCYCLE}")
        for k, e in enumerate(c.eqids):
            ev = c.event(e)
            blk = tot[k * N:(k + 1) * N]
            nd = ev["node"]
            ds = np.abs(ev["slip"][:, -1] - blk[nd, 2]).max()
            dsh = np.abs(ev["shear_stress"][:, -1] - blk[nd, 0]).max() / np.abs(blk[nd, 0]).max()
            drt = np.abs(ev["rupture_time"] - blk[nd, 4]).max()
            integ = ev["slip_rate_t"].sum(axis=1) * ev["dt"]
            smax = np.abs(ev["slip_t"][:, -1]).max()
            dk = np.abs(integ - ev["slip_t"][:, -1]).max() / smax
            missed = np.setdiff1d(np.where(blk[:, 4] < 999)[0], nd).size
            dts = abs(ev["time_yr"] - iv[:k + 1].sum())
            m0_stf = MU_DEFAULT * WIDTH_DEFAULT * float((ev["slip"][:, -1] * ev["length"]).sum())
            m0_cat = MU_DEFAULT * WIDTH_DEFAULT * float((blk[:, 2] * [c.nodes["length"]][0]).sum())
            dm = abs(m0_stf - m0_cat) / m0_cat
            print(f"  eq{e}: {ev['nout']} nodes, {ev['nt']} samples | slip {ds:.1e} m, "
                  f"shear {dsh:.1e}, rupt {drt:.1e} | int(v)dt {dk:.1e} | missed {missed} | "
                  f"time {dts:.1e} yr | M0 {dm:.1e}")
            if ds > 1e-5 * max(smax, 1.0): fails.append(f"eq{e}: final slip mismatch {ds}")
            if dsh > 1e-6: fails.append(f"eq{e}: shear mismatch {dsh}")
            if drt > 1e-6: fails.append(f"eq{e}: rupture time mismatch {drt}")
            if dk > 1e-5: fails.append(f"eq{e}: int(slip rate) != slip ({dk})")
            if missed: fails.append(f"eq{e}: {missed} ruptured node(s) missing")
            if dts > 1e-9: fails.append(f"eq{e}: time_yr off by {dts}")
            if dm > 1e-4: fails.append(f"eq{e}: moment differs by {dm}")

        if fails:
            print("FAIL:\n  " + "\n  ".join(fails))
            return 1
        print("PASS")
        return 0
    finally:
        if a.keep:
            print(f"case kept at {case}")
        else:
            shutil.rmtree(case, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""Single entry point for the test system.

    python3 test_system/run.py unit         # pytest, pure Python, seconds
    python3 test_system/run.py regression   # convention guards, one per past incident
    python3 test_system/run.py smoke        # build + mesh + one cycle (minutes)
    python3 test_system/run.py e2e          # bit-exact solver checks against frozen references
    python3 test_system/run.py ci           # unit + regression + smoke (what CI runs)
    python3 test_system/run.py all          # unit + regression + smoke + e2e (release gate)

With no tier, runs unit + regression. Several tiers may be given; they run in
order and the script exits non-zero if any fails.

e2e compares bit for bit against references produced on the project's own
machine at OMP_NUM_THREADS=1. The build uses -march=native, so on different
hardware those comparisons can fail without any regression; run e2e on the
machine that made the references, not on a CI runner.
"""

from __future__ import annotations

import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
PY = sys.executable

TIERS = {
    "unit": [[PY, "-m", "pytest", os.path.join(HERE, "unit"), "-q"]],
    "regression": [[PY, os.path.join(HERE, "test_conventions.py")]],
    "smoke": [[PY, os.path.join(HERE, "smoke.py")]],
    "e2e": [[PY, os.path.join(HERE, "verify_xianshuihe.py")],
            [PY, os.path.join(HERE, "verify_stf.py")]],
}
GROUPS = {"ci": ["unit", "regression", "smoke"],
          "all": ["unit", "regression", "smoke", "e2e"]}


def main(argv: list[str]) -> int:
    asked = argv or ["unit", "regression"]
    order = []
    for a in asked:
        for t in GROUPS.get(a, [a]):
            if t not in TIERS:
                sys.exit(f"unknown tier {a!r}; choose from {', '.join([*TIERS, *GROUPS])}")
            if t not in order:
                order.append(t)
    results = []
    for t in order:
        t0 = time.time()
        print(f"\n===== {t} =====", flush=True)
        ok = all(subprocess.call(cmd, cwd=ROOT) == 0 for cmd in TIERS[t])
        results.append((t, ok, time.time() - t0))
    print("\n===== summary =====")
    for t, ok, dt in results:
        print(f"  {t:11s} {'PASS' if ok else 'FAIL'}  {dt:6.1f} s")
    return 0 if all(ok for _, ok, _ in results) else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

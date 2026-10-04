"""snapshot_case.py: whole cycles only, restart segments merged in order."""
import os
import subprocess
import sys

import numpy as np

from conftest import ROOT

N = 5   # fault nodes: 2 + 3


def write_segment(d, tag, cycles, partial_rows=0):
    rows = []
    for c in cycles:
        rows += [f"{c}.0 0 0 0 0" for _ in range(N)]
    rows += ["9.9 9 9 9 9"] * partial_rows
    (d / f"totalop.txt{tag}").write_text("\n".join(rows) + "\n")
    (d / f"interval.txt{tag}").write_text("".join(f" {10.0 * c}\n" for c in cycles))
    (d / f"cyclelog.txt{tag}").write_text(f" {tag} {cycles[-1]}\n")


def test_merge_and_truncate(tmp_path):
    case = tmp_path / "case"
    case.mkdir()
    (case / "meshGeneralInfo.txt").write_text("2\n2 3\n")
    (case / "FE_Global.txt").write_text("\n".join(["3", "2", "4", " ", "0.01", "200", "3 4"]) + "\n")
    (case / "user_defined_params.py").write_text("par.icstart, par.icend = 3, 4\n")
    write_segment(case, 1, [1, 2], partial_rows=3)          # segment 1 ends mid-cycle
    write_segment(case, 3, [3])                             # restart segment
    dest = tmp_path / "snap"
    r = subprocess.run([sys.executable, os.path.join(ROOT, "scripts", "snapshot_case.py"),
                        str(case), str(dest), "--quiet"], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    tot = np.loadtxt(dest / "totalop.txt1")
    assert tot.shape == (3 * N, 5)
    assert list(tot[::N, 0]) == [1.0, 2.0, 3.0]             # cycles in order, partial dropped
    assert list(np.loadtxt(dest / "interval.txt1")) == [10.0, 20.0, 30.0]
    assert (dest / "cyclelog.txt1").read_text().split() == ["1", "3"]
    assert (dest / "FE_Global.txt").read_text().splitlines()[6] == "1 3"
    assert "par.icstart, par.icend = 1, 3" in (dest / "user_defined_params.py").read_text()

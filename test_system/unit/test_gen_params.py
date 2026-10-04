"""docs/user/gen_params.py: the parameter reference comes from the code."""
import subprocess
import sys
import os

import gen_params

from conftest import ROOT


def test_rows_read_defaults():
    rows = {n.strip("`"): (d, m) for n, d, m in gen_params.rows()}
    assert rows["fric_fs"][0] == "`0.5`"
    assert rows["icstart"][0] == "`1`" and rows["icend"][0] == "`200`"   # tuple assignment split
    assert rows["exe"][0] == "`run_eqdyna2d_<VERSION>`"
    assert all(m for _, m in rows.values()), "a default has no comment to document it"


def test_committed_page_is_current():
    r = subprocess.run([sys.executable, os.path.join(ROOT, "docs", "user", "gen_params.py"), "--check"],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stdout
